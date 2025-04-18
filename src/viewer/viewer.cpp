#include "viewer/viewer.h"
#include <iostream>

MeshViewer::MeshViewer(int width, int height)
    : width_(width), height_(height) {
    // 初始化 GLFW
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return;
    }

    // 设置 OpenGL 版本和核心模式
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    // 创建窗口
    window_ = glfwCreateWindow(width, height, "ARAP Deformation", NULL, NULL);
    if (!window_) {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return;
    }

    glfwMakeContextCurrent(window_);

    // 初始化 GLEW
    glewExperimental = GL_TRUE; 
    GLenum err = glewInit();
    if (err != GLEW_OK) {
        std::cerr << "GLEW initialization failed: " 
                  << glewGetErrorString(err) << std::endl;
        glfwDestroyWindow(window_);
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    // 初始化成员变量
    shader_ = nullptr;
    mesh_ = nullptr;
    solver_ = nullptr;
    VAO_ = 0;
    VBO_ = 0;
    EBO_ = 0;
    selected_VAO_ = 0;
    selected_VBO_ = 0;
    camera_distance_ = 1.25f;         // 减小相机距离，使模型看起来更大
    camera_theta_ = M_PI * 0.5f;    
    camera_phi_ = M_PI * 0.4f;      // 稍微调整俯视角度
    camera_pos_ = Eigen::Vector3f(0.0f, 0.0f, 5.0f);
    camera_target_ = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
    updateCamera();
    // 初始化 OpenGL 设置
    initGL();
}

void MeshViewer::handleInput() {
    static bool s_pressed_last = false;
    if (glfwGetKey(window_, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        glfwSetWindowShouldClose(window_, true);
    }

    // 切换选择模式
    if (glfwGetKey(window_, GLFW_KEY_S) == GLFW_PRESS) {
        if (deform_mode_) {
            deform_mode_ = false;
            std::cout << "Deformation mode deactivated" << std::endl;
        }
        selection_mode_ = true;
        std::cout << "Selection mode activated" << std::endl;
    } 

    if (glfwGetKey(window_, GLFW_KEY_D) == GLFW_PRESS) {
        if (selection_mode_) {
            selection_mode_ = false;
            std::cout << "Selection mode deactivated" << std::endl;
        }
        deform_mode_ = true;
        std::cout << "Deformation mode activated" << std::endl;
    }

    if (glfwGetKey(window_, GLFW_KEY_Q) == GLFW_PRESS) {
        selection_mode_ = false;
        deform_mode_ = false;
        std::cout << "All modes deactivated" << std::endl;
    }
    
    // 清除选择
    if (glfwGetKey(window_, GLFW_KEY_C) == GLFW_PRESS) {
        selected_vertices_.clear();
        std::cout << "Selection cleared" << std::endl;
    }

    // 按F键让模型回到视野中央
    if (glfwGetKey(window_, GLFW_KEY_F) == GLFW_PRESS) {
        camera_target_ = center_;
        updateCamera();
        std::cout << "Model recentered" << std::endl;
    }

    // 按S键保存当前网格为obj
    bool s_pressed = glfwGetKey(window_, GLFW_KEY_E) == GLFW_PRESS;
    if (s_pressed && !s_pressed_last) {
        if (mesh_) {
            if (mesh_->save("output.obj")) {
                std::cout << "Mesh saved to output.obj" << std::endl;
            } else {
                std::cerr << "Failed to save mesh!" << std::endl;
            }
        }
    }
    s_pressed_last = s_pressed;
}

void MeshViewer::renderUI() {
    // TODO: 实现 UI 渲染
}

void MeshViewer::run(Mesh& mesh) {
    mesh_ = &mesh;
    initializeMesh();

    solver_ = new ARAPSolver(*mesh_);
    solver_->setWeight(1.0, 1000.0); // 设置权重
    
    while (!glfwWindowShouldClose(window_)) {
        handleInput();
        renderMesh();
        renderUI();
        
        glfwSwapBuffers(window_);
        glfwPollEvents();
    }
}

namespace {
    // 鼠标回调函数的用户数据结构
    struct MouseState {
        bool left_pressed = false;
        bool right_pressed = false;
        float last_x = 0.0f;
        float last_y = 0.0f;
        MeshViewer* viewer = nullptr;
    };
}

void MeshViewer::initGL() {
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
    
    try {
        // 添加路径检查
        std::string vert_path = "./shaders/mesh.vert";  // 使用显式的相对路径
        std::string frag_path = "./shaders/mesh.frag";
        
        // 检查文件是否存在和可读
        std::ifstream vert_test(vert_path);
        std::ifstream frag_test(frag_path);
        
        if (!vert_test.is_open()) {
            throw std::runtime_error("Cannot open vertex shader: " + vert_path);
        }
        if (!frag_test.is_open()) {
            throw std::runtime_error("Cannot open fragment shader: " + frag_path);
        }
        
        shader_ = new Shader(vert_path.c_str(), frag_path.c_str());
    }
    catch (const std::exception& e) {
        std::cerr << "Shader initialization failed: " << e.what() << std::endl;
        exit(EXIT_FAILURE);  // 如果着色器加载失败，直接退出
    }
    
    // 初始化相机
    camera_pos_ = Eigen::Vector3f(0.0f, 0.0f, 5.0f);
    camera_target_ = Eigen::Vector3f(0.0f, 0.0f, 0.0f);

    // 初始化鼠标状态
    MouseState* state = new MouseState();
    state->viewer = this;
    glfwSetWindowUserPointer(window_, state);
    
    // 设置鼠标回调
    glfwSetCursorPosCallback(window_, [](GLFWwindow* window, double xpos, double ypos) {
        MouseState* state = static_cast<MouseState*>(glfwGetWindowUserPointer(window));
        if (!state->viewer) return;

        // 只有在非变形模式下允许相机旋转
        if (state->left_pressed && !state->viewer->deform_mode_) {
            float dx = xpos - state->last_x;
            float dy = ypos - state->last_y;
            state->viewer->camera_theta_ += dx * 0.01;
            state->viewer->camera_phi_ = std::max(0.1, 
                std::min(M_PI - 0.1, state->viewer->camera_phi_ + dy * 0.01));
        }

        if (state->right_pressed) {
            float dy = ypos - state->last_y;
            state->viewer->camera_distance_ = std::max(0.1f, 
                state->viewer->camera_distance_ + dy * 0.1f);
        }

        state->last_x = xpos;
        state->last_y = ypos;
    });
    
    glfwSetMouseButtonCallback(window_, [](GLFWwindow* window, int button, int action, int mods) {
        MouseState* state = static_cast<MouseState*>(glfwGetWindowUserPointer(window));
        if (!state->viewer) return;
        
        if (button == GLFW_MOUSE_BUTTON_LEFT) {
            state->left_pressed = (action == GLFW_PRESS);
            
            double xpos, ypos;
            glfwGetCursorPos(window, &xpos, &ypos);

            if (state->viewer->selection_mode_ && action == GLFW_PRESS) {
                state->viewer->handleSelection(xpos, ypos);
            }
            else if (state->viewer->deform_mode_) {
                if (action == GLFW_PRESS) {
                    // 仅记录拖动起点，不做变形
                    state->viewer->is_dragging_ = true;
                    state->viewer->drag_start_x_ = xpos;
                    state->viewer->drag_start_y_ = ypos;
                    auto vh = state->viewer->selected_vertices_.front();
                    auto point = state->viewer->mesh_->getMesh().point(vh);
                    state->viewer->drag_start_world_ = Eigen::Vector3f(point[0], point[1], point[2]);
                    state->viewer->drag_initial_positions_.clear();
                    for (const auto& vh : state->viewer->selected_vertices_) {
                        auto p = state->viewer->mesh_->getMesh().point(vh);
                        state->viewer->drag_initial_positions_[vh.idx()] = Eigen::Vector3d(p[0], p[1], p[2]);
                    }
                } else if (action == GLFW_RELEASE && state->viewer->is_dragging_) {
                    // 只在松开时做一次变形
                    state->viewer->is_dragging_ = false;
                    state->viewer->handleDeformation(xpos, ypos);
                }
            }
        }
        if (button == GLFW_MOUSE_BUTTON_RIGHT) {
            state->right_pressed = (action == GLFW_PRESS);
        }
        
        if (action == GLFW_PRESS) {
            double xpos, ypos;
            glfwGetCursorPos(window, &xpos, &ypos);
            state->last_x = xpos;
            state->last_y = ypos;
        }
    });
}

void MeshViewer::setupMesh() {
    if (!mesh_) {
        std::cerr << "No mesh loaded" << std::endl;
        return;
    }

    auto& trimesh = mesh_->getMesh();
    trimesh.request_vertex_normals();
    if (!trimesh.has_vertex_normals()) {
        trimesh.update_normals();
    }

    vertices_.clear();
    for (auto vh : trimesh.vertices()) {
        auto point = trimesh.point(vh);
        if (!std::isfinite(point[0]) || !std::isfinite(point[1]) || !std::isfinite(point[2])) {
            std::cerr << "Vertex " << vh.idx() << " is NaN or Inf!" << std::endl;
        }
        auto normal = trimesh.normal(vh);
        float x = (point[0] - center_.x()) / max_extent_;
        float y = (point[1] - center_.y()) / max_extent_;
        float z = (point[2] - center_.z()) / max_extent_;
        vertices_.push_back(x);
        vertices_.push_back(y);
        vertices_.push_back(z);
        vertices_.push_back(normal[0]);
        vertices_.push_back(normal[1]);
        vertices_.push_back(normal[2]);
    }
    
    // 收集索引数据
    indices_.clear();
    for (auto fh : trimesh.faces()) {
        for (auto fv_it = trimesh.cfv_begin(fh); fv_it.is_valid(); ++fv_it) {
            indices_.push_back(fv_it->idx());
        }
    }

    // 更新 VBO 数据
    if (!VAO_) glGenVertexArrays(1, &VAO_);
    if (!VBO_) glGenBuffers(1, &VBO_);
    if (!EBO_) glGenBuffers(1, &EBO_);
    
    glBindVertexArray(VAO_);
    
    glBindBuffer(GL_ARRAY_BUFFER, VBO_);
    glBufferData(GL_ARRAY_BUFFER, vertices_.size() * sizeof(float), 
                 vertices_.data(), GL_DYNAMIC_DRAW); // 改为 DYNAMIC_DRAW
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO_);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices_.size() * sizeof(unsigned int),
                 indices_.data(), GL_STATIC_DRAW);
    
    // 设置顶点属性指针
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), 
                         (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
}

void MeshViewer::renderMesh() {
    int width, height;
    glfwGetFramebufferSize(window_, &width, &height);
    glViewport(0, 0, width, height);
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    shader_->use();
    updateMatrices();
    
    // 设置光照参数
    shader_->setVec3("lightPos", 3.0f, 5.0f, 3.0f);      // 调整光源位置
    shader_->setVec3("viewPos", camera_pos_.x(), camera_pos_.y(), camera_pos_.z());
    shader_->setVec3("lightColor", 1.0f, 1.0f, 1.0f);    // 白色光源
    shader_->setVec3("objectColor", 0.8f, 0.8f, 1.0f);   // 更亮的物体颜色
    
    glBindVertexArray(VAO_);
    glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, 0);

    // 先绘制实体
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    shader_->setVec3("objectColor", 0.8f, 0.8f, 1.0f);   // 物体颜色
    glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, 0);
    
    // 再绘制线框
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glEnable(GL_POLYGON_OFFSET_LINE);
    glPolygonOffset(-1.0f, -1.0f);  // 避免z-fighting
    
    shader_->setVec3("objectColor", 0.2f, 0.2f, 0.2f);   // 线框颜色
    glLineWidth(1.0f);  // 设置线宽
    glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, 0);
    
    // 恢复默认设置
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDisable(GL_POLYGON_OFFSET_LINE);

    // 渲染选中的点
    if (!selected_vertices_.empty() && selected_VAO_) {
        updateSelectedVertices();
        
        glEnable(GL_PROGRAM_POINT_SIZE);
        shader_->use();
        shader_->setVec3("objectColor", 1.0f, 0.0f, 0.0f);  // 红色
        shader_->setVec3("lightColor", 1.0f, 1.0f, 1.0f);   // 不受光照影响
        
        glBindVertexArray(selected_VAO_);
        glPointSize(10.0f);
        glDrawArrays(GL_POINTS, 0, selected_vertices_.size());
        
        glDisable(GL_PROGRAM_POINT_SIZE);
    }
    GLenum err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cout << "selcted_VAO_:" << selected_VAO_ << std::endl;
        std::cout << "selcted_VBO_:" << selected_VBO_ << std::endl;
        std::cerr << "OpenGL error(line356): " << err << std::endl;
    }
}

void MeshViewer::initializeMesh() {
    if (!mesh_) return;
    
    // 计算边界框
    Eigen::Vector3f bbox_min = Eigen::Vector3f::Constant(std::numeric_limits<float>::max());
    Eigen::Vector3f bbox_max = Eigen::Vector3f::Constant(std::numeric_limits<float>::lowest());
    
    auto& trimesh = mesh_->getMesh();
    // 首次计算边界框
    for (auto vh : trimesh.vertices()) {
        auto point = trimesh.point(vh);
        bbox_min.x() = std::min(bbox_min.x(), static_cast<float>(point[0]));
        bbox_min.y() = std::min(bbox_min.y(), static_cast<float>(point[1]));
        bbox_min.z() = std::min(bbox_min.z(), static_cast<float>(point[2]));
        
        bbox_max.x() = std::max(bbox_max.x(), static_cast<float>(point[0]));
        bbox_max.y() = std::max(bbox_max.y(), static_cast<float>(point[1]));
        bbox_max.z() = std::max(bbox_max.z(), static_cast<float>(point[2]));
    }
    
    // 计算模型中心和大小
    center_ = (bbox_max + bbox_min) * 0.5f;
    max_extent_ = (bbox_max - bbox_min).maxCoeff();
    
    // 在 initializeMesh() 里保存原始坐标
    original_max_extent_ = max_extent_;
    original_vertex_positions_.clear();
    for (auto vh : trimesh.vertices()) {
        auto point = trimesh.point(vh);
        original_vertex_positions_.emplace_back(point[0], point[1], point[2]);
    }
    
    is_initialized_ = true;
    setupMesh();

    if (max_extent_ < 1e-8) {
        std::cerr << "Warning: max_extent_ is too small (" << max_extent_ << "), model may disappear!" << std::endl;
    }
}

void MeshViewer::updateCamera() {
    // 根据球面坐标计算相机位置
    float x = camera_distance_ * sin(camera_phi_) * cos(camera_theta_);
    float y = camera_distance_ * cos(camera_phi_);
    float z = camera_distance_ * sin(camera_phi_) * sin(camera_theta_);
    
    camera_pos_ = Eigen::Vector3f(x, y, z) + camera_target_;
}

Eigen::Matrix4f MeshViewer::getViewMatrix() {
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();
    
    Eigen::Vector3f up(0.0f, 1.0f, 0.0f);
    Eigen::Vector3f f = (camera_target_ - camera_pos_).normalized();
    Eigen::Vector3f s = f.cross(up).normalized();
    Eigen::Vector3f u = s.cross(f);
    
    view.block<3,1>(0,0) = s;
    view.block<3,1>(0,1) = u;
    view.block<3,1>(0,2) = -f;
    view(0,3) = -s.dot(camera_pos_);
    view(1,3) = -u.dot(camera_pos_);
    view(2,3) = f.dot(camera_pos_);
    
    return view;
}

Eigen::Matrix4f MeshViewer::getProjectionMatrix() {
    int width, height;
    glfwGetFramebufferSize(window_, &width, &height);
    float aspect = static_cast<float>(width) / static_cast<float>(height);

    float fov = 45.0f;  // 视场角
    float near = 0.1f;  // 近平面
    float far = 100.0f;  // 远平面
    
    Eigen::Matrix4f projection = Eigen::Matrix4f::Zero();
    float tanHalfFovy = tan(fov / 2.0f * M_PI / 180.0f);
    
    projection(0,0) = 1.0f / (aspect * tanHalfFovy);
    projection(1,1) = 1.0f / tanHalfFovy;
    projection(2,2) = -(far + near) / (far - near);
    projection(2,3) = -2.0f * far * near / (far - near);
    projection(3,2) = -1.0f;
    
    return projection;
}

void MeshViewer::updateMatrices() {
    updateCamera();
    
    Eigen::Matrix4f model = Eigen::Matrix4f::Identity();
    Eigen::Matrix4f view = getViewMatrix();
    Eigen::Matrix4f projection = getProjectionMatrix();
    
    shader_->use();
    shader_->setMat4("model", model.data());
    shader_->setMat4("view", view.data());
    shader_->setMat4("projection", projection.data());
}

void MeshViewer::handleSelection(double xpos, double ypos) {
    if (!selection_mode_) return;
    last_click_x_ = xpos;
    last_click_y_ = ypos;

    int width, height;
    glfwGetFramebufferSize(window_, &width, &height);

    auto& trimesh = mesh_->getMesh();
    float min_dist = 15.0f; // 最大像素距离
    Mesh::VertexHandle nearest_vh;
    Eigen::Matrix4f model = Eigen::Matrix4f::Identity();
    Eigen::Matrix4f view = getViewMatrix();
    Eigen::Matrix4f proj = getProjectionMatrix();

    for (auto vh : trimesh.vertices()) {
        auto p = trimesh.point(vh);
        Eigen::Vector4f pos((p[0] - center_.x()) / max_extent_,
                            (p[1] - center_.y()) / max_extent_,
                            (p[2] - center_.z()) / max_extent_, 1.0f);
        Eigen::Vector4f clip = proj * view * model * pos;
        clip /= clip.w();

        float sx = (clip.x() * 0.5f + 0.5f) * width;
        float sy = (0.5f - clip.y() * 0.5f) * height;

        float dist = std::sqrt((sx - xpos) * (sx - xpos) + (sy - ypos) * (sy - ypos));
        if (dist < min_dist) {
            min_dist = dist;
            nearest_vh = vh;
        }
    }

    if (nearest_vh.is_valid()) {
        auto it = std::find(selected_vertices_.begin(), selected_vertices_.end(), nearest_vh);
        if (it == selected_vertices_.end()) {
            selected_vertices_.push_back(nearest_vh);
            std::cout << "Selected vertex " << nearest_vh.idx() << std::endl;
        } else {
            selected_vertices_.erase(it);
            std::cout << "Deselected vertex " << nearest_vh.idx() << std::endl;
        }
        updateSelectedVertices();
    }
}

Mesh::VertexHandle MeshViewer::findNearestVertex(const Eigen::Vector3f& origin, 
    const Eigen::Vector3f& direction) {
        if (!mesh_) return Mesh::VertexHandle();
    
        auto& trimesh = mesh_->getMesh();
        float min_t = std::numeric_limits<float>::max();  // 记录最近交点的t值
    Mesh::VertexHandle nearest_vertex;
    bool found_intersection = false;
    
    // 首先找到最近的交点
    for (auto fh : trimesh.faces()) {
        auto fv_it = trimesh.cfv_begin(fh);
        auto v0 = *fv_it; ++fv_it;
        auto v1 = *fv_it; ++fv_it;
        auto v2 = *fv_it;
        
        auto p0 = trimesh.point(v0);
        auto p1 = trimesh.point(v1);
        auto p2 = trimesh.point(v2);
        
        // 转换到归一化空间
        Eigen::Vector3f v0_pos((p0[0] - center_.x()) / max_extent_,
                              (p0[1] - center_.y()) / max_extent_,
                              (p0[2] - center_.z()) / max_extent_);
        Eigen::Vector3f v1_pos((p1[0] - center_.x()) / max_extent_,
                              (p1[1] - center_.y()) / max_extent_,
                              (p1[2] - center_.z()) / max_extent_);
        Eigen::Vector3f v2_pos((p2[0] - center_.x()) / max_extent_,
                              (p2[1] - center_.y()) / max_extent_,
                              (p2[2] - center_.z()) / max_extent_);
        
        float t;
        if (rayTriangleIntersect(origin, direction, v0_pos, v1_pos, v2_pos, t)) {
            // 只处理最近的交点
            if (t < min_t) {
                min_t = t;
                Eigen::Vector3f intersection = origin + direction * t;
                
                // 在这个三角形的顶点中找最近的
                float min_vertex_dist = std::numeric_limits<float>::max();
                
                std::array<std::pair<float, Mesh::VertexHandle>, 3> distances = {
                    {std::make_pair((v0_pos - intersection).squaredNorm(), v0),
                     std::make_pair((v1_pos - intersection).squaredNorm(), v1),
                     std::make_pair((v2_pos - intersection).squaredNorm(), v2)}
                };
                
                auto min_dist = std::min_element(distances.begin(), distances.end());
                nearest_vertex = min_dist->second;
                found_intersection = true;
            }
        }
    }
    
    if (found_intersection) {
        std::cout << "Found intersection at t = " << min_t 
                  << ", vertex: " << nearest_vertex.idx() << std::endl;
        return nearest_vertex;
    }
    
    return Mesh::VertexHandle();
}

bool MeshViewer::rayTriangleIntersect(const Eigen::Vector3f& orig, 
const Eigen::Vector3f& dir,
const Eigen::Vector3f& v0, 
const Eigen::Vector3f& v1, 
const Eigen::Vector3f& v2,
float& t) {
    // Möller–Trumbore intersection algorithm
    const float EPSILON = 1e-7;  // 减小 EPSILON
    Eigen::Vector3f edge1 = v1 - v0;
    Eigen::Vector3f edge2 = v2 - v0;
    Eigen::Vector3f pvec = dir.cross(edge2);
    float det = edge1.dot(pvec);
    
    // 先检查行列式，确保射线与三角形不平行
    if (std::abs(det) < EPSILON) {
        return false;
    }
    
    float inv_det = 1.0f / det;
    Eigen::Vector3f tvec = orig - v0;
    float u = tvec.dot(pvec) * inv_det;
    
    if (u < -EPSILON || u > 1.0f + EPSILON) return false;
    
    Eigen::Vector3f qvec = tvec.cross(edge1);
    float v = dir.dot(qvec) * inv_det;
    
    if (v < -EPSILON || u + v > 1.0f + EPSILON) return false;
    
    t = edge2.dot(qvec) * inv_det;
    
    // 计算三角形的法线
    Eigen::Vector3f normal = edge1.cross(edge2).normalized();
    float view_dot_normal = dir.dot(normal);
    
    // 调试输出
    static int count = 0;
    if (count++ % 1000 == 0) {
        std::cout << "\nRay-triangle intersection test:" << std::endl;
        std::cout << "det: " << det << std::endl;
        std::cout << "view_dot_normal: " << view_dot_normal << std::endl;
        std::cout << "t: " << t << std::endl;
        std::cout << "u: " << u << ", v: " << v << std::endl;
    }
    
    return t > EPSILON;
}

void MeshViewer::updateSelectedVertices() {
    std::vector<float> selected_points;
    const auto& trimesh = mesh_->getMesh();
    
    for (const auto& vh : selected_vertices_) {
        if (!trimesh.is_valid_handle(vh)) {
            std::cerr << "Invalid vertex handle in selected_vertices_!" << std::endl;
            continue;
        }
        auto point = trimesh.point(vh);
        // 转换到归一化空间
        float x = (point[0] - center_.x()) / max_extent_;
        float y = (point[1] - center_.y()) / max_extent_;
        float z = (point[2] - center_.z()) / max_extent_;
        
        selected_points.push_back(x);
        selected_points.push_back(y);
        selected_points.push_back(z);
    }

    if (selected_points.empty()) {
        glBindVertexArray(0);
        return;
    }
    
    if (!selected_VAO_) glGenVertexArrays(1, &selected_VAO_);
    if (!selected_VBO_) glGenBuffers(1, &selected_VBO_);
    
    
    glBindVertexArray(selected_VAO_);
    glBindBuffer(GL_ARRAY_BUFFER, selected_VBO_);
    glBufferData(GL_ARRAY_BUFFER, selected_points.size() * sizeof(float),
                 selected_points.data(), GL_DYNAMIC_DRAW);
    
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glDisableVertexAttribArray(1);

    GLenum err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cout << "selcted_VAO_:" << selected_VAO_ << std::endl;
        std::cout << "selcted_VBO_:" << selected_VBO_ << std::endl;
        std::cerr << "OpenGL error(line658): " << err << std::endl;
    }
    if (selected_points.empty() || selected_points.size() % 3 != 0)
        glBindVertexArray(0);
}

MeshViewer::~MeshViewer() {
    // 清理资源
    if (shader_) {
        delete shader_;
        shader_ = nullptr;
    }
    
    if (VAO_) {
        glDeleteVertexArrays(1, &VAO_);
        VAO_ = 0;
    }
    if (VBO_) {
        glDeleteBuffers(1, &VBO_);
        VBO_ = 0;
    }
    if (EBO_) {
        glDeleteBuffers(1, &EBO_);
        EBO_ = 0;
    }

    if (window_) {
        MouseState* state = static_cast<MouseState*>(glfwGetWindowUserPointer(window_));
        delete state;
        glfwSetWindowUserPointer(window_, nullptr);
        
        glfwDestroyWindow(window_);
        window_ = nullptr;
    }
    if (selected_VAO_) {
        glDeleteVertexArrays(1, &selected_VAO_);
        selected_VAO_ = 0;
    }
    if (selected_VBO_) {
        glDeleteBuffers(1, &selected_VBO_);
        selected_VBO_ = 0;
    }
    glfwTerminate();
}

Eigen::Vector3f MeshViewer::screenToWorld(double xpos, double ypos, float depth) {
    // 获取视口大小
    int width, height;
    glfwGetFramebufferSize(window_, &width, &height);
    
    // 转换到NDC空间
    float x = (2.0f * xpos) / width - 1.0f;
    float y = 1.0f - (2.0f * ypos) / height;
    float z = 2.0f * depth - 1.0f;
    
    // NDC到世界空间
    Eigen::Vector4f ndc(x, y, z, 1.0f);
    Eigen::Matrix4f view_proj_inv = (getProjectionMatrix() * getViewMatrix()).inverse();
    Eigen::Vector4f world = view_proj_inv * ndc;
    return world.head<3>() / world.w();
}

void MeshViewer::updateMeshBuffers() {
    // 重新计算模型中心和最大包围盒
    if (!mesh_) return;

    Eigen::Vector3f bbox_min = Eigen::Vector3f::Constant(std::numeric_limits<float>::max());
    Eigen::Vector3f bbox_max = Eigen::Vector3f::Constant(std::numeric_limits<float>::lowest());
    auto& trimesh = mesh_->getMesh();
    for (auto vh : trimesh.vertices()) {
        auto point = trimesh.point(vh);
        bbox_min.x() = std::min(bbox_min.x(), static_cast<float>(point[0]));
        bbox_min.y() = std::min(bbox_min.y(), static_cast<float>(point[1]));
        bbox_min.z() = std::min(bbox_min.z(), static_cast<float>(point[2]));
        bbox_max.x() = std::max(bbox_max.x(), static_cast<float>(point[0]));
        bbox_max.y() = std::max(bbox_max.y(), static_cast<float>(point[1]));
        bbox_max.z() = std::max(bbox_max.z(), static_cast<float>(point[2]));
    }
    center_ = (bbox_max + bbox_min) * 0.5f;
    max_extent_ = (bbox_max - bbox_min).maxCoeff();

    // 让相机目标点始终对准模型中心
    float max_allowed_extent = original_max_extent_ * 10.0f;
    float used_extent = std::min(max_extent_, max_allowed_extent);
    camera_target_ = center_;
    camera_distance_ = used_extent * 1.25f / std::tan(45.0f * 0.5f * M_PI / 180.0f);
    updateCamera();
    setupMesh();
}

void MeshViewer::handleDeformation(double xpos, double ypos) {
    if (!solver_ || selected_vertices_.empty()) return;

    // 保存变形前的center、max_extent和相机参数
    Eigen::Vector3f center_before = center_;
    float max_extent_before = max_extent_;
    Eigen::Vector3f camera_pos_before = camera_pos_;
    Eigen::Vector3f camera_target_before = camera_target_;
    float camera_distance_before = camera_distance_;
    float camera_theta_before = camera_theta_;
    float camera_phi_before = camera_phi_;

    // 用变形前的相机参数计算view/proj
    auto getViewMatrixBefore = [&]() {
        Eigen::Matrix4f view = Eigen::Matrix4f::Identity();
        Eigen::Vector3f up(0.0f, 1.0f, 0.0f);
        Eigen::Vector3f f = (camera_target_before - camera_pos_before).normalized();
        Eigen::Vector3f s = f.cross(up).normalized();
        Eigen::Vector3f u = s.cross(f);
        view.block<3,1>(0,0) = s;
        view.block<3,1>(0,1) = u;
        view.block<3,1>(0,2) = -f;
        view(0,3) = -s.dot(camera_pos_before);
        view(1,3) = -u.dot(camera_pos_before);
        view(2,3) = f.dot(camera_pos_before);
        return view;
    };
    auto getProjectionMatrixBefore = [&]() {
        int width, height;
        glfwGetFramebufferSize(window_, &width, &height);
        float aspect = static_cast<float>(width) / static_cast<float>(height);
        float fov = 45.0f;
        float near = 0.1f;
        float far = 100.0f;
        Eigen::Matrix4f projection = Eigen::Matrix4f::Zero();
        float tanHalfFovy = tan(fov / 2.0f * M_PI / 180.0f);
        projection(0,0) = 1.0f / (aspect * tanHalfFovy);
        projection(1,1) = 1.0f / tanHalfFovy;
        projection(2,2) = -(far + near) / (far - near);
        projection(2,3) = -2.0f * far * near / (far - near);
        projection(3,2) = -1.0f;
        return projection;
    };

    Eigen::Matrix4f model = Eigen::Matrix4f::Identity();
    Eigen::Matrix4f view = getViewMatrixBefore();
    Eigen::Matrix4f proj = getProjectionMatrixBefore();
    Eigen::Vector4f pos((drag_start_world_.x() - center_before.x()) / max_extent_before,
                        (drag_start_world_.y() - center_before.y()) / max_extent_before,
                        (drag_start_world_.z() - center_before.z()) / max_extent_before, 1.0f);
    Eigen::Vector4f clip = proj * view * model * pos;
    clip /= clip.w();
    float start_depth = (clip.z() * 0.5f + 0.5f);

    // 当前鼠标位置反投影到世界坐标（用变形前的参数和相机参数）
    // 这里也要用变形前的view/proj
    auto screenToWorldBefore = [&](double xpos, double ypos, float depth) {
        int width, height;
        glfwGetFramebufferSize(window_, &width, &height);
        float x = (2.0f * xpos) / width - 1.0f;
        float y = 1.0f - (2.0f * ypos) / height;
        float z = 2.0f * depth - 1.0f;
        Eigen::Vector4f ndc(x, y, z, 1.0f);

        Eigen::Matrix4f proj_inv = getProjectionMatrixBefore().inverse();
        Eigen::Matrix4f view_inv = getViewMatrixBefore().inverse();
        Eigen::Vector4f view_space = proj_inv * ndc;
        view_space /= view_space.w();
        Eigen::Vector4f world_normalized = view_inv * view_space;
        world_normalized /= world_normalized.w();

        // 反归一化
        Eigen::Vector3f world = world_normalized.head<3>() * max_extent_before + center_before;
        return world;
    };

    // 反投影拖动起点和终点的屏幕坐标到世界坐标
    Eigen::Vector3f world0 = screenToWorldBefore(drag_start_x_, drag_start_y_, start_depth);
    Eigen::Vector3f world1 = screenToWorldBefore(xpos, ypos, start_depth);
    Eigen::Vector3d delta = (world1 - world0).cast<double>();
    std::cout << "world0: " << world0.transpose() << std::endl;
    std::cout << "world1: " << world1.transpose() << std::endl;
    std::cout << "delta: " << delta.transpose() << std::endl;
    for (const auto& vh : selected_vertices_) {
        Eigen::Vector3d new_pos = drag_initial_positions_[vh.idx()] + delta;
        mesh_->setHandleConstraint(vh, new_pos);
    }

    auto& trimesh = mesh_->getMesh();
    std::vector<int> anchor_indices = {0, int(trimesh.n_vertices())/3, int(trimesh.n_vertices())-1};
    for (int idx : anchor_indices) {
        if (idx >= 0 && idx < original_vertex_positions_.size()) {
            auto vh = Mesh::VertexHandle(idx);
            mesh_->setHandleConstraint(vh, original_vertex_positions_[idx]);
        }
    }

    if (solver_->solve(10)) {
        std::cout << "Deformation step completed" << std::endl;
        std::cout << "max_extent_: " << max_extent_ << std::endl;
        updateMeshBuffers();
    }
}