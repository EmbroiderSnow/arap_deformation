#include "viewer/viewer.h"
#include <iostream>

MeshViewer::MeshViewer(int width, int height) {
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
    VAO_ = 0;
    VBO_ = 0;
    EBO_ = 0;
    camera_distance_ = 1.5f;         // 减小相机距离，使模型看起来更大
    camera_theta_ = M_PI * 0.5f;    
    camera_phi_ = M_PI * 0.4f;      // 稍微调整俯视角度
    camera_pos_ = Eigen::Vector3f(0.0f, 0.0f, 5.0f);
    camera_target_ = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
    updateCamera();
    // 初始化 OpenGL 设置
    initGL();
}

void MeshViewer::handleInput() {
    if (glfwGetKey(window_, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        glfwSetWindowShouldClose(window_, true);
    }
}

void MeshViewer::renderUI() {
    // TODO: 实现 UI 渲染
}

void MeshViewer::run(Mesh& mesh) {
    mesh_ = &mesh;
    setupMesh();
    
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
        
        if (state->left_pressed) {
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
        
        if (button == GLFW_MOUSE_BUTTON_LEFT) {
            state->left_pressed = (action == GLFW_PRESS);
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
    
    // 计算法线
    trimesh.request_vertex_normals();
    if (!trimesh.has_vertex_normals()) {
        trimesh.update_normals();
    }

    // 收集顶点数据
    vertices_.clear();
    for (auto vh : trimesh.vertices()) {
        auto point = trimesh.point(vh);
        auto normal = trimesh.normal(vh);
        
        // 位置
        vertices_.push_back(point[0]);
        vertices_.push_back(point[1]);
        vertices_.push_back(point[2]);
        
        // 法线
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

    // 计算模型的边界框
    Eigen::Vector3f bbox_min = Eigen::Vector3f::Constant(std::numeric_limits<float>::max());
    Eigen::Vector3f bbox_max = Eigen::Vector3f::Constant(std::numeric_limits<float>::lowest());
    
    for (size_t i = 0; i < vertices_.size(); i += 6) {
        bbox_min.x() = std::min(bbox_min.x(), vertices_[i]);
        bbox_min.y() = std::min(bbox_min.y(), vertices_[i + 1]);
        bbox_min.z() = std::min(bbox_min.z(), vertices_[i + 2]);
        
        bbox_max.x() = std::max(bbox_max.x(), vertices_[i]);
        bbox_max.y() = std::max(bbox_max.y(), vertices_[i + 1]);
        bbox_max.z() = std::max(bbox_max.z(), vertices_[i + 2]);
    }
    
    // 计算模型中心和大小
    Eigen::Vector3f center = (bbox_max + bbox_min) * 0.5f;
    float max_extent = (bbox_max - bbox_min).maxCoeff();
    
    // 归一化模型大小并居中
    for (size_t i = 0; i < vertices_.size(); i += 6) {
        vertices_[i] = (vertices_[i] - center.x()) / max_extent;
        vertices_[i + 1] = (vertices_[i + 1] - center.y()) / max_extent;
        vertices_[i + 2] = (vertices_[i + 2] - center.z()) / max_extent;
    }
    
    // 创建并绑定VAO和VBO
    glGenVertexArrays(1, &VAO_);
    glGenBuffers(1, &VBO_);
    glGenBuffers(1, &EBO_);
    
    glBindVertexArray(VAO_);
    
    glBindBuffer(GL_ARRAY_BUFFER, VBO_);
    glBufferData(GL_ARRAY_BUFFER, vertices_.size() * sizeof(float), 
                 vertices_.data(), GL_STATIC_DRAW);
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO_);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices_.size() * sizeof(unsigned int),
                 indices_.data(), GL_STATIC_DRAW);
    
    // 设置顶点属性指针
    // 位置属性
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    // 法线属性
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

MeshViewer::~MeshViewer() {
    // 清理资源
    if (shader_) {
        delete shader_;
        shader_ = nullptr;
    }
    
    if (VAO_) {
        glDeleteVertexArrays(1, &VAO_);
    }
    if (VBO_) {
        glDeleteBuffers(1, &VBO_);
    }
    if (EBO_) {
        glDeleteBuffers(1, &EBO_);
    }

    if (window_) {
        MouseState* state = static_cast<MouseState*>(glfwGetWindowUserPointer(window_));
        delete state;
        glfwSetWindowUserPointer(window_, nullptr);
        
        glfwDestroyWindow(window_);
        window_ = nullptr;
    }
    glfwTerminate();
}