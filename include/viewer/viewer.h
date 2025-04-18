#pragma once
#include "mesh/mesh.h"
#include "deformation/arap_solver.h"
#include "shader.h"
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <vector>
#include <optional>

class MeshViewer {
public:
    MeshViewer(int width = 1280, int height = 720);
    ~MeshViewer();

    void run(Mesh& mesh);

private:
    int width_ = 1280;
    int height_ = 720;
    void initGL();
    void renderMesh();
    void initializeMesh();
    void handleInput();
    void renderUI();

    GLFWwindow* window_;
    std::vector<Mesh::VertexHandle> selected_vertices_;
    std::optional<Mesh::VertexHandle> nearest_vertex_;
    bool selection_mode_ = false;
    bool dragging_ = false;

    // OpenGL 渲染相关
    GLuint VAO_, VBO_, EBO_;
    Shader* shader_;
    
    // 相机参数
    Eigen::Vector3f camera_pos_;
    Eigen::Vector3f camera_target_;
    float camera_distance_ = 5.0f;
    float camera_theta_ = 0.0f;
    float camera_phi_ = 0.0f;
    
    // 渲染数据
    std::vector<float> vertices_;
    std::vector<unsigned int> indices_;
    
    // 网格引用
    Mesh* mesh_;
    
    // 辅助函数
    void setupMesh();
    void updateCamera();
    Eigen::Matrix4f getViewMatrix();
    Eigen::Matrix4f getProjectionMatrix();
    void updateMatrices();

    // 添加选择相关函数
    void handleSelection(double xpos, double ypos);
    Mesh::VertexHandle findNearestVertex(const Eigen::Vector3f& origin, const Eigen::Vector3f& direction);
    bool rayTriangleIntersect(const Eigen::Vector3f& orig, 
                             const Eigen::Vector3f& dir,
                             const Eigen::Vector3f& v0,
                             const Eigen::Vector3f& v1,
                             const Eigen::Vector3f& v2,
                             float& t);

    Eigen::Vector3f center_;      // 模型中心点
    float max_extent_ = 1.0f;     // 模型最大尺寸
    float original_max_extent_ = 1.0f; // 原始模型最大尺寸
    GLuint selected_VAO_, selected_VBO_;  // 用于渲染选中点的缓冲区
    void updateSelectedVertices();         // 更新选中点的缓冲区

    // 变形相关
    bool deform_mode_ = false;          // 变形模式标志
    bool is_dragging_ = false;          // 拖拽状态
    ARAPSolver* solver_ = nullptr;      // ARAP求解器
    std::map<int, Eigen::Vector3d> drag_initial_positions_;
    Eigen::Vector3f drag_start_;        // 拖拽起点
    Eigen::Vector3f drag_current_;      // 当前拖拽点

    // 新增函数
    void handleDeformation(double xpos, double ypos);
    Eigen::Vector3f screenToWorld(double xpos, double ypos, float depth);
    void updateMeshBuffers();           // 更新网格缓冲区
    std::vector<Eigen::Vector3d> original_vertex_positions_;
    bool is_initialized_ = false;  // 添加初始化标志

    double last_click_x_ = 0.0;
    double last_click_y_ = 0.0;
    Eigen::Vector3f drag_start_world_;
    double drag_start_x_ = 0.0, drag_start_y_ = 0.0;
};