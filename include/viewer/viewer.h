#pragma once
#include "mesh/mesh.h"
#include "shader.h"
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <vector>

class MeshViewer {
public:
    MeshViewer(int width = 1280, int height = 720);
    ~MeshViewer();

    void run(Mesh& mesh);

private:
    void initGL();
    void renderMesh();
    void handleInput();
    void renderUI();

    GLFWwindow* window_;
    std::vector<Mesh::VertexHandle> selected_vertices_;
    bool dragging_ = false;
    Eigen::Vector3d drag_start_;

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
};