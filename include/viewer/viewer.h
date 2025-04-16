#pragma once
#include "mesh/mesh.h"
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
};