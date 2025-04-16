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
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW" << std::endl;
        return;
    }

    // 初始化 OpenGL 设置
    initGL();
}

MeshViewer::~MeshViewer() {
    if (window_) {
        glfwDestroyWindow(window_);
    }
    glfwTerminate();
}

void MeshViewer::initGL() {
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
}

void MeshViewer::renderMesh() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // TODO: 实现网格渲染
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
    while (!glfwWindowShouldClose(window_)) {
        handleInput();
        renderMesh();
        renderUI();

        glfwSwapBuffers(window_);
        glfwPollEvents();
    }
}