#include "deformation/arap_solver.h"

// 构造函数：初始化权重和构建Laplacian矩阵
ARAPSolver::ARAPSolver(Mesh& mesh)
    : mesh_(mesh), w_rot_(1.0), w_pos_(1000.0) {
    buildWeightMatrix();
}

// 设置能量项权重的方法
void ARAPSolver::setWeight(double w_rot, double w_pos) {
    w_rot_ = w_rot;
    w_pos_ = w_pos;
}

void ARAPSolver::buildWeightMatrix() {
    const auto& trimesh = mesh_.getMesh();
    int n_vertices = trimesh.n_vertices();
    
    // Initialize Laplacian matrix
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(n_vertices * 7); // Approximate size

    // TODO: Build cotangent Laplacian matrix
    // We'll implement this next
}