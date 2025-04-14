#pragma once
#include "mesh/mesh.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>

class ARAPSolver {
public:
    explicit ARAPSolver(Mesh& mesh);

    // Set weight for rotation and position constraints
    void setWeight(double w_rot = 1.0, double w_pos = 1000.0);
    
    // Main solve function
    bool solve(int maxIterations = 5);

private:
    // Build cotangent weight matrix
    void buildWeightMatrix();
    
    // Local step: estimate rotations
    void estimateRotations();
    
    // Global step: solve for positions
    bool solvePosition();

    // Calculate cotangent weight for an edge
    double computeCotWeight(const Mesh::TriMesh::EdgeHandle& eh) const;
    
    // Get vertex position as Eigen vector
    Eigen::Vector3d getPosition(const Mesh::TriMesh::VertexHandle& vh) const;

    Mesh& mesh_;
    Eigen::SparseMatrix<double> L_;        // Laplacian matrix
    std::vector<Eigen::Matrix3d> rotations_; // Per-vertex rotations
    double w_rot_, w_pos_;                 // Weights for energy terms
};