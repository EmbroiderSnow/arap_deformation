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

    // Cell data for local step
    struct Cell {
        std::vector<Mesh::TriMesh::VertexHandle> neighbors;
        std::vector<double> weights;
        Eigen::Matrix3d rotation;
        int index;
    };

    // Initialize per-vertex cell data
    void initCellData();
    
    // Solve system using Eigen's sparse solver
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver_;
    std::vector<Cell> cells_;
    Eigen::MatrixXd positions_;  // Current positions
    Eigen::MatrixXd initial_positions_;  // Original positions
};