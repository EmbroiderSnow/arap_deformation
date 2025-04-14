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

double ARAPSolver::computeCotWeight(const Mesh::TriMesh::EdgeHandle& eh) const {
    auto& trimesh = mesh_.getMesh();
    double weight = 0.0;

    // Get halfedges
    auto heh = trimesh.halfedge_handle(eh, 0);
    auto heh_opp = trimesh.halfedge_handle(eh, 1);

    // For each halfedge, compute cotangent weight
    for (auto curr_heh : {heh, heh_opp}) {
        if (!trimesh.is_boundary(curr_heh)) {
            // Get vertices
            auto v0 = trimesh.from_vertex_handle(curr_heh);
            auto v1 = trimesh.to_vertex_handle(curr_heh);
            auto v2 = trimesh.to_vertex_handle(trimesh.next_halfedge_handle(curr_heh));

            // Get positions
            Eigen::Vector3d p0 = getPosition(v0);
            Eigen::Vector3d p1 = getPosition(v1);
            Eigen::Vector3d p2 = getPosition(v2);

            // Compute vectors
            Eigen::Vector3d e1 = p2 - p0;
            Eigen::Vector3d e2 = p2 - p1;

            // Compute cotangent
            double cos_angle = e1.dot(e2) / (e1.norm() * e2.norm());
            double sin_angle = (e1.cross(e2)).norm() / (e1.norm() * e2.norm());
            weight += cos_angle / sin_angle;
        }
    }

    return weight * 0.5;
}

Eigen::Vector3d ARAPSolver::getPosition(const Mesh::TriMesh::VertexHandle& vh) const {
    const auto& trimesh = mesh_.getMesh();
    const auto& p = trimesh.point(vh);
    return Eigen::Vector3d(p[0], p[1], p[2]);
}

void ARAPSolver::buildWeightMatrix() {
    const auto& trimesh = mesh_.getMesh();
    int n_vertices = trimesh.n_vertices();
    
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(n_vertices * 7);

    // Build cotangent Laplacian matrix
    for (auto vh : trimesh.vertices()) {
        double sum_weights = 0.0;
        int idx = vh.idx();

        // For each outgoing halfedge
        for (auto heh : trimesh.voh_range(vh)) {
            auto eh = trimesh.edge_handle(heh);
            auto vh_target = trimesh.to_vertex_handle(heh);
            int idx_target = vh_target.idx();

            // Compute cotangent weight
            double weight = computeCotWeight(eh);
            sum_weights += weight;

            // Add off-diagonal element
            triplets.emplace_back(idx, idx_target, -weight);
        }

        // Add diagonal element
        triplets.emplace_back(idx, idx, sum_weights);
    }

    // Create sparse matrix
    L_.resize(n_vertices, n_vertices);
    L_.setFromTriplets(triplets.begin(), triplets.end());
}

void ARAPSolver::initCellData() {
    auto& trimesh = mesh_.getMesh();
    int n_vertices = trimesh.n_vertices();
    
    cells_.resize(n_vertices);
    positions_ = Eigen::MatrixXd::Zero(n_vertices, 3);
    initial_positions_ = Eigen::MatrixXd::Zero(n_vertices, 3);
    rotations_.resize(n_vertices, Eigen::Matrix3d::Identity());

    // Initialize positions
    for (auto vh : trimesh.vertices()) {
        int idx = vh.idx();
        Eigen::Vector3d p = getPosition(vh);
        positions_.row(idx) = p;
        initial_positions_.row(idx) = p;
        
        Cell& cell = cells_[idx];
        cell.index = idx;

        // Get neighbors and weights
        for (auto voh_it = trimesh.voh_iter(vh); voh_it.is_valid(); ++voh_it) {
            auto eh = trimesh.edge_handle(*voh_it);
            auto vh_target = trimesh.to_vertex_handle(*voh_it);
            
            cell.neighbors.push_back(vh_target);
            cell.weights.push_back(computeCotWeight(eh));
        }
    }
}

void ARAPSolver::estimateRotations() {
    // For each vertex
    for (auto& cell : cells_) {
        Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
        
        // Compute covariance matrix
        for (size_t i = 0; i < cell.neighbors.size(); ++i) {
            int j = cell.neighbors[i].idx();
            double wij = cell.weights[i];
            
            Eigen::Vector3d pi = initial_positions_.row(cell.index);
            Eigen::Vector3d pj = initial_positions_.row(j);
            Eigen::Vector3d qi = positions_.row(cell.index);
            Eigen::Vector3d qj = positions_.row(j);
            
            covariance += wij * (pj - pi) * (qj - qi).transpose();
        }
        
        // SVD decomposition
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(covariance, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix3d U = svd.matrixU();
        Eigen::Matrix3d V = svd.matrixV();
        
        // Compute rotation
        cell.rotation = V * U.transpose();
        if (cell.rotation.determinant() < 0) {
            V.col(2) = -V.col(2);
            cell.rotation = V * U.transpose();
        }
        
        rotations_[cell.index] = cell.rotation;
    }
}

bool ARAPSolver::solvePosition() {
    const auto& trimesh = mesh_.getMesh();
    int n_vertices = trimesh.n_vertices();
    
    // Build right-hand side
    Eigen::MatrixXd b = Eigen::MatrixXd::Zero(n_vertices, 3);
    
    // For each vertex
    for (const auto& cell : cells_) {
        int i = cell.index;
        
        // Add rotated differential coordinates
        for (size_t j = 0; j < cell.neighbors.size(); ++j) {
            int neighbor_idx = cell.neighbors[j].idx();
            double wij = cell.weights[j];
            
            Eigen::Vector3d diff = initial_positions_.row(neighbor_idx) - initial_positions_.row(i);
            b.row(i) += 0.5 * wij * (rotations_[i] + rotations_[neighbor_idx]) * diff;
        }
    }
    
    // Add position constraints
    for (const auto& constraint : mesh_.getHandleConstraints()) {
        int idx = constraint.first.idx();
        b.row(idx) = w_pos_ * constraint.second;
        L_.coeffRef(idx, idx) += w_pos_;
    }
    
    // Solve the system
    solver_.compute(L_);
    if (solver_.info() != Eigen::Success) {
        return false;
    }
    
    positions_ = solver_.solve(b);
    return solver_.info() == Eigen::Success;
}

bool ARAPSolver::solve(int maxIterations) {
    // Initialize cell data if not done
    if (cells_.empty()) {
        initCellData();
    }
    
    // Local-global optimization
    for (int iter = 0; iter < maxIterations; ++iter) {
        // Local step: estimate rotations
        estimateRotations();
        
        // Global step: solve for positions
        if (!solvePosition()) {
            std::cerr << "Failed to solve global step at iteration " << iter << std::endl;
            return false;
        }
        
        // Update mesh vertices
        auto& trimesh = mesh_.getMesh();
        for (auto vh : trimesh.vertices()) {
            int idx = vh.idx();
            OpenMesh::Vec3d new_pos(positions_(idx, 0), 
                                  positions_(idx, 1), 
                                  positions_(idx, 2));
            trimesh.set_point(vh, new_pos);
        }
    }
    
    return true;
}