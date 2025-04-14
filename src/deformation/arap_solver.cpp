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
    const auto& trimesh = mesh_.getMesh();
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