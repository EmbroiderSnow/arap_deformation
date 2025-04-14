#include "mesh/mesh.h"
#include <algorithm>
#include <iostream>

bool Mesh::load(const std::string& filename) {
    if (!OpenMesh::IO::read_mesh(mesh_, filename)) {
        std::cerr << "Error: Cannot read mesh from file " << filename << std::endl;
        return false;
    }
    mesh_.request_vertex_status();
    mesh_.request_face_status();
    mesh_.request_edge_status();
    return true;
}

bool Mesh::save(const std::string& filename) const {
    if (!OpenMesh::IO::write_mesh(mesh_, filename)) {
        std::cerr << "Error: Cannot write mesh to file " << filename << std::endl;
        return false;
    }
    return true;
}

void Mesh::selectVertex(VertexHandle vh) {
    if (!mesh_.is_valid_handle(vh)) {
        std::cerr << "Error: Invalid vertex handle." << std::endl;
        return;
    }

    mesh_.status(vh).set_selected(true);
    if (std::find(selected_vertices_.begin(), selected_vertices_.end(), vh) == selected_vertices_.end()) {
        selected_vertices_.push_back(vh);
    }
}

void Mesh::unselectVertex(VertexHandle vh) {
    if (!mesh_.is_valid_handle(vh)) {
        std::cerr << "Error: Invalid vertex handle." << std::endl;
        return;
    }

    mesh_.status(vh).set_selected(false);
    auto it = std::find(selected_vertices_.begin(), selected_vertices_.end(), vh);
    if (it != selected_vertices_.end()) {
        selected_vertices_.erase(it);
    }
}

void Mesh::clearSelection() {
    for (auto vh : selected_vertices_) {
        mesh_.status(vh).set_selected(false);
    }
    selected_vertices_.clear();
}

bool Mesh::isSelected(VertexHandle vh) const {
    if (!mesh_.is_valid_handle(vh)) {
        std::cerr << "Error: Invalid vertex handle." << std::endl;
        return false;
    }
    return mesh_.status(vh).selected();
}

void Mesh::setHandleConstraint(VertexHandle vh, const Eigen::Vector3d& position) {
    if (!mesh_.is_valid_handle(vh)) {
        std::cerr << "Error: Invalid vertex handle." << std::endl;
        return;
    }
    handle_constraints_[vh] = position;
}

void Mesh::clearConstraints() {
    handle_constraints_.clear();
}

std::map<Mesh::VertexHandle, Eigen::Vector3d> Mesh::getHandleConstraints() const {
    return handle_constraints_;
}

std::vector<Mesh::VertexHandle> Mesh::getSelectedVertices() const {
    return selected_vertices_;
}