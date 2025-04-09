#include "mesh/mesh.h"

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