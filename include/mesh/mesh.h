#pragma once
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <Eigen/Dense>

struct MeshTraits : public OpenMesh::DefaultTraits {
    using Point = OpenMesh::Vec3d;
    using Normal = OpenMesh::Vec3d;
    VertexAttributes(OpenMesh::Attributes::Status);
    FaceAttributes(OpenMesh::Attributes::Status);
    EdgeAttributes(OpenMesh::Attributes::Status);
};

class Mesh {
public:
    using TriMesh = OpenMesh::TriMesh_ArrayKernelT<MeshTraits>;

    bool load(const std::string& filename);
    bool save(const std::string& filename) const;

    TriMesh& getMesh() { return mesh_; }
    const TriMesh& getMesh() const { return mesh_; }

private:
    TriMesh mesh_;
};