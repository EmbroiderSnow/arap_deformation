#pragma once
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <Eigen/Dense>
#include <vector>
#include <map>

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
    using VertexHandle = TriMesh::VertexHandle;

    bool load(const std::string& filename);
    bool save(const std::string& filename) const;

    TriMesh& getMesh() { return mesh_; }
    const TriMesh& getMesh() const { return mesh_; }

    // Handle selection methods
    void selectVertex(VertexHandle vh);
    void unselectVertex(VertexHandle vh);
    void clearSelection();
    bool isSelected(VertexHandle vh) const;

    // Handle constraints
    void setHandleConstraint(VertexHandle vh, const Eigen::Vector3d& position);
    void clearConstraints();
    std::map<VertexHandle, Eigen::Vector3d> getHandleConstraints() const;

    // Get selected vertices
    std::vector<VertexHandle> getSelectedVertices() const;

private:
    TriMesh mesh_;
    std::vector<VertexHandle> selected_vertices_;
    std::map<VertexHandle, Eigen::Vector3d> handle_constraints_;
};