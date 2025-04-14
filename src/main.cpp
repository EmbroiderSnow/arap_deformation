#include "mesh/mesh.h"
#include "deformation/arap_solver.h"
#include <iostream>
#define DEBUG

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <mesh_file>" << std::endl;
        return 1;
    }

    Mesh mesh;
    if (!mesh.load(argv[1])) {
        std::cerr << "Failed to load mesh." << std::endl;
        return 1;
    }

    std::cout << "Loaded mesh with "
              << mesh.getMesh().n_vertices() << " vertices, "
              << mesh.getMesh().n_faces() << " faces, and "
              << mesh.getMesh().n_edges() << " edges." << std::endl;

    #ifdef DEBUG
    std::cout << "Debug mode is enabled." << std::endl;
    // 测试顶点选择
    auto& trimesh = mesh.getMesh();
    Mesh::VertexHandle vh = trimesh.vertex_handle(0);
    mesh.selectVertex(vh);
    
    if (!mesh.isSelected(vh)) {
        std::cerr << "Vertex selection failed" << std::endl;
        return 1;
    }

    // Set a test constraint
    Eigen::Vector3d new_pos(1.0, 0.0, 0.0);
    mesh.setHandleConstraint(vh, new_pos);
    
    // Initialize and run solver
    ARAPSolver solver(mesh);
    if (!solver.solve(5)) {
        std::cerr << "ARAP solving failed" << std::endl;
        return 1;
    }
    
    // Save deformed mesh
    if (!mesh.save("deformed.obj")) {
        std::cerr << "Failed to save deformed mesh" << std::endl;
        return 1;
    }
    
    std::cout << "All tests passed!" << std::endl;
    #endif

    return 0;
}