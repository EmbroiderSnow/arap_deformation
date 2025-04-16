#include "mesh/mesh.h"
#include "deformation/arap_solver.h"
#include "viewer/viewer.h"
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

    // ...existing code...

    #ifdef DEBUG
    MeshViewer viewer;
    viewer.run(mesh);  // This will handle the interaction loop
    #endif

    return 0;
}