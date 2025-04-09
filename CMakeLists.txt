cmake_minimum_required(VERSION 3.10)
project(arap-deformation)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(
    ${project_SOURCE_DIR}/include
    ${project_SOURCE_DIR}/extern/OpenMesh/src
    ${project_SOURCE_DIR}/extern/eigen
)

file(
    GLOB_RECURSE
    SOURCES
    "${project_SOURCE_DIR}/src/*.cpp"
    "${project_SOURCE_DIR}/src/*.h"
)

add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} 
    OpenMeshCore 
    OpenMeshTools
)