cmake_minimum_required(VERSION 3.10)
project(arap-deformation)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(extern/OpenMesh)
add_subdirectory(extern/eigen)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include  # 使用 CMAKE_CURRENT_SOURCE_DIR 确保路径正确
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/OpenMesh/src
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/eigen
)

file(GLOB_RECURSE SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
)

add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} 
    OpenMeshCore 
    OpenMeshTools
)

# Enable C++11 for Eigen
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_11)