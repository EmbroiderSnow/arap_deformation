cmake_minimum_required(VERSION 3.10)
project(arap-deformation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(extern/OpenMesh)
add_subdirectory(extern/eigen)

# Find required packages
find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)
find_package(GLEW REQUIRED)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include  # 使用 CMAKE_CURRENT_SOURCE_DIR 确保路径正确
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/OpenMesh/src
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/eigen
    ${GLFW3_INCLUDE_DIRS}
    ${GLEW_INCLUDE_DIRS}
)

# Add shader directory to resources
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/shaders DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(GLOB_RECURSE SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
)

add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} 
    OpenMeshCore 
    OpenMeshTools
    OpenGL::GL
    glfw
    GLEW::GLEW
)

# Enable C++11 for Eigen
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_11)