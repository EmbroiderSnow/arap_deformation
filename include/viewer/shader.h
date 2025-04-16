#pragma once
#include <GL/glew.h>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

class Shader {
public:
    // 程序ID
    GLuint ID;

    // 构造器和析构器
    Shader(const char* vertexPath, const char* fragmentPath);
    ~Shader();

    // 使用/激活程序
    void use();
    // uniform工具函数
    void setBool(const std::string &name, bool value) const;
    void setInt(const std::string &name, int value) const;
    void setFloat(const std::string &name, float value) const;
    void setMat4(const std::string &name, const float* value) const;
    void setVec3(const std::string &name, float x, float y, float z) const;

private:
    void checkCompileErrors(GLuint shader, std::string type);
};