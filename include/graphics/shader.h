#pragma once
#include <string>
#include <utility>

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "glad/gl.h"

#include "util/filesystem.h"

namespace graphics {
class Shader final {
 public:
    Shader(util::fs::path shader, GLenum shaderType);
    Shader(const Shader&) = delete;
    Shader& operator=(const Shader&) = delete;
    operator GLuint() const { return id; }
    ~Shader();

    GLuint getID() const;

 private:
    std::string loadFromFile(util::fs::path filename);
    GLuint id;
};

class Program final {
 public:
    Program();
    Program(const Program&) = delete;
    Program& operator=(const Program&) = delete;
    ~Program();

    template <typename... Args>
    void attachShader(Args&&... args) {
        (glAttachShader(this->id, args), ...);
    }
    void link();
    template <typename... Args>
    void attachLinkShader(Args&&... args) {
        this->attachShader(std::forward<Args>(args)...);
        this->link();
    }

    GLuint getID() const;
    int getUniformLocation(const char* name) const;
    void use() const;
    void setUniform(const char* name, int i1);
    void setUniform(const char* name, const Eigen::Matrix4f& mat4);
    void setUniform(const char* name, const Eigen::Affine3f& mat4);
    void setUniform(const char* name, const Eigen::Vector3f& vec3);
    void setUniform(const char* name, const Eigen::Vector4f& vec4);

 private:
    GLuint id;
};
}  // namespace graphics
