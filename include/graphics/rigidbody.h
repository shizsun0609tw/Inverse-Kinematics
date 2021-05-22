#pragma once
#include <memory>

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "glad/gl.h"

#include "texture.h"

namespace graphics {
class Program;
class Rigidbody {
 public:
    Rigidbody() noexcept;
    Rigidbody(const Rigidbody& other) noexcept;
    Rigidbody(Rigidbody&& other) noexcept;
    // You need this for alignment otherwise it may crash
    // Ref: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Rigidbody& operator=(const Rigidbody& other) noexcept;
    Rigidbody& operator=(Rigidbody&& other) noexcept;
    void setModelMatrix(const Eigen::Affine3f& _modelMatrix);
    void setTexture(const std::shared_ptr<TextureBase>& texture);
    void setTexture(const Eigen::Vector4f& color);

    void initialize();

    virtual void generateVertices() = 0;
    virtual void render(Program* shaderProgram) = 0;

 protected:
    virtual ~Rigidbody();
    virtual void bindVBO() const = 0;
    GLuint vao = 0;

    Eigen::Affine3f modelMatrix = Eigen::Affine3f::Identity();
    Eigen::Matrix4f inverseTransposeModel = Eigen::Matrix4f::Identity();
    Eigen::Vector4f baseColor = Eigen::Vector4f::UnitX();
    std::shared_ptr<TextureBase> texture = nullptr;
};
}  // namespace graphics
