#pragma once
#include <memory>

#include "buffer.h"
#include "rigidbody.h"
namespace graphics {

class Cylinder final : public Rigidbody {
 public:
    Cylinder() noexcept;
    Cylinder(const Cylinder&) noexcept;
    Cylinder(Cylinder&&) noexcept;
    // You need this for alignment otherwise it may crash
    // Ref: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Cylinder& operator=(const Cylinder&) noexcept;
    Cylinder& operator=(Cylinder&&) noexcept;

    void bindVBO() const override;
    void generateVertices() override;
    void render(Program* shaderProgram) override;

 private:
    std::shared_ptr<Buffer<1, GL_ARRAY_BUFFER>> vbo = nullptr;
    std::shared_ptr<Buffer<2, GL_ELEMENT_ARRAY_BUFFER>> ebo = nullptr;
    static std::weak_ptr<Buffer<1, GL_ARRAY_BUFFER>> vbo_weak;
    static std::weak_ptr<Buffer<2, GL_ELEMENT_ARRAY_BUFFER>> ebo_weak;
    static bool isInitialized;
};
}  // namespace graphics
