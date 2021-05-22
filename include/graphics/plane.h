#pragma once
#include <memory>

#include "buffer.h"
#include "rigidbody.h"

class Program;

namespace graphics {
class Plane final : public Rigidbody {
 public:
    Plane() noexcept;
    Plane(const Plane&) noexcept;
    Plane(Plane&&) noexcept;
    // You need this for alignment otherwise it may crash
    // Ref: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Plane& operator=(const Plane&) noexcept;
    Plane& operator=(Plane&&) noexcept;
    void bindVBO() const override;
    void generateVertices() override;
    void render(Program* shaderProgram) override;

 private:
    std::shared_ptr<Buffer<1, GL_ARRAY_BUFFER>> vbo = nullptr;
    static std::weak_ptr<Buffer<1, GL_ARRAY_BUFFER>> vbo_weak;
    static bool isInitialized;
};
}  // namespace graphics
