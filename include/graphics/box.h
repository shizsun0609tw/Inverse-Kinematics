#pragma once
#include <memory>

#include "buffer.h"
#include "rigidbody.h"

namespace graphics {
class Box final : public Rigidbody {
 public:
    Box() noexcept;
    Box(const Box&) noexcept;
    Box(Box&&) noexcept;
    // You need this for alignment otherwise it may crash
    // Ref: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Box& operator=(const Box&) noexcept;
    Box& operator=(Box&&) noexcept;

    void bindVBO() const override;
    void generateVertices() override;
    void render(Program* shaderProgram) override;

 private:
    void initialize();
    std::shared_ptr<Buffer<1, GL_ARRAY_BUFFER>> vbo = nullptr;
    std::shared_ptr<Buffer<1, GL_ELEMENT_ARRAY_BUFFER>> ebo = nullptr;
    static std::weak_ptr<Buffer<1, GL_ARRAY_BUFFER>> vbo_weak;
    static std::weak_ptr<Buffer<1, GL_ELEMENT_ARRAY_BUFFER>> ebo_weak;
    static bool isInitialized;
};
}  // namespace graphics
