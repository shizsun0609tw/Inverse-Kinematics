#pragma once

#include <memory>

#include "Eigen/Core"
#include "graphics/sphere.h"

namespace graphics {
class Program;
}

namespace kinematics {
class Ball {
 public:
    Ball() noexcept;
    // no copy constructor
    Ball(const Ball&) = delete;
    Ball(Ball&&) noexcept;
    // You need this for alignment otherwise it may crash
    // Ref: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Ball& operator=(const Ball&) = delete;
    Ball& operator=(Ball&&) noexcept;
    void render(graphics::Program* program);
    Eigen::Vector4d& getCurrentPosition();
    // Calculate and set position
    void setModelMatrix();
    void setCurrentPosition(const Eigen::Vector4d& pos);

 private:
    Eigen::Vector4d current_position = Eigen::Vector4d::Zero();
    std::unique_ptr<graphics::Sphere> graphics;
};
}  // namespace kinematics
