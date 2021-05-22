#include "simulation/ball.h"

#include <utility>

namespace kinematics {
Ball::Ball() noexcept
    : current_position(-0.0692501, 3.85358, -1.63441, 0.0), graphics(std::make_unique<graphics::Sphere>()) {
    graphics->setTexture(Eigen::Vector4f(0.0f, 0.5f, 1.0f, 0.0f));
}

Ball::Ball(Ball&& other) noexcept
    : current_position(std::move(other.current_position)), graphics(std::move(other.graphics)) {}

Ball& Ball::operator=(Ball&& other) noexcept {
    if (this != &other) {
        current_position = std::move(other.current_position);
        graphics = std::move(other.graphics);
    }
    return *this;
}

void Ball::render(graphics::Program* program) { graphics->render(program); }

Eigen::Vector4d& Ball::getCurrentPosition() { return current_position; }

void Ball::setModelMatrix() {
    Eigen::Affine3d trans = Eigen::Affine3d::Identity();
    trans.translate(current_position.head<3>());
    trans.scale(Eigen::Vector3d(0.125, 0.125, 0.125));
    graphics->setModelMatrix(trans.cast<float>());
}

void Ball::setCurrentPosition(const Eigen::Vector4d& pos) { current_position = pos; }
}  // namespace kinematics
