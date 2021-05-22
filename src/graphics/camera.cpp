#include "graphics/camera.h"

namespace graphics {
float Camera::getAspectRatio() const { return aspectRatio; }
Eigen::Vector4f Camera::getPosition() const { return position; }
Eigen::Vector4f Camera::getCenter() const { return center; }
Eigen::Matrix4f Camera::getViewMatrix() const { return viewMatrix; }
Eigen::Matrix4f Camera::getProjectionMatrix() const { return projectionMatrix; }
Eigen::Matrix4f Camera::getViewWithProjectionMatrix() const { return VPMatrix; }

void Camera::setAspectRatio(int width, int height) { aspectRatio = static_cast<float>(width) / height; }
void Camera::setPosition(const Eigen::Vector4f& _position) { position = _position; }
void Camera::setCenter(const Eigen::Vector4f& _center) { center = _center; }

void Camera::update() {
    updateViewMatrix();
    updateProjectionMatrix();
    VPMatrix = projectionMatrix * viewMatrix;
}
}  // namespace graphics
