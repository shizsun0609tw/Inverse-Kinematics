#pragma once

#include "Eigen/Core"

namespace graphics {
class Camera {
 public:
    // You need this for alignment otherwise it may crash
    // Ref: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    float getAspectRatio() const;
    Eigen::Vector4f getCenter() const;
    Eigen::Vector4f getPosition() const;
    Eigen::Matrix4f getViewMatrix() const;
    Eigen::Matrix4f getProjectionMatrix() const;
    Eigen::Matrix4f getViewWithProjectionMatrix() const;

    void setAspectRatio(int width, int height);
    void setPosition(const Eigen::Vector4f& _position);
    void setCenter(const Eigen::Vector4f& _center);

    void update();

 protected:
    virtual ~Camera() = default;
    // For calulating view and projection matrices
    Eigen::Vector4f position = Eigen::Vector4f::Zero();
    Eigen::Vector4f up = Eigen::Vector4f::UnitY();
    Eigen::Vector4f center = Eigen::Vector4f::UnitZ();
    // Cache for those matrices
    Eigen::Matrix4f viewMatrix = Eigen::Matrix4f::Identity();
    Eigen::Matrix4f projectionMatrix = Eigen::Matrix4f::Identity();
    Eigen::Matrix4f VPMatrix = Eigen::Matrix4f::Identity();
    float fieldOfView = 60.0f;
    float nearClipPlane = 0.01f;
    float farClipPlane = 100.0f;
    float aspectRatio = 1.0f;
    // Updater
    virtual void updateViewMatrix() = 0;
    virtual void updateProjectionMatrix() = 0;
};
}  // namespace graphics
