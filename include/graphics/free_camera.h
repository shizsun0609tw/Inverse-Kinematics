#pragma once
#include "Eigen/StdVector"

#include "camera.h"

struct GLFWwindow;

namespace graphics {
class FreeCamera final : public Camera {
 public:
    FreeCamera();
    // You need this for alignment otherwise it may crash
    // Ref: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    float* getMoveSpeedPointer();
    float* getMouseSensitivityPointer();
    void reset();
    void moveSight(GLFWwindow* window);
    void moveCamera(GLFWwindow* window);

 protected:
    void updateViewMatrix() override;
    void updateProjectionMatrix() override;

 private:
    float moveSpeed = 30, mouseSensitivity = 0.1f;
    float lastx = 0, lasty = 0;
    float yaw = 135, pitch = -45;
    float lastFrameTime = 0.0f;
    Eigen::Vector4f right = Eigen::Vector4f::Zero();
};
}  // namespace graphics

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(graphics::FreeCamera)
