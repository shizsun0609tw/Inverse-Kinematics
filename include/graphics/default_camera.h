#pragma once
#include "Eigen/StdVector"

#include "camera.h"

namespace graphics {
class DefaultCamera final : public Camera {
 public:
    // You need this for alignment otherwise it may crash
    // Ref: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    float* getCameraRotationAnglePointer();
    float* getCameraRotationRadiusPointer();
    float* getCameraYOffsetPointer();

 protected:
    void updateViewMatrix() override;
    void updateProjectionMatrix() override;

 private:
    float cameraRotationAngle = 60.0f;
    float cameraYOffset = 0.2f;
    float cameraRotationRadius = 17.0f;
};
}  // namespace graphics

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(graphics::DefaultCamera)
