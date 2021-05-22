#include "graphics/default_camera.h"

#include "util/helper.h"

namespace graphics {
float* DefaultCamera::getCameraRotationAnglePointer() { return &cameraRotationAngle; }
float* DefaultCamera::getCameraYOffsetPointer() { return &cameraYOffset; }
float* DefaultCamera::getCameraRotationRadiusPointer() { return &cameraRotationRadius; }

void DefaultCamera::updateProjectionMatrix() {
    projectionMatrix = util::perspective(fieldOfView, aspectRatio, nearClipPlane, farClipPlane);
}

void DefaultCamera::updateViewMatrix() {
    position[0] = center[0] + cameraRotationRadius * cosf(util::toRadian(cameraRotationAngle));
    position[1] = (center[1] + cameraYOffset) * cameraRotationRadius;
    position[2] = center[2] + cameraRotationRadius * sinf(util::toRadian(cameraRotationAngle));
    viewMatrix = util::lookAt(position, center, up);
}
}  // namespace graphics
