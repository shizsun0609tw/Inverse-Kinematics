#include "graphics/free_camera.h"

#include "GLFW/glfw3.h"

#include "util/helper.h"

namespace graphics {

float* FreeCamera::getMouseSensitivityPointer() { return &mouseSensitivity; }
float* FreeCamera::getMoveSpeedPointer() { return &moveSpeed; }

FreeCamera::FreeCamera() { position = Eigen::Vector4f(8.0f, 21.0f, -9.0f, 0.0f); }

void FreeCamera::moveCamera(GLFWwindow* window) {
    float currentFrameTime = static_cast<float>(glfwGetTime());
    float deltaTime = currentFrameTime - lastFrameTime;
    lastFrameTime = currentFrameTime;
    float speed = moveSpeed * deltaTime;
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        position += center * speed;
    else if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        position -= center * speed;
    else if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        position -= right * speed;
    else if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        position += right * speed;
}

void FreeCamera::moveSight(GLFWwindow* window) {
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);

    float x = static_cast<float>(xpos);
    float y = static_cast<float>(ypos);

    if (lastx == 0 && lasty == 0) {
        lastx = x;
        lasty = y;
    } else {
        float dx = x - lastx;
        float dy = lasty - y;
        lastx = x;
        lasty = y;
        dx *= mouseSensitivity;
        dy *= mouseSensitivity;
        yaw += dx;
        pitch += dy;
        pitch = std::clamp(pitch, -89.9f, 89.9f);
    }
}

void FreeCamera::updateProjectionMatrix() {
    projectionMatrix = util::perspective(fieldOfView, aspectRatio, nearClipPlane, farClipPlane);
}

void FreeCamera::reset() {
    lastx = 0;
    lasty = 0;
}

void FreeCamera::updateViewMatrix() {
    center[0] = cosf(util::toRadian(yaw)) * cosf(util::toRadian(pitch));
    center[1] = sinf(util::toRadian(pitch));
    center[2] = sinf(util::toRadian(yaw)) * cosf(util::toRadian(pitch));
    center.normalize();

    right = center.cross3(Eigen::Vector4f::UnitY());
    right.normalize();
    up = right.cross3(center);
    up.normalize();
    viewMatrix = util::lookAt(position, position + center, up);
}
}  // namespace graphics
