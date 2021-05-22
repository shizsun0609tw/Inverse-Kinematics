#pragma once
#include "Eigen/Core"
#include "Eigen/Geometry"

namespace util {
// Math constant PI
constexpr double PI = 3.141592653589793238;
// Degree to radian
constexpr float toRadian(float angle) { return static_cast<float>(angle * PI / 180.0); }
// Degree to radian
constexpr double toRadian(double angle) { return angle * PI / 180.0; }
// Radian to degree
Eigen::Vector4d toDegree(const Eigen::Vector4d& angle);
// Degree to radian
Eigen::Vector4d toRadian(const Eigen::Vector4d& angle);
// Generate perspective matrix
Eigen::Matrix4f perspective(float fovy, float aspect, float zNear, float zFar);
// Generate lookat(view) matrix
Eigen::Matrix4f lookAt(const Eigen::Vector4f& eye, const Eigen::Vector4f& center, const Eigen::Vector4f& up);
// Generate orthogonal projection matrix
Eigen::Matrix4f ortho(float left, float right, float bottom, float top, float zNear, float zFar);
// Rotate along X axis first then Y axis then Z axis
Eigen::Quaterniond rotateDegreeZYX(double x, double y, double z);
// Rotate along X axis first then Y axis then Z axis
Eigen::Quaterniond rotateDegreeZYX(const Eigen::Vector4d& rotation);
// Rotate along Z axis first then Y axis then X axis
Eigen::Quaterniond rotateDegreeXYZ(double x, double y, double z);
// Rotate along Z axis first then Y axis then X axis
Eigen::Quaterniond rotateDegreeXYZ(const Eigen::Vector4d& rotation);
// Rotate along X axis first then Y axis then Z axis
Eigen::Quaterniond rotateRadianZYX(double x, double y, double z);
// Rotate along X axis first then Y axis then Z axis
Eigen::Quaterniond rotateRadianZYX(const Eigen::Vector4d& rotation);
// Rotate along Z axis first then Y axis then X axis
Eigen::Quaterniond rotateRadianXYZ(double x, double y, double z);
// Rotate along Z axis first then Y axis then X axis
Eigen::Quaterniond rotateRadianXYZ(const Eigen::Vector4d& rotation);
}  // namespace util
