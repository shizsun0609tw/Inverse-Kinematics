#include "util/helper.h"
#include <cmath>

namespace util {
namespace {
Eigen::AngleAxisd rotateRadianX(double x) { return Eigen::AngleAxisd(x, Eigen::Vector3d::UnitX()); }
Eigen::AngleAxisd rotateRadianY(double y) { return Eigen::AngleAxisd(y, Eigen::Vector3d::UnitY()); }
Eigen::AngleAxisd rotateRadianZ(double z) { return Eigen::AngleAxisd(z, Eigen::Vector3d::UnitZ()); }
}  // namespace

Eigen::Vector4d toDegree(const Eigen::Vector4d& angle) { return angle * 180.0 / PI; }
Eigen::Vector4d toRadian(const Eigen::Vector4d& angle) { return angle * PI / 180.0; }
// https://github.com/OpenGP/htrack/blob/master/util/eigen_opengl_helpers.h

Eigen::Matrix4f perspective(float fovy, float aspect, float zNear, float zFar) {
    Eigen::Projective3f tr;
    tr.matrix().setZero();
    float tan_half_fovy = tanf(toRadian(fovy) / 2.0f);
    tr(0, 0) = 1.0f / (aspect * tan_half_fovy);
    tr(1, 1) = 1.0f / tan_half_fovy;
    tr(2, 2) = -(zFar + zNear) / (zFar - zNear);
    tr(3, 2) = -1.0f;
    tr(2, 3) = -(2.0f * zFar * zNear) / (zFar - zNear);
    return tr.matrix();
}
Eigen::Matrix4f lookAt(const Eigen::Vector4f& eye, const Eigen::Vector4f& center, const Eigen::Vector4f& up) {
    Eigen::Vector4f f = center - eye;
    f.normalize();
    Eigen::Vector4f u = up;
    u.normalize();
    Eigen::Vector4f s = f.cross3(u);
    s.normalize();
    u = s.cross3(f);
    Eigen::Matrix4f mat;
    mat.row(0) = s;
    mat.row(1) = u;
    mat.row(2) = -f;
    mat.row(3) = Eigen::Vector4f::UnitW();
    mat(0, 3) = -s.dot(eye);
    mat(1, 3) = -u.dot(eye);
    mat(2, 3) = f.dot(eye);
    return mat;
}
Eigen::Matrix4f ortho(float left, float right, float bottom, float top, float zNear, float zFar) {
    Eigen::Matrix4f mat = Eigen::Matrix4f::Identity();
    mat(0, 0) = 2.0f / (right - left);
    mat(1, 1) = 2.0f / (top - bottom);
    mat(2, 2) = -2.0f / (zFar - zNear);
    mat(3, 0) = -(right + left) / (right - left);
    mat(3, 1) = -(top + bottom) / (top - bottom);
    mat(3, 2) = -(zFar + zNear) / (zFar - zNear);
    return mat;
}
Eigen::Quaterniond rotateDegreeZYX(double x, double y, double z) {
    return rotateRadianZYX(toRadian(x), toRadian(y), toRadian(z));
}
Eigen::Quaterniond rotateDegreeZYX(const Eigen::Vector4d& rotation) { return rotateRadianZYX(toRadian(rotation)); }
Eigen::Quaterniond rotateDegreeXYZ(double x, double y, double z) {
    return rotateRadianXYZ(toRadian(x), toRadian(y), toRadian(z));
}
Eigen::Quaterniond rotateDegreeXYZ(const Eigen::Vector4d& rotation) { return rotateRadianXYZ(toRadian(rotation)); }
Eigen::Quaterniond rotateRadianZYX(double x, double y, double z) {
    return rotateRadianZ(z) * rotateRadianY(y) * rotateRadianX(x);
}
Eigen::Quaterniond rotateRadianZYX(const Eigen::Vector4d& rotation) {
    return rotateRadianZ(rotation[2]) * rotateRadianY(rotation[1]) * rotateRadianX(rotation[0]);
}
Eigen::Quaterniond rotateRadianXYZ(double x, double y, double z) {
    return rotateRadianX(x) * rotateRadianY(y) * rotateRadianZ(z);
}
Eigen::Quaterniond rotateRadianXYZ(const Eigen::Vector4d& rotation) {
    return rotateRadianX(rotation[0]) * rotateRadianY(rotation[1]) * rotateRadianZ(rotation[2]);
}
}  // namespace util
