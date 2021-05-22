#pragma once

#include <vector>

#include "Eigen/Core"
#include "Eigen/StdVector"

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Vector4d)

namespace acclaim {
struct Posture final {
    Posture() noexcept;
    explicit Posture(const std::size_t size) noexcept;
    Posture(const Posture &) noexcept;
    Posture(Posture &&) noexcept;

    Posture &operator=(const Posture &) noexcept;
    Posture &operator=(Posture &&) noexcept;
    // You need this for alignment otherwise it may crash
    // Ref: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    std::vector<Eigen::Vector4d> bone_rotations;
    std::vector<Eigen::Vector4d> bone_translations;
};
}  // namespace acclaim
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(acclaim::Posture)
