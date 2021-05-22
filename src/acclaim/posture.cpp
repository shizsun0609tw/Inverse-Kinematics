#include "acclaim/posture.h"

#include <utility>

namespace acclaim {
Posture::Posture() noexcept {}

Posture::Posture(const std::size_t size) noexcept
    : bone_rotations(size, Eigen::Vector4d::Zero()), bone_translations(size, Eigen::Vector4d::Zero()) {}

Posture::Posture(const Posture &other) noexcept
    : bone_rotations(other.bone_rotations), bone_translations(other.bone_translations) {}

Posture::Posture(Posture &&other) noexcept
    : bone_rotations(std::move(other.bone_rotations)), bone_translations(std::move(other.bone_translations)) {}

Posture &Posture::operator=(const Posture &other) noexcept {
    if (this != &other) {
        bone_rotations = other.bone_rotations;
        bone_translations = other.bone_translations;
    }
    return *this;
}
Posture &Posture::operator=(Posture &&other) noexcept {
    if (this != &other) {
        bone_rotations = std::move(other.bone_rotations);
        bone_translations = std::move(other.bone_translations);
    }
    return *this;
}
}  // namespace acclaim
