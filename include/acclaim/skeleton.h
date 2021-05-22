#pragma once
#include <string>
#include <vector>

#include "bone.h"
#include "graphics/cylinder.h"
#include "util/filesystem.h"

namespace graphics {
class Program;
}

namespace acclaim {
class Skeleton final {
 public:
    // Root always has index 0
    static constexpr int root_idx() noexcept { return 0; }
    Skeleton(const util::fs::path &file_name, const double _scale) noexcept;
    Skeleton(const Skeleton &) noexcept;
    Skeleton(Skeleton &&) noexcept;

    Skeleton &operator=(const Skeleton &) noexcept;
    Skeleton &operator=(Skeleton &&) noexcept;
    // You need this for alignment otherwise it may crash
    // Ref: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    // get skeleton's scale
    double getScale() const;
    // get total bones in the skeleton
    int getBoneNum() const;
    // get total movable bones in the skeleton
    int getMovableBoneNum() const;
    // get specific bone by its name
    Bone *getBonePointer(const std::string &name);
    // get specific bone by its index
    Bone *getBonePointer(const int bone_idx);
    // set bone's color (for rendering)
    void setBoneColor(const Eigen::Vector4f &boneColor) const;
    // set bone's model matrices (for rendering)
    void setModelMatrices() const;
    // render the bone
    void render(graphics::Program *program) const;

 private:
    bool readASFFile(const util::fs::path &file_name);
    // This function sets sibling or child for parent bone
    // If parent bone does not have a child,
    // then child is set as parent's child
    // else child is set as a sibling of parent's already existing child
    void setBoneHierarchy(Bone *parent, Bone *child);
    // Transform the direction vector (dir),
    // which is defined in character's global coordinate system in the ASF file,
    // to local coordinate
    void computeLocalDirection();
    // Calculate rotation from each bone local coordinate system to the coordinate system of its parent
    // store it in rot_parent_current variable for each bone
    void computeLocalRotation();
    // setup graphics
    void setBoneGraphics();

    double scale = 0.2;
    int movableBones = 1;
    std::vector<Bone> bones = std::vector<Bone>(1);
    mutable std::vector<graphics::Cylinder> bone_graphics;
};
}  // namespace acclaim
