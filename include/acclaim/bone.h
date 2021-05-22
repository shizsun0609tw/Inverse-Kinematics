#pragma once
#include <string>

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Eigen/StdVector"

namespace acclaim {
// Bone segment names used in ASF file
// this structure defines the property of each bone segment, including its
// connection to other bones, DOF (degrees of freedom), relative orientation and
// distance to the outboard bone
struct Bone final {
    // You need this for alignment otherwise it may crash
    // Ref: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    // Note: Affine3d is actually a 4*4 matrix!
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    // Pointer to the sibling (branch bone) in the hierarchy tree
    Bone *sibling = nullptr;
    // Pointer to the child (outboard bone) in the hierarchy tree
    Bone *child = nullptr;
    // Pointer to the parent (inboard bone) in the hierarchy tree
    Bone *parent = nullptr;
    // Bone index
    int idx = 0;
    // Unit vector describes the direction from local origin to
    // the origin of the child bone
    // Notice: stored in local coordinate system of the bone
    Eigen::Vector4d dir = Eigen::Vector4d::Zero();
    // Bone name
    std::string name = "";
    // Bone length
    double length = 0.0;
    // Orientation of each bone's local coordinate
    // system as specified in ASF file (axis field)
    Eigen::Vector4d axis = Eigen::Vector4d::Zero();
    // number of bone's degrees of freedom
    int dof = 0;
    // degree of freedom mask in x, y, z axis (local)
    bool dofrx = false, dofry = false, dofrz = false;  // Rotate
    bool doftx = false, dofty = false, doftz = false;  // Translate
    // Rotation matrix from parent to child
    Eigen::Affine3d rot_parent_current = Eigen::Affine3d::Identity();
    // Initial rotation and scaling for bone
    Eigen::Affine3d global_facing = Eigen::Affine3d::Identity();
    // Bone's start pos in global position
    Eigen::Vector4d start_position = Eigen::Vector4d::Zero();
    // Bone's end pos in global position
    Eigen::Vector4d end_position = Eigen::Vector4d::Zero();
    // Bone's rotation in global position
    Eigen::Affine3d rotation = Eigen::Affine3d::Identity();
};
}  // namespace acclaim
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(acclaim::Bone)
