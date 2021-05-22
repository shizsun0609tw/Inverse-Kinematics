#pragma once
#include <vector>

#include "acclaim/posture.h"

namespace acclaim {
struct Bone;
}
namespace kinematics {
// Apply forward kinematics to skeleton
void forwardSolver(const acclaim::Posture& posture, acclaim::Bone* bone);

Eigen::VectorXd pseudoInverseLinearSolver(const Eigen::Matrix4Xd& Jacobian, const Eigen::Vector4d& target);

bool inverseJacobianIKSolver(const Eigen::Vector4d& target_pos, acclaim::Bone* start_bone, acclaim::Bone* end_bone,
                             acclaim::Posture& posture);
}  // namespace kinematics
