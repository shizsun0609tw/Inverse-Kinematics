#include "simulation/kinematics.h"

#include "Eigen/Dense"

#include "acclaim/bone.h"
#include "util/helper.h"

#include <iostream>

namespace kinematics {
Eigen::VectorXd pseudoInverseLinearSolver(const Eigen::Matrix4Xd& Jacobian, const Eigen::Vector4d& target) {
    // J_plus * V = theta
	// J_plus = JT * (J * JT)^-1 (row linear independent)
    Eigen::Matrix3Xd J = Jacobian.topRows(3);
    Eigen::Vector3d V = target.head<3>();
	Eigen::MatrixXd J_plus = J.transpose() * (J * J.transpose()).inverse();
    Eigen::VectorXd theta = J_plus * V;
	
    return theta;
}

/**
 * @brief Perform inverse kinematics (IK)
 *
 * @param target_pos The position where `end_bone` will move to.
 * @param start_bone This bone is the last bone you can move while doing IK
 * @param end_bone This bone will try to reach `target_pos`
 * @param posture The original AMC motion's reference, you need to modify this
 *
 * @return True if IK is stable (HW3 bonus)
 */
bool inverseJacobianIKSolver(const Eigen::Vector4d& target_pos, acclaim::Bone* start_bone, acclaim::Bone* end_bone,
                             acclaim::Posture& posture) {
    constexpr int max_iteration = 1000;
    constexpr double epsilon = 1E-3;
    constexpr double step = 0.1;
    // Since bone stores in bones[i] that i == bone->idx, we can use bone - bone->idx to find bones[0] which is root.
    acclaim::Bone* root_bone = start_bone - start_bone->idx;
    size_t bone_num = 0;
    // HINT:
    // calculate number of bones need to move to perform IK, store in `bone_num`
    // a.k.a. how may bones from end_bone to its parent than to start_bone (include both side)
    acclaim::Bone* current_bone = end_bone;
	while (current_bone != NULL)
	{
        bone_num++;
        if (current_bone == start_bone) break;
        current_bone = current_bone->parent;
	}

    Eigen::Vector4d past_error = target_pos - end_bone->end_position;
    Eigen::Matrix4Xd Jacobian(4, 3 * bone_num);
    Jacobian.setZero();
    for (int iter = 0; iter < max_iteration; ++iter) {
        forwardSolver(posture, root_bone);
        Eigen::Vector4d desiredVector = target_pos - end_bone->end_position;
        if (desiredVector.norm() < epsilon) {
            break;
        }

        // HINT:
        // Calculate Jacobian, store in `Jacobian`
        int count = 0;
        current_bone = end_bone;
		while (current_bone != NULL) 
		{
            Eigen::Vector4d alpha = Eigen::Vector4d::Zero();
            Eigen::Vector4d r = end_bone->end_position - current_bone->start_position;
            
			alpha = current_bone->dofrx == true ? (current_bone->rotation * Eigen::Vector4d(1, 0, 0, 0)).normalized() : Eigen::Vector4d::Zero();
            Jacobian.col(count * 3 + 0) = alpha.cross3(r);
			
            alpha = current_bone->dofry == true ? (current_bone->rotation * Eigen::Vector4d(0, 1, 0, 0)).normalized() : Eigen::Vector4d::Zero();
            Jacobian.col(count * 3 + 1) = alpha.cross3(r);

            alpha = current_bone->dofrz == true ? (current_bone->rotation * Eigen::Vector4d(0, 0, 1, 0)).normalized() : Eigen::Vector4d::Zero();
            Jacobian.col(count * 3 + 2) = alpha.cross3(r);
            
			count++;
            if (current_bone == start_bone) break;
            current_bone = current_bone->parent;
		}

        Eigen::VectorXd deltatheta = step * pseudoInverseLinearSolver(Jacobian, desiredVector);

        // HINT:
        // Change `posture.bone_rotation` based on deltatheta
        count = 0;
        current_bone = end_bone;
		while (current_bone != NULL)
		{
            posture.bone_rotations[current_bone->idx] += Eigen::Vector4d(deltatheta[count * 3 + 0], deltatheta[count * 3 + 1], deltatheta[count * 3 + 2], 0);

			count++;
            if (current_bone == start_bone) break;
            current_bone = current_bone->parent;
		}
    }
    // TODO (Bonus)
    // Return IK is stable?
    // i.e. not swinging its hand in air
    if (past_error.norm() > (target_pos - end_bone->end_position).norm()) return false;
    return true;
}
}  // namespace kinematics
