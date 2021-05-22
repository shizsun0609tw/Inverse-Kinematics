#include "acclaim/skeleton.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>

#include "util/helper.h"

namespace acclaim {

Skeleton::Skeleton(const util::fs::path &file_name, const double _scale) noexcept : scale(_scale) {
    bones.reserve(64);
    // Initializaton of root bone
    bones[0].name = std::string("root");
    bones[0].idx = root_idx();
    bones[0].parent = nullptr;
    bones[0].sibling = nullptr;
    bones[0].child = nullptr;
    bones[0].dir.setZero();
    bones[0].axis.setZero();
    bones[0].length = 0;
    bones[0].dof = 6;
    bones[0].dofrx = true;
    bones[0].dofry = true;
    bones[0].dofrz = true;
    bones[0].doftx = true;
    bones[0].dofty = true;
    bones[0].doftz = true;
    // build hierarchy and read in each bone's DOF information
    readASFFile(file_name);
    computeLocalDirection();
    computeLocalRotation();

    setBoneGraphics();
}

Skeleton::Skeleton(const Skeleton &other) noexcept
    : scale(other.scale), movableBones(other.movableBones), bones(other.bones), bone_graphics(other.bone_graphics) {
    for (std::size_t i = 0; i < bones.size(); ++i) {
        if (bones[i].parent != nullptr) {
            bones[i].parent = &bones[other.bones[i].parent->idx];
        }
        if (bones[i].child != nullptr) {
            bones[i].child = &bones[other.bones[i].child->idx];
        }
        if (bones[i].sibling != nullptr) {
            bones[i].sibling = &bones[other.bones[i].sibling->idx];
        }
    }
}

Skeleton::Skeleton(Skeleton &&other) noexcept
    : scale(other.scale),
      movableBones(other.movableBones),
      bones(std::move(other.bones)),
      bone_graphics(std::move(other.bone_graphics)) {}

Skeleton &Skeleton::operator=(const Skeleton &other) noexcept {
    if (this != &other) {
        scale = other.scale;
        movableBones = other.movableBones;
        bones = other.bones;
        // We need to reset all pointer in bones
        for (std::size_t i = 0; i < bones.size(); ++i) {
            if (bones[i].parent != nullptr) {
                bones[i].parent = &bones[other.bones[i].parent->idx];
            }
            if (bones[i].child != nullptr) {
                bones[i].child = &bones[other.bones[i].child->idx];
            }
            if (bones[i].sibling != nullptr) {
                bones[i].sibling = &bones[other.bones[i].sibling->idx];
            }
        }
        bone_graphics = other.bone_graphics;
    }
    return *this;
}

Skeleton &Skeleton::operator=(Skeleton &&other) noexcept {
    if (this != &other) {
        scale = other.scale;
        movableBones = other.movableBones;
        bones = std::move(other.bones);
        bone_graphics = std::move(other.bone_graphics);
    }
    return *this;
}

double Skeleton::getScale() const { return scale; }

int Skeleton::getBoneNum() const { return static_cast<int>(bones.size()); }

int Skeleton::getMovableBoneNum() const { return movableBones; }

Bone *Skeleton::getBonePointer(const std::string &name) {
    for (size_t i = 0; i < bones.size(); ++i) {
        if (name == bones[i].name) {
            return &bones[i];
        }
    }
    return nullptr;
}

Bone *Skeleton::getBonePointer(const int bone_idx) { return &bones[bone_idx]; }

void Skeleton::setBoneColor(const Eigen::Vector4f &boneColor) const {
    for (size_t i = 0; i < bones.size(); ++i) bone_graphics[i].setTexture(boneColor);
}

void Skeleton::setModelMatrices() const {
    for (size_t i = 0; i < bone_graphics.size(); ++i) {
        auto &&bone = bones[i];
        Eigen::Vector4d trans = 0.5 * (bone.start_position + bone.end_position);
        Eigen::Affine3d model = bone.global_facing;
        model.prerotate(bone.rotation.rotation()).pretranslate(trans.head<3>());
        bone_graphics[i].setModelMatrix(model.cast<float>());
    }
}

void Skeleton::render(graphics::Program *program) const {
    for (size_t i = 0; i < bone_graphics.size(); ++i) {
        bone_graphics[i].render(program);
    }
}

bool Skeleton::readASFFile(const util::fs::path &file_name) {
    std::ifstream input_stream(file_name);
    if (!input_stream) {
        std::cerr << "Failed to open " << file_name << std::endl;
        return false;
    }
    // ignore header information
    std::string str, keyword;
    while (true) {
        std::getline(input_stream, str);
        if (str.compare(0, 9, ":bonedata") == 0) {
            break;
        }
    }
    // Ignore begin
    input_stream.ignore(1024, '\n');
    bool done = false;
    while (!done) {
        auto &&current_bone = bones.emplace_back();
        while (true) {
            input_stream >> keyword;
            if (keyword == "end") {
                break;
            }
            // finish read bone data, start setup bone hierarchy
            if (keyword == ":hierarchy") {
                done = true;
                bones.pop_back();
                break;
            }
            // id of bone
            if (keyword == "id") {
                input_stream >> current_bone.idx;
                continue;
            }
            // name of the bone
            if (keyword == "name") {
                input_stream >> current_bone.name;
                continue;
            }
            // this line describes the bone's direction vector in global coordinate
            // it will later be converted to local coorinate system
            if (keyword == "direction") {
                input_stream >> current_bone.dir[0] >> current_bone.dir[1] >> current_bone.dir[2];
                continue;
            }
            // length of the bone
            if (keyword == "length") {
                input_stream >> current_bone.length;
                current_bone.length *= scale;
                continue;
            }
            // this line describes the orientation of bone's local coordinate
            // system relative to the world coordinate system
            if (keyword == "axis") {
                input_stream >> current_bone.axis[0] >> current_bone.axis[1] >> current_bone.axis[2];
                continue;
            }
            // this line describes the bone's dof
            if (keyword == "dof") {
                ++movableBones;
                std::string token;
                current_bone.dof = 0;
                std::getline(input_stream, str);
                std::istringstream iss(str);
                while (iss >> token) {
                    if (token.compare(0, 2, "rx") == 0) {
                        current_bone.dofrx = true;
                        ++current_bone.dof;
                    } else if (token.compare(0, 2, "ry") == 0) {
                        current_bone.dofry = true;
                        ++current_bone.dof;
                    } else if (token.compare(0, 2, "rz") == 0) {
                        current_bone.dofrz = true;
                        ++current_bone.dof;
                    } else if (token.compare(0, 2, "tx") == 0) {
                        current_bone.doftx = true;
                        ++current_bone.dof;
                    } else if (token.compare(0, 2, "ty") == 0) {
                        current_bone.dofty = true;
                        ++current_bone.dof;
                    } else if (token.compare(0, 2, "tz") == 0) {
                        current_bone.doftz = true;
                        ++current_bone.dof;
                    } else {
                        std::cerr << "Unknown token: " << token << std::endl;
                    }
                }
            }
        }
    }
    // skip "begin" line
    input_stream.ignore(1024, '\n');
    input_stream.ignore(1024, '\n');
    // Assign parent/child relationship to the bones
    while (true) {
        // read next line
        std::getline(input_stream, str);
        std::istringstream iss(str);
        iss >> keyword;
        // check if we are done
        if (keyword == "end") break;
        // parse this line, it contains parent followed by children
        Bone *parent = this->getBonePointer(keyword);
        while (iss >> keyword) {
            this->setBoneHierarchy(parent, getBonePointer(keyword));
        }
    }
    std::cout << bones.size() << " bones in " << file_name.string() << " are read" << std::endl;
    input_stream.close();
    return true;
}

void Skeleton::setBoneHierarchy(Bone *parent, Bone *child) {
    if (parent == nullptr) {
        std::cerr << "inbord bone is undefined" << std::endl;
        return;
    }
    child->parent = parent;
    if (parent->child == nullptr) {
        parent->child = child;
    } else {
        Bone *current = parent->child;
        while (current->sibling != nullptr) current = current->sibling;
        current->sibling = child;
    }
}

void Skeleton::computeLocalDirection() {
    for (size_t i = 1; i < bones.size(); ++i) {
        // Transform dir vector into local coordinate system
        bones[i].dir = Eigen::Affine3d(util::rotateDegreeXYZ(-bones[i].axis)) * bones[i].dir;
    }
}

void Skeleton::computeLocalRotation() {
    Eigen::Quaterniond rootRotation = util::rotateDegreeZYX(bones[0].axis);
    bones[0].rot_parent_current = rootRotation;
    // Compute rot_parent_current for all other bones
    for (size_t i = 0; i < bones.size(); ++i) {
        if (bones[i].child != nullptr) {
            // compute child
            bones[i].child->rot_parent_current =
                util::rotateDegreeXYZ(-bones[i].axis) * util::rotateDegreeZYX(bones[i].child->axis);
            // compute siblings
            for (Bone *current = bones[i].child->sibling; current != nullptr; current = current->sibling) {
                current->rot_parent_current =
                    util::rotateDegreeXYZ(-bones[i].axis) * util::rotateDegreeZYX(current->axis);
            }
        }
    }
}

void Skeleton::setBoneGraphics() {
    bone_graphics.resize(bones.size());
    for (size_t i = 0; i < bones.size(); ++i) {
        bone_graphics[i].setTexture(Eigen::Vector4f(0.6f, 0.6f, 0.0f, 1.0f));
        auto &&bone = bones[i];
        Eigen::Vector4d rotaion_axis = Eigen::Vector4d::UnitZ().cross3(bone.dir);
        double dot_val = Eigen::Vector4d::UnitZ().dot(bone.dir);
        double cross_val = rotaion_axis.norm();
        rotaion_axis.normalize();
        double theta = atan2(cross_val, dot_val);
        bone.global_facing = Eigen::AngleAxisd(theta, rotaion_axis.head<3>());
        bone.global_facing.scale(Eigen::Vector3d(1.0, 1.0, bone.length));
    }
}

}  // namespace acclaim
