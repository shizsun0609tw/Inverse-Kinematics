#include "graphics/rigidbody.h"

#include <utility>

namespace graphics {
Rigidbody::Rigidbody() noexcept { glGenVertexArrays(1, &vao); }

Rigidbody::Rigidbody(const Rigidbody& other) noexcept
    : vao(0),
      modelMatrix(other.modelMatrix),
      inverseTransposeModel(other.inverseTransposeModel),
      baseColor(other.baseColor),
      texture(other.texture) {
    glGenVertexArrays(1, &vao);
}

Rigidbody::Rigidbody(Rigidbody&& other) noexcept
    : vao(other.vao),
      modelMatrix(std::move(other.modelMatrix)),
      inverseTransposeModel(std::move(other.inverseTransposeModel)),
      baseColor(std::move(other.baseColor)),
      texture(std::move(other.texture)) {
    other.vao = 0;
}

Rigidbody::~Rigidbody() { glDeleteVertexArrays(1, &vao); }

Rigidbody& Rigidbody::operator=(const Rigidbody& other) noexcept {
    if (this != &other) {
        modelMatrix = other.modelMatrix;
        inverseTransposeModel = other.inverseTransposeModel;
        baseColor = other.baseColor;
        texture = other.texture;
    }
    return *this;
}
Rigidbody& Rigidbody::operator=(Rigidbody&& other) noexcept {
    if (this != &other) {
        vao = other.vao;
        other.vao = 0;
        modelMatrix = std::move(other.modelMatrix);
        inverseTransposeModel = std::move(other.inverseTransposeModel);
        baseColor = std::move(other.baseColor);
        texture = std::move(other.texture);
    }
    return *this;
}

void Rigidbody::setModelMatrix(const Eigen::Affine3f& _modelMatrix) {
    modelMatrix = _modelMatrix;
    inverseTransposeModel = modelMatrix.inverse().matrix().transpose();
}

void Rigidbody::setTexture(const std::shared_ptr<TextureBase>& _texture) { texture = _texture; }

void Rigidbody::setTexture(const Eigen::Vector4f& color) { baseColor = color; }

void Rigidbody::initialize() {
    generateVertices();
    glBindVertexArray(vao);
    bindVBO();
    // Vertex
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), nullptr);
    // normal
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), reinterpret_cast<void*>(3 * sizeof(GLfloat)));
    // Texture
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), reinterpret_cast<void*>(6 * sizeof(GLfloat)));
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}
}  // namespace graphics
