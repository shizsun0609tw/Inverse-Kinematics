#include "graphics/box.h"

#include "glad/gl.h"

#include "graphics/shader.h"

namespace graphics {

bool Box::isInitialized = false;
std::weak_ptr<Buffer<1, GL_ARRAY_BUFFER>> Box::vbo_weak;
std::weak_ptr<Buffer<1, GL_ELEMENT_ARRAY_BUFFER>> Box::ebo_weak;

Box::Box() noexcept {
    if (!(vbo = Box::vbo_weak.lock())) {
        vbo = std::make_shared<Buffer<1, GL_ARRAY_BUFFER>>();
        vbo_weak = vbo;
    }
    if (!(ebo = Box::ebo_weak.lock())) {
        ebo = std::make_shared<Buffer<1, GL_ELEMENT_ARRAY_BUFFER>>();
        ebo_weak = ebo;
    }
    this->initialize();
}

Box::Box(const Box& other) noexcept : Rigidbody(other), vbo(other.vbo), ebo(other.ebo) { this->initialize(); }

Box::Box(Box&& other) noexcept
    : Rigidbody(std::forward<Box>(other)), vbo(std::move(other.vbo)), ebo(std::move(other.ebo)) {}

Box& Box::operator=(const Box& other) noexcept {
    if (this != &other) {
        Rigidbody::operator=(other);
        vbo = other.vbo;
        ebo = other.ebo;
        this->initialize();
    }
    return *this;
}

Box& Box::operator=(Box&& other) noexcept {
    if (this != &other) {
        Rigidbody::operator=(std::forward<Box>(other));
        vbo = std::move(other.vbo);
        ebo = std::move(other.ebo);
    }
    return *this;
}

void Box::initialize() {
    generateVertices();
    glBindVertexArray(vao);
    vbo->bind();
    // Vertex
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), nullptr);
    ebo->bind();
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    // Unbind all
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void Box::bindVBO() const { vbo->bind(); }

void Box::generateVertices() {
    if (!isInitialized) {
        GLfloat vertices[] = {-1.0f, -1.0f, -1.0f, 1.0f, -1.0f, -1.0f, 1.0f, 1.0f, -1.0f, -1.0f, 1.0f, -1.0f,
                              -1.0f, -1.0f, 1.0f,  1.0f, -1.0f, 1.0f,  1.0f, 1.0f, 1.0f,  -1.0f, 1.0f, 1.0f};
        GLuint indices[] = {0, 1, 3, 3, 1, 2, 1, 5, 2, 2, 5, 6, 5, 4, 6, 6, 4, 7,
                            4, 0, 7, 7, 0, 3, 3, 2, 7, 7, 2, 6, 4, 5, 0, 0, 5, 1};
        Box::vbo->bind();
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        Box::ebo->bind();
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        Box::isInitialized = true;
    }
}

void Box::render(Program* shaderProgram) {
    glDepthFunc(GL_LEQUAL);  // Put on the background
    glBindVertexArray(vao);
    shaderProgram->setUniform("skybox", texture->getIndex());
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, nullptr);
    glBindVertexArray(0);
    glDepthFunc(GL_LESS);
}
}  // namespace graphics
