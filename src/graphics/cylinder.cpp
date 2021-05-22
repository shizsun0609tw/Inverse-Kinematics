#include "graphics/cylinder.h"

#include <vector>

#include "glad/gl.h"

#include "graphics/configs.h"
#include "graphics/shader.h"
#include "util/helper.h"

namespace graphics {

namespace {
GLsizei nIndices = 0;
}

bool Cylinder::isInitialized = false;
std::weak_ptr<Buffer<1, GL_ARRAY_BUFFER>> Cylinder::vbo_weak;
std::weak_ptr<Buffer<2, GL_ELEMENT_ARRAY_BUFFER>> Cylinder::ebo_weak;

Cylinder::Cylinder() noexcept {
    if (!(vbo = Cylinder::vbo_weak.lock())) {
        vbo = std::make_shared<Buffer<1, GL_ARRAY_BUFFER>>();
        vbo_weak = vbo;
    }
    if (!(ebo = Cylinder::ebo_weak.lock())) {
        ebo = std::make_shared<Buffer<2, GL_ELEMENT_ARRAY_BUFFER>>();
        ebo_weak = ebo;
    }
    initialize();
    glBindVertexArray(vao);
    ebo->bind();
    glBindVertexArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

Cylinder::Cylinder(const Cylinder& other) noexcept : Rigidbody(other), vbo(other.vbo), ebo(other.ebo) {
    initialize();
    glBindVertexArray(vao);
    ebo->bind();
    glBindVertexArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

Cylinder::Cylinder(Cylinder&& other) noexcept
    : Rigidbody(std::forward<Cylinder>(other)), vbo(std::move(other.vbo)), ebo(std::move(other.ebo)) {}

Cylinder& Cylinder::operator=(const Cylinder& other) noexcept {
    if (this != &other) {
        Rigidbody::operator=(other);
        vbo = other.vbo;
        ebo = other.ebo;
        initialize();
        glBindVertexArray(vao);
        ebo->bind();
        glBindVertexArray(0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    }
    return *this;
}

Cylinder& Cylinder::operator=(Cylinder&& other) noexcept {
    if (this != &other) {
        Rigidbody::operator=(std::forward<Cylinder>(other));
        vbo = std::move(other.vbo);
        ebo = std::move(other.ebo);
    }
    return *this;
}

void Cylinder::bindVBO() const { vbo->bind(); }

void Cylinder::generateVertices() {
    // http://www.songho.ca/opengl/gl_cylinder.html#cylinderconst
    if (!Cylinder::isInitialized) {
        std::vector<GLfloat> vertices;
        std::vector<GLuint> indices;

        vertices.reserve(2 * 8 * (2 * g_CylinderSectors + 1) + 2);

        float sectorStep = static_cast<float>(2.0 * util::PI / g_CylinderSectors);
        float sectorAngle = 0;  // radian
        std::vector<float> unitCircle((g_CylinderSectors + 1) * 3);

        std::vector<float> unitCircleVertices;
        int currentPos = -1;
        for (int i = 0; i <= g_CylinderSectors; ++i) {
            unitCircle[++currentPos] = cosf(sectorAngle);  // x
            unitCircle[++currentPos] = sinf(sectorAngle);  // y
            unitCircle[++currentPos] = 0;                  // z
            sectorAngle += sectorStep;
        }
        // put side vertices to arrays
        for (int i = 0; i < 2; ++i) {
            float h = -g_CylinderHeight / 2.0f + i * g_CylinderHeight;  // z value; -h/2 to h/2
            float t = 1.0f - i;                                         // vertical tex coord; 1 to 0

            for (int j = 0, k = 0; j <= g_CylinderSectors; ++j, k += 3) {
                float ux = unitCircle[k];
                float uy = unitCircle[k + 1];
                float uz = unitCircle[k + 2];
                float s = static_cast<float>(j) / g_CylinderSectors;
                vertices.insert(vertices.end(), {ux * g_CylinderRadius, uy * g_CylinderRadius, h, ux, uy, uz, s, t});
            }
        }
        // the starting index for the base/top surface
        // NOTE: it is used for generating indices later
        int baseCenterIndex = static_cast<int>(vertices.size()) / 8;
        int topCenterIndex = baseCenterIndex + g_CylinderSectors + 1;  // include center vertex

        // put base and top vertices to arrays
        for (int i = 0; i < 2; ++i) {
            float h = -g_CylinderHeight / 2.0f + i * g_CylinderHeight;  // z value; -h/2 to h/2
            float nz = -1.0f + 2.0f * i;                                // z value of normal; -1 to 1

            // center point
            vertices.insert(vertices.end(), {0, 0, h, 0, 0, nz, 0.5f, 0.5f});

            for (int j = 0, k = 0; j < g_CylinderSectors; ++j, k += 3) {
                float ux = unitCircle[k];
                float uy = unitCircle[k + 1];
                float s = -ux * 0.5f + 0.5f;
                float t = -uy * 0.5f + 0.5f;
                vertices.insert(vertices.end(), {ux * g_CylinderRadius, uy * g_CylinderRadius, h, 0, 0, nz, s, t});
            }
        }

        int k1 = 0;                      // 1st vertex index at base
        int k2 = g_CylinderSectors + 1;  // 1st vertex index at top

        // indices for the side surface
        for (int i = 0; i < g_CylinderSectors; ++i, ++k1, ++k2) {
            // 2 triangles per sector
            // k1 => k1+1 => k2
            indices.push_back(k1);
            indices.push_back(k1 + 1);
            indices.push_back(k2);

            // k2 => k1+1 => k2+1
            indices.push_back(k2);
            indices.push_back(k1 + 1);
            indices.push_back(k2 + 1);
        }

        // indices for the base surface
        // NOTE: baseCenterIndex and topCenterIndices are pre-computed during vertex generation
        //      please see the previous code snippet
        for (int i = 0, k = baseCenterIndex + 1; i < g_CylinderSectors; ++i, ++k) {
            if (i < g_CylinderSectors - 1) {
                indices.push_back(baseCenterIndex);
                indices.push_back(k + 1);
                indices.push_back(k);
            } else {
                indices.push_back(baseCenterIndex);
                indices.push_back(baseCenterIndex + 1);
                indices.push_back(k);
            }
        }

        // indices for the top surface
        for (int i = 0, k = topCenterIndex + 1; i < g_CylinderSectors; ++i, ++k) {
            if (i < g_CylinderSectors - 1) {
                indices.push_back(topCenterIndex);
                indices.push_back(k);
                indices.push_back(k + 1);
            } else {
                indices.push_back(topCenterIndex);
                indices.push_back(k);
                indices.push_back(topCenterIndex + 1);
            }
        }
        vbo->bind();
        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(GLfloat), vertices.data(), GL_STATIC_DRAW);
        ebo->bind();
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), indices.data(), GL_STATIC_DRAW);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        nIndices = static_cast<GLsizei>(indices.size());
        Cylinder::isInitialized = true;
    }
}

void Cylinder::render(Program* shaderProgram) {
    if (texture) {
        shaderProgram->setUniform("useTexture", 1);
        shaderProgram->setUniform("diffuseTexture", texture->getIndex());
    } else {
        shaderProgram->setUniform("useTexture", 0);
        shaderProgram->setUniform("baseColor", baseColor);
    }
    shaderProgram->setUniform("model", modelMatrix);
    shaderProgram->setUniform("invtransmodel", inverseTransposeModel);
    glBindVertexArray(vao);
    glDrawElements(GL_TRIANGLES, nIndices, GL_UNSIGNED_INT, nullptr);
    glBindVertexArray(0);
}
}  // namespace graphics
