#include "graphics/sphere.h"
#include <vector>

#include "glad/gl.h"

#include "graphics/configs.h"
#include "graphics/shader.h"
#include "util/helper.h"

namespace graphics {
namespace {
constexpr GLsizei nFilled = 6 * g_SphereStacks * g_SphereSectors - 6 * g_SphereStacks;
constexpr GLsizei nWire = 4 * g_SphereStacks * g_SphereSectors - 2 * g_SphereStacks;
}  // namespace

bool Sphere::isInitialized = false;
std::weak_ptr<Buffer<1, GL_ARRAY_BUFFER>> Sphere::vbo_weak;
std::weak_ptr<Buffer<2, GL_ELEMENT_ARRAY_BUFFER>> Sphere::ebo_weak;

Sphere::Sphere(RenderMode m) noexcept : mode(m) {
    if (!(vbo = Sphere::vbo_weak.lock())) {
        vbo = std::make_shared<Buffer<1, GL_ARRAY_BUFFER>>();
        vbo_weak = vbo;
    }
    if (!(ebo = Sphere::ebo_weak.lock())) {
        ebo = std::make_shared<Buffer<2, GL_ELEMENT_ARRAY_BUFFER>>();
        ebo_weak = ebo;
    }
    initialize();
    glBindVertexArray(vao);
    ebo->bind(static_cast<int>(mode));
    glBindVertexArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

Sphere::Sphere(const Sphere& other) noexcept : Rigidbody(other), mode(other.mode), vbo(other.vbo), ebo(other.ebo) {
    initialize();
}

Sphere::Sphere(Sphere&& other) noexcept
    : Rigidbody(std::forward<Sphere>(other)), mode(other.mode), vbo(std::move(other.vbo)), ebo(std::move(other.ebo)) {}

Sphere& Sphere::operator=(const Sphere& other) noexcept {
    if (this != &other) {
        Rigidbody::operator=(other);
        mode = other.mode;
        vbo = other.vbo;
        ebo = other.ebo;
        initialize();
    }
    return *this;
}

Sphere& Sphere::operator=(Sphere&& other) noexcept {
    if (this != &other) {
        Rigidbody::operator=(std::forward<Sphere>(other));
        mode = other.mode;
        vbo = std::move(other.vbo);
        ebo = std::move(other.ebo);
    }
    return *this;
}

void Sphere::generateVertices() {
    // http://www.songho.ca/opengl/gl_sphere.html#sphere
    if (!Sphere::isInitialized) {
        std::vector<GLfloat> vertices;
        std::vector<GLuint> fillIndices, wireIndices;
        vertices.reserve(8 * (g_SphereStacks + 1) * (g_SphereSectors + 1));
        fillIndices.reserve(6 * g_SphereStacks * g_SphereSectors);
        wireIndices.reserve(4 * g_SphereStacks * g_SphereSectors);

        float x, y, z, xy;                                    //  position
        float nx, ny, nz, lengthInv = 1.0f / g_SphereRadius;  //  normal
        float s, t;                                           //  texCoord

        constexpr float sectorStep = static_cast<float>(2.0 * util::PI / g_SphereSectors);
        constexpr float stackStep = static_cast<float>(util::PI / g_SphereStacks);
        float sectorAngle, stackAngle;

        for (int i = 0; i <= g_SphereStacks; ++i) {
            stackAngle = static_cast<float>(util::PI / 2.0 - i * stackStep);  // [pi/2, -pi/2]
            xy = cosf(stackAngle);                                            // r * cos(u)
            z = sinf(stackAngle);                                             // r * sin(u)

            for (int j = 0; j <= g_SphereSectors; ++j) {
                sectorAngle = j * sectorStep;  // [0, 2pi]

                x = xy * cosf(sectorAngle);  // r * cos(u) * cos(v)
                y = xy * sinf(sectorAngle);  // r * cos(u) * sin(v)

                // normalized vertex normal (nx, ny, nz)
                nx = x * lengthInv;
                ny = y * lengthInv;
                nz = z * lengthInv;

                // vertex tex coord (s, t) range between [0, 1]
                s = static_cast<float>(j) / g_SphereSectors;
                t = static_cast<float>(i) / g_SphereStacks;

                vertices.insert(vertices.end(), {x, y, z, nx, ny, nz, s, t});
            }
        }

        unsigned int k1, k2;  // EBO index
        for (int i = 0; i < g_SphereStacks; ++i) {
            k1 = i * (g_SphereSectors + 1);  // beginning of current stack
            k2 = k1 + g_SphereSectors + 1;   // beginning of next stack
            for (int j = 0; j < g_SphereSectors; ++j, ++k1, ++k2) {
                wireIndices.insert(wireIndices.end(), {k1, k2});
                if (i != 0) {
                    wireIndices.insert(wireIndices.end(), {k1, k1 + 1});
                    fillIndices.insert(fillIndices.end(), {k1, k2, k1 + 1});
                }
                // k1+1 => k2 => k2+1
                if (i != (g_SphereStacks - 1)) {
                    fillIndices.insert(fillIndices.end(), {k1 + 1, k2, k2 + 1});
                }
            }
        }
        vbo->bind();
        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(GLfloat), vertices.data(), GL_STATIC_DRAW);
        ebo->bind(static_cast<int>(RenderMode::WIREFRAME));
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, wireIndices.size() * sizeof(GLuint), wireIndices.data(), GL_STATIC_DRAW);
        ebo->bind(static_cast<int>(RenderMode::FILLED));
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, fillIndices.size() * sizeof(GLuint), fillIndices.data(), GL_STATIC_DRAW);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        assert(wireIndices.size() == nWire);
        assert(fillIndices.size() == nFilled);
        Sphere::isInitialized = true;
    }
}
void Sphere::bindVBO() const { vbo->bind(); }

void Sphere::render(Program* shaderProgram) {
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
    if (mode == RenderMode::FILLED) {
        glDrawElements(GL_TRIANGLES, nFilled, GL_UNSIGNED_INT, nullptr);
    } else {
        glDrawElements(GL_LINES, nWire, GL_UNSIGNED_INT, nullptr);
    }
    glBindVertexArray(0);
}
}  // namespace graphics
