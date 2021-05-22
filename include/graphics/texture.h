#pragma once
#include <array>

#include "Eigen/Dense"
#include "glad/gl.h"

#include "util/filesystem.h"

namespace graphics {

class TextureBase {
 public:
    TextureBase() noexcept;
    TextureBase(const TextureBase&) = delete;
    TextureBase(TextureBase&&) = default;
    TextureBase& operator=(const TextureBase&) = delete;
    TextureBase& operator=(TextureBase&&) = default;

    GLuint getIndex() const;

 protected:
    virtual ~TextureBase();
    // This is for recording which id is free to assign to a new texture
    static GLuint freeIndex;
    GLuint id, index;
};

class Texture final : public TextureBase {
 public:
    explicit Texture(const char* fileName);
    explicit Texture(util::fs::path filePath);

 private:
    void loadTexture(const char* fileName);
};

class ShadowMapTexture final : public TextureBase {
 public:
    ShadowMapTexture(unsigned int size);
    ~ShadowMapTexture();

    unsigned int getShadowSize() const;

    void bindFrameBuffer() const;
    void unbindFrameBuffer() const;

 private:
    GLuint depthMapFBO;
    unsigned int shadowSize;
};

class CubeTexture final : public TextureBase {
 public:
    explicit CubeTexture(const util::fs::path& name);
    explicit CubeTexture(const std::array<util::fs::path, 6>& names);

 private:
    void loadTexture(const std::array<const char*, 6>& fileName);
};
}  // namespace graphics
