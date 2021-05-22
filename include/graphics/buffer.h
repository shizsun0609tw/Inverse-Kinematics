#pragma once
#include <array>

#include "glad/gl.h"

namespace graphics {

template <int size_, int buffer_type>
class Buffer final {
    static_assert(size_ > 0, "size must be positive");

 public:
    // Allocates the buffer
    Buffer() { glGenBuffers(size_, buffer.data()); }
    // You should not copy this class
    Buffer(const Buffer &) = delete;
    Buffer &operator=(const Buffer &) = delete;
    // You can move this class
    Buffer(Buffer &&) = default;
    Buffer &operator=(Buffer &&) = default;
    // Default bind to 0
    void bind(int num = 0) { glBindBuffer(buffer_type, buffer[num]); }
    // Deallocates the buffer
    ~Buffer() { glDeleteBuffers(size_, buffer.data()); }

 private:
    std::array<GLuint, size_> buffer{};
};
}  // namespace graphics
