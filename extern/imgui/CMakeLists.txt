cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
project(imgui
    VERSION         1.82
    DESCRIPTION     "Bloat-free Graphical User interface for C++ with minimal dependencies"
    HOMEPAGE_URL    "https://github.com/ocornut/imgui"
    LANGUAGES       CXX
)

add_library(imgui STATIC
    src/imgui.cpp
    src/imgui_draw.cpp
    src/imgui_impl_glfw.cpp
    src/imgui_impl_opengl3.cpp
    src/imgui_tables.cpp
    src/imgui_widgets.cpp
)

add_dependencies(imgui glfw glad)
target_link_libraries(imgui
    PRIVATE glfw
    PRIVATE glad
)
# public and internal headers
target_include_directories(imgui
    PUBLIC include
    PRIVATE src
)
# Use c++17
target_compile_features(imgui PRIVATE cxx_std_17)
set_target_properties(imgui PROPERTIES CMAKE_CXX_EXTENSIONS OFF)
target_compile_definitions(imgui PUBLIC IMGUI_IMPL_OPENGL_LOADER_GLAD2)
