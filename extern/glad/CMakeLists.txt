cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
project(glad
    VERSION         2
    DESCRIPTION     "Multi-Language Vulkan/GL/GLES/EGL/GLX/WGL Loader-Generator based on the official specs."
    HOMEPAGE_URL    "https://github.com/Dav1dde/glad/tree/glad2"
    LANGUAGES       C
)
add_library(glad INTERFACE)
target_include_directories(glad INTERFACE include)
if (MSVC)
    target_compile_options(glad INTERFACE "/wd4551")
endif()
