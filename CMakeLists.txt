cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(InverseKinematics LANGUAGES CXX)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin/$<0:>)
# Set to Release by default
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
if (NOT (CMAKE_BUILD_TYPE OR CMAKE_CONFIGURATION_TYPES))
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release" "RelWithDebInfo" "MinSizeRel" "Debug")
endif ()
# Detect some compiler flags
include(CheckCXXCompilerFlag)
include(CheckIPOSupported)
if (NOT COMPILER_FLAG_TEST_COMPLETE)
    check_ipo_supported(RESULT COMPILER_SUPPORT_IPO)
    check_cxx_compiler_flag("-march=native" COMPILER_SUPPORT_MARCH_NATIVE)
    check_cxx_compiler_flag("-xHost" COMPILER_SUPPORT_xHOST)
    check_cxx_compiler_flag("/QxHost" COMPILER_SUPPORT_QxHOST)
    set(COMPILER_FLAG_TEST_COMPLETE TRUE)
endif()
# Link time optimization
if (MSVC)
	add_compile_options("$<$<CONFIG:Release>:/GL>")
	add_link_options("$<$<CONFIG:Release>:/LTCG:incremental>")
elseif(COMPILER_SUPPORT_IPO)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
# Visual Studio need this for multithread compile
if (MSVC)
    add_compile_options("/MP")
endif()
# SIMD support
if (COMPILER_SUPPORT_MARCH_NATIVE)
    add_compile_options("-march=native")
elseif(COMPILER_SUPPORT_xHOST)
    add_compile_options("-xHost")
elseif(COMPILER_SUPPORT_QxHOST)
    add_compile_options("/QxHost")
elseif(MSVC)
	add_compile_options("/arch:AVX2")
endif()
# Softbody simulation part
add_executable(InverseKinematics
    ${CMAKE_CURRENT_SOURCE_DIR}/src/acclaim/motion.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/acclaim/posture.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/acclaim/skeleton.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/graphics/box.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/graphics/camera.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/graphics/cylinder.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/graphics/default_camera.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/graphics/free_camera.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/graphics/plane.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/graphics/rigidbody.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/graphics/shader.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/graphics/sphere.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/graphics/texture.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/simulation/ball.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/simulation/kinematics.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/util/filesystem.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/util/helper.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/InverseKinematics/main.cpp
)
# Base include files
target_include_directories(InverseKinematics PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
# Use C++17
target_compile_features(InverseKinematics PRIVATE cxx_std_17)
target_compile_options(InverseKinematics PRIVATE "-Wall")
target_compile_options(InverseKinematics PRIVATE "-Wextra")
target_compile_options(InverseKinematics PRIVATE "-Wpedantic")
# Use std C++17 not GNU C++17
set_target_properties(InverseKinematics PROPERTIES CMAKE_CXX_EXTENSIONS OFF)
target_compile_definitions(InverseKinematics PRIVATE GLFW_INCLUDE_NONE)
# Add third-party libraries
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extern)
# Link those third-party libraries
target_link_libraries(InverseKinematics
    PRIVATE eigen
    PRIVATE glad
    PRIVATE glfw
    PRIVATE imgui
    PRIVATE stb
    PRIVATE FK
)
