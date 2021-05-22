cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
project(eigen
    VERSION         3.3.9
    DESCRIPTION     "Eigen is a C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms."
    HOMEPAGE_URL    "https://eigen.tuxfamily.org"
    LANGUAGES       CXX
)

add_library(eigen INTERFACE)
target_include_directories(eigen INTERFACE include)
target_compile_definitions(eigen INTERFACE EIGEN_MPL2_ONLY)
target_compile_definitions(eigen INTERFACE EIGEN_NO_DEBUG)
target_compile_definitions(eigen INTERFACE EIGEN_DONT_PARALLELIZE)
