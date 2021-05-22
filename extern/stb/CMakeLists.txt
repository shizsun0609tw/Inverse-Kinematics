cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
project(stb
    DESCRIPTION     "stb single-file public domain libraries for C/C++"
    HOMEPAGE_URL    "https://github.com/nothings/stb"
    LANGUAGES       C
)
add_library(stb INTERFACE)
target_include_directories(stb INTERFACE include)
