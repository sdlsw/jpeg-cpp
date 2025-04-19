cmake_minimum_required(VERSION 3.30)

# Enable support of `import std`
# See https://www.kitware.com/import-std-in-cmake-3-30/
set(
    CMAKE_EXPERIMENTAL_CXX_IMPORT_STD
    "0e5b6991-d74f-4b3d-a41c-cf096e0b2508"
)

project(jpeg-impl VERSION 1.0 LANGUAGES CXX)

# Compiler options.
add_library(compiler-flags INTERFACE)
target_compile_features(compiler-flags INTERFACE cxx_std_23)
set(msvc "$<COMPILE_LANG_AND_ID:CXX,MSVC>")
target_compile_options(compiler-flags INTERFACE "$<${msvc}:-W3>")

add_subdirectory(msg)
add_subdirectory(jpeg)

add_executable(main main.cpp)
set_target_properties(main PROPERTIES OUTPUT_NAME "jpeg")
target_link_libraries(main PUBLIC msg jpeg compiler-flags)
