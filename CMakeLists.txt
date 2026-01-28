cmake_minimum_required(VERSION 3.20)

project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

add_executable(rearrangement src/main.cpp)
target_compile_options(rearrangement PRIVATE -O3 -std=c++20 -static)
install(TARGETS rearrangement RUNTIME DESTINATION ${SKBUILD_SCRIPTS_DIR})
install(PROGRAMS src/correct_micro_homology.awk src/naapam_parse.sh src/naapam_align.sh DESTINATION ${SKBUILD_SCRIPTS_DIR})
