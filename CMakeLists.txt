cmake_minimum_required(VERSION 3.20)
set(CMAKE_C_COMPILER clang)
project(gridgen C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Fall back on debug build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# Header files
include_directories(include/)

# Fix for pkg-config
set(ENV{PKG_CONFIG_PATH} "/usr/local/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")

# Dependencies
find_package(GSL REQUIRED)
find_package(PkgConfig REQUIRED)
find_package(LibXml2)
pkg_check_modules(CRITERION criterion)

# Source files
set(LIB_SOURCES
    source/errors.c
    source/config.c
    source/aux_functions.c
    source/file_operations.c
    source/geometry.c
    source/jacobian_ders.c
    source/newton.c
    source/grid_solvers.c
)

# Place executable in bin/
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Main executable
add_executable(gridexec ${LIB_SOURCES} source/main.c)
target_link_libraries(gridexec
    GSL::gsl GSL::gslcblas m
)

# Test executable
if(CRITERION_FOUND)
    file(GLOB TEST_SOURCES tests/*.c)

    add_executable(tests_runner ${LIB_SOURCES} ${TEST_SOURCES})
    target_include_directories(tests_runner PRIVATE ${CRITERION_INCLUDE_DIRS})
    target_link_libraries(tests_runner
        ${CRITERION_LIBRARIES} pthread
        GSL::gsl GSL::gslcblas m
    )
    target_link_options(tests_runner PRIVATE -no-pie)

    if(LibXml2_FOUND)
        target_include_directories(tests_runner PRIVATE ${LIBXML2_INCLUDE_DIR})
        target_link_libraries(tests_runner ${LIBXML2_LIBRARIES})
        target_compile_definitions(tests_runner PRIVATE HAVE_LIBXML2)
    endif()
else()
    message(WARNING "Criterion not found - tests_runner target skipped.")
endif()
