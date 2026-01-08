
#cmake_minimum_required(VERSION 3.22)
#project(radt LANGUAGES CXX)

# Set C++ standard
#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable verbose output
set(CMAKE_VERBOSE_MAKEFILE ON)

# Find MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# HYPRE installation path
set(HYPRE_DIR "/home/mike/proj/hypre/src/hypre")
#set(GOOGLETEST_DIR "/home/mike/proj/googletest/build")
set(GOOGLETEST_DIR "/home/mike/proj/googletest")


#set(GOOGLETEST_DIR "/home/mike/proj/radhypre_classkrymgdif/exdifkry3/tests/googletest")
include_directories(${HYPRE_DIR}/include
                    ${GOOGLETEST_DIR}/googletest/include
                    ${GOOGLETEST_DIR}/googlemock/include)
link_directories(${HYPRE_DIR}/lib
                 ${GOOGLETEST_DIR}/build/lib)


add_subdirectory(${GOOGLETEST_DIR} ${CMAKE_BINARY_DIR}/gtest_build)




                 

find_library(HYPRE_LIB NAMES HYPRE PATHS ${HYPRE_DIR}/lib NO_DEFAULT_PATH)

if(NOT HYPRE_LIB)
    message(FATAL_ERROR "Could not find HYPRE library in ${HYPRE_DIR}/lib")
endif()
message(STATUS "Found HYPRE library: ${HYPRE_LIB}")


# Include your own headers
include_directories(include)

# Source files
#set(SOURCES
#    hyptest1.cpp
#    hyptest1_unittest.cpp
#)
enable_testing()
# Create executable
#add_executable(hyptest ${SOURCES})

# Link libraries
#target_link_libraries(hyptest ${MPI_CXX_LIBRARIES} HYPRE gtest gtest_main PRIVATE radt_lib)

# Integration tests
add_executable(integration_tests
    TopHatTest.cpp
    PipeRadiationTest.cpp
    MockPipeRadiationTest.cpp
    PipeTwoMaterials_MockedHeating.cpp
    MarshakWaveTest.cpp
)
target_link_libraries(integration_tests  ${MPI_CXX_LIBRARIES} radt_lib ${HYPRE_LIB} gmock gmock_main gtest gtest_main)
add_test(NAME integration_tests COMMAND integration_tests)

