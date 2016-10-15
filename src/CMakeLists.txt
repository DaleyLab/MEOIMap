# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

#   You may need to set these for your system
#SET(CMAKE_MODULE_PATH "/...pathTo.../SeqAn/seqan/util/cmake")
#SET(SEQAN_INCLUDE_PATH "/...pathTo.../SeqAn/seqan/include")
#SET(CMAKE_CXX_COMPILER "/...pathTo.../g++")

#Specify use of PATH current c and c++ compilers
find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

# Name of project and that it is C++ only.
project (MEOIMap CXX)

# Add srcs
file( GLOB SRC *.cpp *.h)

# Set flags for compilation
# Optimize and profile
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -pg -std=c++11")

# required for the find seqan module to get a checkpoint
# obv make sure you have actually set this above 
SET(CXX11_FOUND "yes")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB)

# Load the SeqAn module and fail if not found.
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Add executable and link against SeqAn dependencies.
add_executable (MEOIMap ${SRC})
target_link_libraries (MEOIMap ${SEQAN_LIBRARIES})
