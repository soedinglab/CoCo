cmake_minimum_required(VERSION 3.4)
project(COCO)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")

set(COCO_SOURCE_DIR src)

include_directories(kseq)

list (APPEND KSIZE_DEFAULT_LIST 64)

add_subdirectory(lib/gatb-core/gatb-core)
# We add the compilation options for the library
add_definitions (${gatb-core-flags})
include_directories (${gatb-core-includes})
# We add the path for extra libraries
link_directories (${gatb-core-extra-libraries-path})

file (GLOB_RECURSE  SOURCE_FILES  ${COCO_SOURCE_DIR}/*)
add_executable(coco main.cpp)
target_link_libraries (coco ${gatb-core-libraries})
