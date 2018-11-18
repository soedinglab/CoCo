cmake_minimum_required(VERSION 3.4)
project(COCO)

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")

include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
add_definitions("-DGIT_SHA1=${GIT_SHA1}")

#set(COCO_SOURCE_DIR src)
set(coco_source_files
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
    )
add_subdirectory(src)
include_directories(kseq)
#include_directories(src)

list (APPEND KSIZE_DEFAULT_LIST 64)

add_subdirectory(lib/gatb-core/gatb-core)
# We add the compilation options for the library
add_definitions (${gatb-core-flags})
include_directories (${gatb-core-includes})
# We add the path for extra libraries
link_directories (${gatb-core-extra-libraries-path})


#add_library(${coco_header_files}
#           ${coco_source_files})

#file (GLOB_RECURSE  SOURCE_FILES  ${COCO_SOURCE_DIR}/*)
add_executable(coco ${coco_source_files} ${coco_header_files})
target_link_libraries (coco ${gatb-core-libraries})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
