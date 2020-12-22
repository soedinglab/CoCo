cmake_minimum_required(VERSION 3.4)
project(CoCo)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")

include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
add_definitions("-DGIT_SHA1=${GIT_SHA1}")

include_directories(lib/kseq)
include_directories(lib/kmer-cnt)

set(coco_source_files ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp src/mask_permuter.cpp src/mask_permuter.h)

list (APPEND KSIZE_DEFAULT_LIST 64)
SET (GATB_CORE_EXCLUDE_TOOLS     1)
SET (GATB_CORE_EXCLUDE_TESTS     1)
SET (GATB_CORE_EXCLUDE_EXAMPLES  1)

add_subdirectory(lib/gatb-core/gatb-core)
add_subdirectory(src)

add_definitions (${gatb-core-flags})
include_directories (${gatb-core-includes})
link_directories (${gatb-core-extra-libraries-path})
set_target_properties(gatbcore-static PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -w" LINK_FLAGS "${CMAKE_CXX_FLAGS} -w")

add_executable(coco ${coco_source_files} ${coco_header_files})
target_link_libraries (coco ${gatb-core-libraries})
