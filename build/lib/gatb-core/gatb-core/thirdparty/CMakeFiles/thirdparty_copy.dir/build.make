# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/annika/COCO

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/annika/COCO/build

# Utility rule file for thirdparty_copy.

# Include the progress variables for this target.
include lib/gatb-core/gatb-core/thirdparty/CMakeFiles/thirdparty_copy.dir/progress.make

thirdparty_copy: lib/gatb-core/gatb-core/thirdparty/CMakeFiles/thirdparty_copy.dir/build.make
	cd /home/annika/COCO/build/lib/gatb-core/gatb-core/thirdparty && /usr/bin/cmake -E make_directory /home/annika/COCO/build/lib/gatb-core/gatb-core/include//json
	cd /home/annika/COCO/build/lib/gatb-core/gatb-core/thirdparty && /usr/bin/cmake -E make_directory /home/annika/COCO/build/lib/gatb-core/gatb-core/include//BooPHF
	cd /home/annika/COCO/build/lib/gatb-core/gatb-core/thirdparty && /usr/bin/cmake -E copy_if_different /home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/BooPHF/BooPHF.h /home/annika/COCO/build/lib/gatb-core/gatb-core/include//BooPHF/
	cd /home/annika/COCO/build/lib/gatb-core/gatb-core/thirdparty && /usr/bin/cmake -E copy_if_different /home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/json/json.hpp /home/annika/COCO/build/lib/gatb-core/gatb-core/include//json/
.PHONY : thirdparty_copy

# Rule to build all files generated by this target.
lib/gatb-core/gatb-core/thirdparty/CMakeFiles/thirdparty_copy.dir/build: thirdparty_copy

.PHONY : lib/gatb-core/gatb-core/thirdparty/CMakeFiles/thirdparty_copy.dir/build

lib/gatb-core/gatb-core/thirdparty/CMakeFiles/thirdparty_copy.dir/clean:
	cd /home/annika/COCO/build/lib/gatb-core/gatb-core/thirdparty && $(CMAKE_COMMAND) -P CMakeFiles/thirdparty_copy.dir/cmake_clean.cmake
.PHONY : lib/gatb-core/gatb-core/thirdparty/CMakeFiles/thirdparty_copy.dir/clean

lib/gatb-core/gatb-core/thirdparty/CMakeFiles/thirdparty_copy.dir/depend:
	cd /home/annika/COCO/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/annika/COCO /home/annika/COCO/lib/gatb-core/gatb-core/thirdparty /home/annika/COCO/build /home/annika/COCO/build/lib/gatb-core/gatb-core/thirdparty /home/annika/COCO/build/lib/gatb-core/gatb-core/thirdparty/CMakeFiles/thirdparty_copy.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/gatb-core/gatb-core/thirdparty/CMakeFiles/thirdparty_copy.dir/depend
