# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/cmake-3.15.3-f32dvhh7wjnswxcvedgni64eatk4lueo/bin/cmake

# The command to remove a file.
RM = /opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/cmake-3.15.3-f32dvhh7wjnswxcvedgni64eatk4lueo/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /lore/adesoa/dev/wdmapp_coupling/XGC_adios2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /lore/adesoa/dev/wdmapp_coupling/XGC_adios2/build

# Utility rule file for ExperimentalMemCheck.

# Include the progress variables for this target.
include CMakeFiles/ExperimentalMemCheck.dir/progress.make

CMakeFiles/ExperimentalMemCheck:
	/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/cmake-3.15.3-f32dvhh7wjnswxcvedgni64eatk4lueo/bin/ctest -D ExperimentalMemCheck

ExperimentalMemCheck: CMakeFiles/ExperimentalMemCheck
ExperimentalMemCheck: CMakeFiles/ExperimentalMemCheck.dir/build.make

.PHONY : ExperimentalMemCheck

# Rule to build all files generated by this target.
CMakeFiles/ExperimentalMemCheck.dir/build: ExperimentalMemCheck

.PHONY : CMakeFiles/ExperimentalMemCheck.dir/build

CMakeFiles/ExperimentalMemCheck.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ExperimentalMemCheck.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ExperimentalMemCheck.dir/clean

CMakeFiles/ExperimentalMemCheck.dir/depend:
	cd /lore/adesoa/dev/wdmapp_coupling/XGC_adios2/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /lore/adesoa/dev/wdmapp_coupling/XGC_adios2 /lore/adesoa/dev/wdmapp_coupling/XGC_adios2 /lore/adesoa/dev/wdmapp_coupling/XGC_adios2/build /lore/adesoa/dev/wdmapp_coupling/XGC_adios2/build /lore/adesoa/dev/wdmapp_coupling/XGC_adios2/build/CMakeFiles/ExperimentalMemCheck.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ExperimentalMemCheck.dir/depend

