# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /tmp/tmp.l9RbY9S9dM

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /tmp/tmp.l9RbY9S9dM/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/challenge.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/challenge.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/challenge.dir/flags.make

CMakeFiles/challenge.dir/main.cpp.o: CMakeFiles/challenge.dir/flags.make
CMakeFiles/challenge.dir/main.cpp.o: ../main.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /tmp/tmp.l9RbY9S9dM/cmake-build-debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/challenge.dir/main.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/challenge.dir/main.cpp.o -c /tmp/tmp.l9RbY9S9dM/main.cpp

CMakeFiles/challenge.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/challenge.dir/main.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /tmp/tmp.l9RbY9S9dM/main.cpp > CMakeFiles/challenge.dir/main.cpp.i

CMakeFiles/challenge.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/challenge.dir/main.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /tmp/tmp.l9RbY9S9dM/main.cpp -o CMakeFiles/challenge.dir/main.cpp.s

CMakeFiles/challenge.dir/main.cpp.o.requires:
.PHONY : CMakeFiles/challenge.dir/main.cpp.o.requires

CMakeFiles/challenge.dir/main.cpp.o.provides: CMakeFiles/challenge.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/challenge.dir/build.make CMakeFiles/challenge.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/challenge.dir/main.cpp.o.provides

CMakeFiles/challenge.dir/main.cpp.o.provides.build: CMakeFiles/challenge.dir/main.cpp.o

# Object files for target challenge
challenge_OBJECTS = \
"CMakeFiles/challenge.dir/main.cpp.o"

# External object files for target challenge
challenge_EXTERNAL_OBJECTS =

challenge: CMakeFiles/challenge.dir/main.cpp.o
challenge: CMakeFiles/challenge.dir/build.make
challenge: CMakeFiles/challenge.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable challenge"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/challenge.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/challenge.dir/build: challenge
.PHONY : CMakeFiles/challenge.dir/build

CMakeFiles/challenge.dir/requires: CMakeFiles/challenge.dir/main.cpp.o.requires
.PHONY : CMakeFiles/challenge.dir/requires

CMakeFiles/challenge.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/challenge.dir/cmake_clean.cmake
.PHONY : CMakeFiles/challenge.dir/clean

CMakeFiles/challenge.dir/depend:
	cd /tmp/tmp.l9RbY9S9dM/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /tmp/tmp.l9RbY9S9dM /tmp/tmp.l9RbY9S9dM /tmp/tmp.l9RbY9S9dM/cmake-build-debug /tmp/tmp.l9RbY9S9dM/cmake-build-debug /tmp/tmp.l9RbY9S9dM/cmake-build-debug/CMakeFiles/challenge.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/challenge.dir/depend

