# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.21

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/pavan/local/bin/cmake

# The command to remove a file.
RM = /home/pavan/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/pavan/academics/sp2022/cse402/project/cse402_project

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/pavan/academics/sp2022/cse402/project/cse402_project/build

# Include any dependencies generated for this target.
include src/openmp/CMakeFiles/helpersomp.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/openmp/CMakeFiles/helpersomp.dir/compiler_depend.make

# Include the progress variables for this target.
include src/openmp/CMakeFiles/helpersomp.dir/progress.make

# Include the compile flags for this target's objects.
include src/openmp/CMakeFiles/helpersomp.dir/flags.make

src/openmp/CMakeFiles/helpersomp.dir/__/helpers/helpers.cpp.o: src/openmp/CMakeFiles/helpersomp.dir/flags.make
src/openmp/CMakeFiles/helpersomp.dir/__/helpers/helpers.cpp.o: ../src/helpers/helpers.cpp
src/openmp/CMakeFiles/helpersomp.dir/__/helpers/helpers.cpp.o: src/openmp/CMakeFiles/helpersomp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pavan/academics/sp2022/cse402/project/cse402_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/openmp/CMakeFiles/helpersomp.dir/__/helpers/helpers.cpp.o"
	cd /home/pavan/academics/sp2022/cse402/project/cse402_project/build/src/openmp && g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/openmp/CMakeFiles/helpersomp.dir/__/helpers/helpers.cpp.o -MF CMakeFiles/helpersomp.dir/__/helpers/helpers.cpp.o.d -o CMakeFiles/helpersomp.dir/__/helpers/helpers.cpp.o -c /home/pavan/academics/sp2022/cse402/project/cse402_project/src/helpers/helpers.cpp

src/openmp/CMakeFiles/helpersomp.dir/__/helpers/helpers.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/helpersomp.dir/__/helpers/helpers.cpp.i"
	cd /home/pavan/academics/sp2022/cse402/project/cse402_project/build/src/openmp && g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/pavan/academics/sp2022/cse402/project/cse402_project/src/helpers/helpers.cpp > CMakeFiles/helpersomp.dir/__/helpers/helpers.cpp.i

src/openmp/CMakeFiles/helpersomp.dir/__/helpers/helpers.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/helpersomp.dir/__/helpers/helpers.cpp.s"
	cd /home/pavan/academics/sp2022/cse402/project/cse402_project/build/src/openmp && g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/pavan/academics/sp2022/cse402/project/cse402_project/src/helpers/helpers.cpp -o CMakeFiles/helpersomp.dir/__/helpers/helpers.cpp.s

# Object files for target helpersomp
helpersomp_OBJECTS = \
"CMakeFiles/helpersomp.dir/__/helpers/helpers.cpp.o"

# External object files for target helpersomp
helpersomp_EXTERNAL_OBJECTS =

src/openmp/libhelpersomp.a: src/openmp/CMakeFiles/helpersomp.dir/__/helpers/helpers.cpp.o
src/openmp/libhelpersomp.a: src/openmp/CMakeFiles/helpersomp.dir/build.make
src/openmp/libhelpersomp.a: src/openmp/CMakeFiles/helpersomp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/pavan/academics/sp2022/cse402/project/cse402_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libhelpersomp.a"
	cd /home/pavan/academics/sp2022/cse402/project/cse402_project/build/src/openmp && $(CMAKE_COMMAND) -P CMakeFiles/helpersomp.dir/cmake_clean_target.cmake
	cd /home/pavan/academics/sp2022/cse402/project/cse402_project/build/src/openmp && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/helpersomp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/openmp/CMakeFiles/helpersomp.dir/build: src/openmp/libhelpersomp.a
.PHONY : src/openmp/CMakeFiles/helpersomp.dir/build

src/openmp/CMakeFiles/helpersomp.dir/clean:
	cd /home/pavan/academics/sp2022/cse402/project/cse402_project/build/src/openmp && $(CMAKE_COMMAND) -P CMakeFiles/helpersomp.dir/cmake_clean.cmake
.PHONY : src/openmp/CMakeFiles/helpersomp.dir/clean

src/openmp/CMakeFiles/helpersomp.dir/depend:
	cd /home/pavan/academics/sp2022/cse402/project/cse402_project/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/pavan/academics/sp2022/cse402/project/cse402_project /home/pavan/academics/sp2022/cse402/project/cse402_project/src/openmp /home/pavan/academics/sp2022/cse402/project/cse402_project/build /home/pavan/academics/sp2022/cse402/project/cse402_project/build/src/openmp /home/pavan/academics/sp2022/cse402/project/cse402_project/build/src/openmp/CMakeFiles/helpersomp.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/openmp/CMakeFiles/helpersomp.dir/depend
