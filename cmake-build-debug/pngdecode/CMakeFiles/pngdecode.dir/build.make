# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.17

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

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2020.2.1\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2020.2.1\bin\cmake\win\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\cmake-build-debug"

# Include any dependencies generated for this target.
include pngdecode/CMakeFiles/pngdecode.dir/depend.make

# Include the progress variables for this target.
include pngdecode/CMakeFiles/pngdecode.dir/progress.make

# Include the compile flags for this target's objects.
include pngdecode/CMakeFiles/pngdecode.dir/flags.make

pngdecode/CMakeFiles/pngdecode.dir/png.cpp.obj: pngdecode/CMakeFiles/pngdecode.dir/flags.make
pngdecode/CMakeFiles/pngdecode.dir/png.cpp.obj: pngdecode/CMakeFiles/pngdecode.dir/includes_CXX.rsp
pngdecode/CMakeFiles/pngdecode.dir/png.cpp.obj: ../pngdecode/png.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object pngdecode/CMakeFiles/pngdecode.dir/png.cpp.obj"
	cd /d C:\Users\FeuFeve\CLIONP~1\Master\S2-PRO~1\RAYTRA~1\CMAKE-~1\PNGDEC~1 && C:\PROGRA~2\mingw32\bin\G__~1.EXE  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\pngdecode.dir\png.cpp.obj -c "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\pngdecode\png.cpp"

pngdecode/CMakeFiles/pngdecode.dir/png.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pngdecode.dir/png.cpp.i"
	cd /d C:\Users\FeuFeve\CLIONP~1\Master\S2-PRO~1\RAYTRA~1\CMAKE-~1\PNGDEC~1 && C:\PROGRA~2\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\pngdecode\png.cpp" > CMakeFiles\pngdecode.dir\png.cpp.i

pngdecode/CMakeFiles/pngdecode.dir/png.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pngdecode.dir/png.cpp.s"
	cd /d C:\Users\FeuFeve\CLIONP~1\Master\S2-PRO~1\RAYTRA~1\CMAKE-~1\PNGDEC~1 && C:\PROGRA~2\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\pngdecode\png.cpp" -o CMakeFiles\pngdecode.dir\png.cpp.s

pngdecode/CMakeFiles/pngdecode.dir/pngdec.cpp.obj: pngdecode/CMakeFiles/pngdecode.dir/flags.make
pngdecode/CMakeFiles/pngdecode.dir/pngdec.cpp.obj: pngdecode/CMakeFiles/pngdecode.dir/includes_CXX.rsp
pngdecode/CMakeFiles/pngdecode.dir/pngdec.cpp.obj: ../pngdecode/pngdec.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object pngdecode/CMakeFiles/pngdecode.dir/pngdec.cpp.obj"
	cd /d C:\Users\FeuFeve\CLIONP~1\Master\S2-PRO~1\RAYTRA~1\CMAKE-~1\PNGDEC~1 && C:\PROGRA~2\mingw32\bin\G__~1.EXE  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\pngdecode.dir\pngdec.cpp.obj -c "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\pngdecode\pngdec.cpp"

pngdecode/CMakeFiles/pngdecode.dir/pngdec.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pngdecode.dir/pngdec.cpp.i"
	cd /d C:\Users\FeuFeve\CLIONP~1\Master\S2-PRO~1\RAYTRA~1\CMAKE-~1\PNGDEC~1 && C:\PROGRA~2\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\pngdecode\pngdec.cpp" > CMakeFiles\pngdecode.dir\pngdec.cpp.i

pngdecode/CMakeFiles/pngdecode.dir/pngdec.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pngdecode.dir/pngdec.cpp.s"
	cd /d C:\Users\FeuFeve\CLIONP~1\Master\S2-PRO~1\RAYTRA~1\CMAKE-~1\PNGDEC~1 && C:\PROGRA~2\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\pngdecode\pngdec.cpp" -o CMakeFiles\pngdecode.dir\pngdec.cpp.s

pngdecode/CMakeFiles/pngdecode.dir/pngenc.cpp.obj: pngdecode/CMakeFiles/pngdecode.dir/flags.make
pngdecode/CMakeFiles/pngdecode.dir/pngenc.cpp.obj: pngdecode/CMakeFiles/pngdecode.dir/includes_CXX.rsp
pngdecode/CMakeFiles/pngdecode.dir/pngenc.cpp.obj: ../pngdecode/pngenc.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object pngdecode/CMakeFiles/pngdecode.dir/pngenc.cpp.obj"
	cd /d C:\Users\FeuFeve\CLIONP~1\Master\S2-PRO~1\RAYTRA~1\CMAKE-~1\PNGDEC~1 && C:\PROGRA~2\mingw32\bin\G__~1.EXE  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\pngdecode.dir\pngenc.cpp.obj -c "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\pngdecode\pngenc.cpp"

pngdecode/CMakeFiles/pngdecode.dir/pngenc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pngdecode.dir/pngenc.cpp.i"
	cd /d C:\Users\FeuFeve\CLIONP~1\Master\S2-PRO~1\RAYTRA~1\CMAKE-~1\PNGDEC~1 && C:\PROGRA~2\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\pngdecode\pngenc.cpp" > CMakeFiles\pngdecode.dir\pngenc.cpp.i

pngdecode/CMakeFiles/pngdecode.dir/pngenc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pngdecode.dir/pngenc.cpp.s"
	cd /d C:\Users\FeuFeve\CLIONP~1\Master\S2-PRO~1\RAYTRA~1\CMAKE-~1\PNGDEC~1 && C:\PROGRA~2\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\pngdecode\pngenc.cpp" -o CMakeFiles\pngdecode.dir\pngenc.cpp.s

pngdecode/CMakeFiles/pngdecode.dir/zdec.cpp.obj: pngdecode/CMakeFiles/pngdecode.dir/flags.make
pngdecode/CMakeFiles/pngdecode.dir/zdec.cpp.obj: pngdecode/CMakeFiles/pngdecode.dir/includes_CXX.rsp
pngdecode/CMakeFiles/pngdecode.dir/zdec.cpp.obj: ../pngdecode/zdec.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object pngdecode/CMakeFiles/pngdecode.dir/zdec.cpp.obj"
	cd /d C:\Users\FeuFeve\CLIONP~1\Master\S2-PRO~1\RAYTRA~1\CMAKE-~1\PNGDEC~1 && C:\PROGRA~2\mingw32\bin\G__~1.EXE  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\pngdecode.dir\zdec.cpp.obj -c "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\pngdecode\zdec.cpp"

pngdecode/CMakeFiles/pngdecode.dir/zdec.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pngdecode.dir/zdec.cpp.i"
	cd /d C:\Users\FeuFeve\CLIONP~1\Master\S2-PRO~1\RAYTRA~1\CMAKE-~1\PNGDEC~1 && C:\PROGRA~2\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\pngdecode\zdec.cpp" > CMakeFiles\pngdecode.dir\zdec.cpp.i

pngdecode/CMakeFiles/pngdecode.dir/zdec.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pngdecode.dir/zdec.cpp.s"
	cd /d C:\Users\FeuFeve\CLIONP~1\Master\S2-PRO~1\RAYTRA~1\CMAKE-~1\PNGDEC~1 && C:\PROGRA~2\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\pngdecode\zdec.cpp" -o CMakeFiles\pngdecode.dir\zdec.cpp.s

pngdecode/CMakeFiles/pngdecode.dir/zenc.cpp.obj: pngdecode/CMakeFiles/pngdecode.dir/flags.make
pngdecode/CMakeFiles/pngdecode.dir/zenc.cpp.obj: pngdecode/CMakeFiles/pngdecode.dir/includes_CXX.rsp
pngdecode/CMakeFiles/pngdecode.dir/zenc.cpp.obj: ../pngdecode/zenc.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object pngdecode/CMakeFiles/pngdecode.dir/zenc.cpp.obj"
	cd /d C:\Users\FeuFeve\CLIONP~1\Master\S2-PRO~1\RAYTRA~1\CMAKE-~1\PNGDEC~1 && C:\PROGRA~2\mingw32\bin\G__~1.EXE  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\pngdecode.dir\zenc.cpp.obj -c "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\pngdecode\zenc.cpp"

pngdecode/CMakeFiles/pngdecode.dir/zenc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pngdecode.dir/zenc.cpp.i"
	cd /d C:\Users\FeuFeve\CLIONP~1\Master\S2-PRO~1\RAYTRA~1\CMAKE-~1\PNGDEC~1 && C:\PROGRA~2\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\pngdecode\zenc.cpp" > CMakeFiles\pngdecode.dir\zenc.cpp.i

pngdecode/CMakeFiles/pngdecode.dir/zenc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pngdecode.dir/zenc.cpp.s"
	cd /d C:\Users\FeuFeve\CLIONP~1\Master\S2-PRO~1\RAYTRA~1\CMAKE-~1\PNGDEC~1 && C:\PROGRA~2\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\pngdecode\zenc.cpp" -o CMakeFiles\pngdecode.dir\zenc.cpp.s

pngdecode/CMakeFiles/pngdecode.dir/zss.cpp.obj: pngdecode/CMakeFiles/pngdecode.dir/flags.make
pngdecode/CMakeFiles/pngdecode.dir/zss.cpp.obj: pngdecode/CMakeFiles/pngdecode.dir/includes_CXX.rsp
pngdecode/CMakeFiles/pngdecode.dir/zss.cpp.obj: ../pngdecode/zss.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object pngdecode/CMakeFiles/pngdecode.dir/zss.cpp.obj"
	cd /d C:\Users\FeuFeve\CLIONP~1\Master\S2-PRO~1\RAYTRA~1\CMAKE-~1\PNGDEC~1 && C:\PROGRA~2\mingw32\bin\G__~1.EXE  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\pngdecode.dir\zss.cpp.obj -c "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\pngdecode\zss.cpp"

pngdecode/CMakeFiles/pngdecode.dir/zss.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pngdecode.dir/zss.cpp.i"
	cd /d C:\Users\FeuFeve\CLIONP~1\Master\S2-PRO~1\RAYTRA~1\CMAKE-~1\PNGDEC~1 && C:\PROGRA~2\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\pngdecode\zss.cpp" > CMakeFiles\pngdecode.dir\zss.cpp.i

pngdecode/CMakeFiles/pngdecode.dir/zss.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pngdecode.dir/zss.cpp.s"
	cd /d C:\Users\FeuFeve\CLIONP~1\Master\S2-PRO~1\RAYTRA~1\CMAKE-~1\PNGDEC~1 && C:\PROGRA~2\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\pngdecode\zss.cpp" -o CMakeFiles\pngdecode.dir\zss.cpp.s

# Object files for target pngdecode
pngdecode_OBJECTS = \
"CMakeFiles/pngdecode.dir/png.cpp.obj" \
"CMakeFiles/pngdecode.dir/pngdec.cpp.obj" \
"CMakeFiles/pngdecode.dir/pngenc.cpp.obj" \
"CMakeFiles/pngdecode.dir/zdec.cpp.obj" \
"CMakeFiles/pngdecode.dir/zenc.cpp.obj" \
"CMakeFiles/pngdecode.dir/zss.cpp.obj"

# External object files for target pngdecode
pngdecode_EXTERNAL_OBJECTS =

pngdecode/libpngdecode.a: pngdecode/CMakeFiles/pngdecode.dir/png.cpp.obj
pngdecode/libpngdecode.a: pngdecode/CMakeFiles/pngdecode.dir/pngdec.cpp.obj
pngdecode/libpngdecode.a: pngdecode/CMakeFiles/pngdecode.dir/pngenc.cpp.obj
pngdecode/libpngdecode.a: pngdecode/CMakeFiles/pngdecode.dir/zdec.cpp.obj
pngdecode/libpngdecode.a: pngdecode/CMakeFiles/pngdecode.dir/zenc.cpp.obj
pngdecode/libpngdecode.a: pngdecode/CMakeFiles/pngdecode.dir/zss.cpp.obj
pngdecode/libpngdecode.a: pngdecode/CMakeFiles/pngdecode.dir/build.make
pngdecode/libpngdecode.a: pngdecode/CMakeFiles/pngdecode.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX static library libpngdecode.a"
	cd /d C:\Users\FeuFeve\CLIONP~1\Master\S2-PRO~1\RAYTRA~1\CMAKE-~1\PNGDEC~1 && $(CMAKE_COMMAND) -P CMakeFiles\pngdecode.dir\cmake_clean_target.cmake
	cd /d C:\Users\FeuFeve\CLIONP~1\Master\S2-PRO~1\RAYTRA~1\CMAKE-~1\PNGDEC~1 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\pngdecode.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
pngdecode/CMakeFiles/pngdecode.dir/build: pngdecode/libpngdecode.a

.PHONY : pngdecode/CMakeFiles/pngdecode.dir/build

pngdecode/CMakeFiles/pngdecode.dir/clean:
	cd /d C:\Users\FeuFeve\CLIONP~1\Master\S2-PRO~1\RAYTRA~1\CMAKE-~1\PNGDEC~1 && $(CMAKE_COMMAND) -P CMakeFiles\pngdecode.dir\cmake_clean.cmake
.PHONY : pngdecode/CMakeFiles/pngdecode.dir/clean

pngdecode/CMakeFiles/pngdecode.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template" "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\pngdecode" "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\cmake-build-debug" "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\cmake-build-debug\pngdecode" "C:\Users\FeuFeve\CLionProjects\Master\S2 - Programmation 3D\raytracer-template\cmake-build-debug\pngdecode\CMakeFiles\pngdecode.dir\DependInfo.cmake" --color=$(COLOR)
.PHONY : pngdecode/CMakeFiles/pngdecode.dir/depend

