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
CMAKE_COMMAND = /opt/apps/cmake/2.8.9/bin/cmake

# The command to remove a file.
RM = /opt/apps/cmake/2.8.9/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /opt/apps/cmake/2.8.9/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /work/00004/michoski/arcOn-r2-11.13/arcOn/ArcOn

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /work/00004/michoski/arcOn-r2-11.13/arcOn/ArcOn

# Include any dependencies generated for this target.
include CMakeFiles/./lib/2d/ArcOn.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/./lib/2d/ArcOn.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/./lib/2d/ArcOn.dir/flags.make

CMakeFiles/./lib/2d/ArcOn.dir/source/Main.cc.o: CMakeFiles/./lib/2d/ArcOn.dir/flags.make
CMakeFiles/./lib/2d/ArcOn.dir/source/Main.cc.o: source/Main.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /work/00004/michoski/arcOn-r2-11.13/arcOn/ArcOn/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/./lib/2d/ArcOn.dir/source/Main.cc.o"
	/opt/apps/intel13/mvapich2/1.9/bin/mpicxx   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/./lib/2d/ArcOn.dir/source/Main.cc.o -c /work/00004/michoski/arcOn-r2-11.13/arcOn/ArcOn/source/Main.cc

CMakeFiles/./lib/2d/ArcOn.dir/source/Main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/./lib/2d/ArcOn.dir/source/Main.cc.i"
	/opt/apps/intel13/mvapich2/1.9/bin/mpicxx  $(CXX_DEFINES) $(CXX_FLAGS) -E /work/00004/michoski/arcOn-r2-11.13/arcOn/ArcOn/source/Main.cc > CMakeFiles/./lib/2d/ArcOn.dir/source/Main.cc.i

CMakeFiles/./lib/2d/ArcOn.dir/source/Main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/./lib/2d/ArcOn.dir/source/Main.cc.s"
	/opt/apps/intel13/mvapich2/1.9/bin/mpicxx  $(CXX_DEFINES) $(CXX_FLAGS) -S /work/00004/michoski/arcOn-r2-11.13/arcOn/ArcOn/source/Main.cc -o CMakeFiles/./lib/2d/ArcOn.dir/source/Main.cc.s

CMakeFiles/./lib/2d/ArcOn.dir/source/Main.cc.o.requires:
.PHONY : CMakeFiles/./lib/2d/ArcOn.dir/source/Main.cc.o.requires

CMakeFiles/./lib/2d/ArcOn.dir/source/Main.cc.o.provides: CMakeFiles/./lib/2d/ArcOn.dir/source/Main.cc.o.requires
	$(MAKE) -f CMakeFiles/./lib/2d/ArcOn.dir/build.make CMakeFiles/./lib/2d/ArcOn.dir/source/Main.cc.o.provides.build
.PHONY : CMakeFiles/./lib/2d/ArcOn.dir/source/Main.cc.o.provides

CMakeFiles/./lib/2d/ArcOn.dir/source/Main.cc.o.provides.build: CMakeFiles/./lib/2d/ArcOn.dir/source/Main.cc.o

# Object files for target ./lib/2d/ArcOn
_/lib/2d/ArcOn_OBJECTS = \
"CMakeFiles/./lib/2d/ArcOn.dir/source/Main.cc.o"

# External object files for target ./lib/2d/ArcOn
_/lib/2d/ArcOn_EXTERNAL_OBJECTS =

./lib/2d/ArcOn: CMakeFiles/./lib/2d/ArcOn.dir/source/Main.cc.o
./lib/2d/ArcOn: CMakeFiles/./lib/2d/ArcOn.dir/build.make
./lib/2d/ArcOn: /work/00004/michoski/arcOn-r2-11.13/arcOn/trigger/lib/libdeal_II.so.8.1.pre
./lib/2d/ArcOn: /work/00004/michoski/arcOn-r2-11.13/arcOn/p4est/FAST/lib/libp4est.so
./lib/2d/ArcOn: /work/00004/michoski/arcOn-r2-11.13/arcOn/p4est/FAST/lib/libsc.so
./lib/2d/ArcOn: /work/00004/michoski/arcOn-r2-11.13/arcOn/petsc-3.4.3/intel-13.0.2.146-MVAPICH2-1.9a2-cxx-opt/lib/libpetsc.so
./lib/2d/ArcOn: /work/00004/michoski/arcOn-r2-11.13/arcOn/petsc-3.4.3/intel-13.0.2.146-MVAPICH2-1.9a2-cxx-opt/lib/libspai.a
./lib/2d/ArcOn: /work/00004/michoski/arcOn-r2-11.13/arcOn/petsc-3.4.3/intel-13.0.2.146-MVAPICH2-1.9a2-cxx-opt/lib/libcmumps.a
./lib/2d/ArcOn: /work/00004/michoski/arcOn-r2-11.13/arcOn/petsc-3.4.3/intel-13.0.2.146-MVAPICH2-1.9a2-cxx-opt/lib/libdmumps.a
./lib/2d/ArcOn: /work/00004/michoski/arcOn-r2-11.13/arcOn/petsc-3.4.3/intel-13.0.2.146-MVAPICH2-1.9a2-cxx-opt/lib/libsmumps.a
./lib/2d/ArcOn: /work/00004/michoski/arcOn-r2-11.13/arcOn/petsc-3.4.3/intel-13.0.2.146-MVAPICH2-1.9a2-cxx-opt/lib/libzmumps.a
./lib/2d/ArcOn: /work/00004/michoski/arcOn-r2-11.13/arcOn/petsc-3.4.3/intel-13.0.2.146-MVAPICH2-1.9a2-cxx-opt/lib/libmumps_common.a
./lib/2d/ArcOn: /work/00004/michoski/arcOn-r2-11.13/arcOn/petsc-3.4.3/intel-13.0.2.146-MVAPICH2-1.9a2-cxx-opt/lib/libpord.a
./lib/2d/ArcOn: /work/00004/michoski/arcOn-r2-11.13/arcOn/petsc-3.4.3/intel-13.0.2.146-MVAPICH2-1.9a2-cxx-opt/lib/libscalapack.a
./lib/2d/ArcOn: /work/00004/michoski/arcOn-r2-11.13/arcOn/petsc-3.4.3/intel-13.0.2.146-MVAPICH2-1.9a2-cxx-opt/lib/libHYPRE.a
./lib/2d/ArcOn: /work/00004/michoski/arcOn-r2-11.13/arcOn/petsc-3.4.3/intel-13.0.2.146-MVAPICH2-1.9a2-cxx-opt/lib/libsuperlu_4.3.a
./lib/2d/ArcOn: /work/00004/michoski/arcOn-r2-11.13/arcOn/petsc-3.4.3/intel-13.0.2.146-MVAPICH2-1.9a2-cxx-opt/lib/libsuperlu_dist_3.3.a
./lib/2d/ArcOn: /opt/apps/intel/13/composer_xe_2013.2.146/mkl/lib/intel64/libmkl_sequential.so
./lib/2d/ArcOn: /usr/lib64/libX11.so
./lib/2d/ArcOn: /work/00004/michoski/arcOn-r2-11.13/arcOn/petsc-3.4.3/intel-13.0.2.146-MVAPICH2-1.9a2-cxx-opt/lib/libparmetis.so
./lib/2d/ArcOn: /work/00004/michoski/arcOn-r2-11.13/arcOn/petsc-3.4.3/intel-13.0.2.146-MVAPICH2-1.9a2-cxx-opt/lib/libmetis.so
./lib/2d/ArcOn: /opt/apps/intel13/mvapich2/1.9/lib/libmpichcxx.so
./lib/2d/ArcOn: /opt/apps/intel/13/composer_xe_2013.2.146/compiler/lib/intel64/libirng.so
./lib/2d/ArcOn: /opt/apps/intel/13/composer_xe_2013.2.146/compiler/lib/intel64/libdecimal.a
./lib/2d/ArcOn: /opt/apps/intel/13/composer_xe_2013.2.146/compiler/lib/intel64/libcilkrts.so
./lib/2d/ArcOn: /opt/apps/intel/13/composer_xe_2013.2.146/compiler/lib/intel64/libirc.so
./lib/2d/ArcOn: /opt/apps/intel/13/composer_xe_2013.2.146/mkl/lib/intel64/libmkl_gf_lp64.so
./lib/2d/ArcOn: /opt/apps/intel/13/composer_xe_2013.2.146/mkl/lib/intel64/libmkl_intel_lp64.so
./lib/2d/ArcOn: /opt/apps/intel/13/composer_xe_2013.2.146/mkl/lib/intel64/libmkl_intel_thread.so
./lib/2d/ArcOn: /opt/apps/intel/13/composer_xe_2013.2.146/mkl/lib/intel64/libmkl_core.so
./lib/2d/ArcOn: /opt/apps/intel/13/composer_xe_2013.2.146/compiler/lib/intel64/libiomp5.so
./lib/2d/ArcOn: /opt/apps/intel13/mvapich2/1.9/lib/libmpichf90.so
./lib/2d/ArcOn: /opt/apps/intel/13/composer_xe_2013.2.146/compiler/lib/intel64/libifport.so
./lib/2d/ArcOn: /opt/apps/intel/13/composer_xe_2013.2.146/compiler/lib/intel64/libifcore.so
./lib/2d/ArcOn: /opt/apps/intel/13/composer_xe_2013.2.146/compiler/lib/intel64/libimf.so
./lib/2d/ArcOn: /usr/lib64/libm.so
./lib/2d/ArcOn: /opt/apps/intel/13/composer_xe_2013.2.146/compiler/lib/intel64/libipgo.a
./lib/2d/ArcOn: /opt/apps/intel/13/composer_xe_2013.2.146/compiler/lib/intel64/libintlc.so
./lib/2d/ArcOn: /opt/apps/intel/13/composer_xe_2013.2.146/compiler/lib/intel64/libsvml.so
./lib/2d/ArcOn: /opt/apps/intel/13/composer_xe_2013.2.146/compiler/lib/intel64/libirc_s.a
./lib/2d/ArcOn: /usr/lib64/libdl.so
./lib/2d/ArcOn: /usr/lib64/libc.so
./lib/2d/ArcOn: /opt/apps/intel13/mvapich2/1.9/lib/libmpich.so
./lib/2d/ArcOn: /opt/apps/intel13/mvapich2/1.9/lib/libopa.so
./lib/2d/ArcOn: /opt/apps/intel13/mvapich2/1.9/lib/libmpl.so
./lib/2d/ArcOn: /opt/ofed/lib64/libibmad.so
./lib/2d/ArcOn: /opt/ofed/lib64/librdmacm.so
./lib/2d/ArcOn: /opt/ofed/lib64/libibumad.so
./lib/2d/ArcOn: /opt/ofed/lib64/libibverbs.so
./lib/2d/ArcOn: /usr/lib64/librt.so
./lib/2d/ArcOn: /opt/apps/limic2/0.5.5/lib/liblimic2.so
./lib/2d/ArcOn: /usr/lib64/libpthread.so
./lib/2d/ArcOn: /usr/lib64/libz.so
./lib/2d/ArcOn: CMakeFiles/./lib/2d/ArcOn.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ./lib/2d/ArcOn"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/./lib/2d/ArcOn.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/./lib/2d/ArcOn.dir/build: ./lib/2d/ArcOn
.PHONY : CMakeFiles/./lib/2d/ArcOn.dir/build

CMakeFiles/./lib/2d/ArcOn.dir/requires: CMakeFiles/./lib/2d/ArcOn.dir/source/Main.cc.o.requires
.PHONY : CMakeFiles/./lib/2d/ArcOn.dir/requires

CMakeFiles/./lib/2d/ArcOn.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/./lib/2d/ArcOn.dir/cmake_clean.cmake
.PHONY : CMakeFiles/./lib/2d/ArcOn.dir/clean

CMakeFiles/./lib/2d/ArcOn.dir/depend:
	cd /work/00004/michoski/arcOn-r2-11.13/arcOn/ArcOn && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /work/00004/michoski/arcOn-r2-11.13/arcOn/ArcOn /work/00004/michoski/arcOn-r2-11.13/arcOn/ArcOn /work/00004/michoski/arcOn-r2-11.13/arcOn/ArcOn /work/00004/michoski/arcOn-r2-11.13/arcOn/ArcOn /work/00004/michoski/arcOn-r2-11.13/arcOn/ArcOn/CMakeFiles/lib/2d/ArcOn.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/./lib/2d/ArcOn.dir/depend

