SET(CMAKE_C_COMPILER "/opt/apps/intel13/mvapich2/1.9/bin/mpicc")
SET(CMAKE_C_COMPILER_ARG1 "")
SET(CMAKE_C_COMPILER_ID "Intel")
SET(CMAKE_C_COMPILER_VERSION "13.1.0.20130121")
SET(CMAKE_C_PLATFORM_ID "Linux")

SET(CMAKE_AR "/usr/bin/ar")
SET(CMAKE_RANLIB "/usr/bin/ranlib")
SET(CMAKE_LINKER "/opt/apps/xalt/0.4.6/bin/ld")
SET(CMAKE_COMPILER_IS_GNUCC )
SET(CMAKE_C_COMPILER_LOADED 1)
SET(CMAKE_COMPILER_IS_MINGW )
SET(CMAKE_COMPILER_IS_CYGWIN )
IF(CMAKE_COMPILER_IS_CYGWIN)
  SET(CYGWIN 1)
  SET(UNIX 1)
ENDIF(CMAKE_COMPILER_IS_CYGWIN)

SET(CMAKE_C_COMPILER_ENV_VAR "CC")

IF(CMAKE_COMPILER_IS_MINGW)
  SET(MINGW 1)
ENDIF(CMAKE_COMPILER_IS_MINGW)
SET(CMAKE_C_COMPILER_ID_RUN 1)
SET(CMAKE_C_SOURCE_FILE_EXTENSIONS c)
SET(CMAKE_C_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
SET(CMAKE_C_LINKER_PREFERENCE 10)

# Save compiler ABI information.
SET(CMAKE_C_SIZEOF_DATA_PTR "8")
SET(CMAKE_C_COMPILER_ABI "ELF")
SET(CMAKE_C_LIBRARY_ARCHITECTURE "")

IF(CMAKE_C_SIZEOF_DATA_PTR)
  SET(CMAKE_SIZEOF_VOID_P "${CMAKE_C_SIZEOF_DATA_PTR}")
ENDIF(CMAKE_C_SIZEOF_DATA_PTR)

IF(CMAKE_C_COMPILER_ABI)
  SET(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_C_COMPILER_ABI}")
ENDIF(CMAKE_C_COMPILER_ABI)

IF(CMAKE_C_LIBRARY_ARCHITECTURE)
  SET(CMAKE_LIBRARY_ARCHITECTURE "")
ENDIF()

SET(CMAKE_C_HAS_ISYSROOT "")


SET(CMAKE_C_IMPLICIT_LINK_LIBRARIES "mpich;opa;mpl;ibmad;rdmacm;ibumad;ibverbs;rt;limic2;pthread;imf;svml;irng;m;ipgo;decimal;cilkrts;stdc++;irc;svml;c;irc_s;dl;c")
SET(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "/opt/ofed/lib64;/lib;/opt/apps/limic2/0.5.5/lib;/opt/apps/intel13/mvapich2/1.9/lib;/opt/apps/intel/13/composer_xe_2013.2.146/compiler/lib/intel64;/usr/lib/gcc/x86_64-redhat-linux/4.4.7;/usr/lib64;/lib64;/usr/lib")



