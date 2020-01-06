set(CMAKE_C_COMPILER "/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/mpich-3.3-diz4f6ieln25ouifyc7ndtqlfksom6nb/bin/mpicc")
set(CMAKE_C_COMPILER_ARG1 "")
set(CMAKE_C_COMPILER_ID "GNU")
set(CMAKE_C_COMPILER_VERSION "7.3.0")
set(CMAKE_C_COMPILER_VERSION_INTERNAL "")
set(CMAKE_C_COMPILER_WRAPPER "")
set(CMAKE_C_STANDARD_COMPUTED_DEFAULT "11")
set(CMAKE_C_COMPILE_FEATURES "c_std_90;c_function_prototypes;c_std_99;c_restrict;c_variadic_macros;c_std_11;c_static_assert")
set(CMAKE_C90_COMPILE_FEATURES "c_std_90;c_function_prototypes")
set(CMAKE_C99_COMPILE_FEATURES "c_std_99;c_restrict;c_variadic_macros")
set(CMAKE_C11_COMPILE_FEATURES "c_std_11;c_static_assert")

set(CMAKE_C_PLATFORM_ID "Linux")
set(CMAKE_C_SIMULATE_ID "")
set(CMAKE_C_SIMULATE_VERSION "")



set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_C_COMPILER_AR "/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-rhel7_4.8.5/gcc-7.3.0-bt47fwrzijla4xzdx4o4au45yljqptsk/bin/gcc-ar")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_C_COMPILER_RANLIB "/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-rhel7_4.8.5/gcc-7.3.0-bt47fwrzijla4xzdx4o4au45yljqptsk/bin/gcc-ranlib")
set(CMAKE_LINKER "/usr/bin/ld")
set(CMAKE_COMPILER_IS_GNUCC 1)
set(CMAKE_C_COMPILER_LOADED 1)
set(CMAKE_C_COMPILER_WORKS TRUE)
set(CMAKE_C_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_C_COMPILER_ENV_VAR "CC")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_C_COMPILER_ID_RUN 1)
set(CMAKE_C_SOURCE_FILE_EXTENSIONS c;m)
set(CMAKE_C_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_C_LINKER_PREFERENCE 10)

# Save compiler ABI information.
set(CMAKE_C_SIZEOF_DATA_PTR "8")
set(CMAKE_C_COMPILER_ABI "ELF")
set(CMAKE_C_LIBRARY_ARCHITECTURE "")

if(CMAKE_C_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_C_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_C_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_C_COMPILER_ABI}")
endif()

if(CMAKE_C_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()

set(CMAKE_C_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_C_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_C_CL_SHOWINCLUDES_PREFIX}")
endif()





set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "mpi;gcc;gcc_s;c;gcc;gcc_s")
set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/mpich-3.3-diz4f6ieln25ouifyc7ndtqlfksom6nb/lib;/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-rhel7_4.8.5/gcc-7.3.0-bt47fwrzijla4xzdx4o4au45yljqptsk/lib64;/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-rhel7_4.8.5/gcc-7.3.0-bt47fwrzijla4xzdx4o4au45yljqptsk/lib/gcc/x86_64-pc-linux-gnu/7.3.0;/lib64;/usr/lib64;/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/adios-1.13.1-ootbumw6qyec5wyqcd7udfgcz5gewx5p/lib;/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/zoltan-3.83-4eacxwksu7jeuy52v6pqfkv53xqnfivu/lib;/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/petsc-3.11.0-264efxjpnpgl5ls3erg652m6jy4qn4a7/lib;/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/hdf5-1.10.5-teac4rs63rxd4bqexluxhmrlurx2ly4l/lib;/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/hypre-2.15.1-uhae7i4c2x5nwcyciv7l5ge3gcozm6ik/lib;/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/superlu-dist-6.1.1-fypwv2tkqn3wdsgn7nblrrbwul7ttnha/lib;/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/parmetis-4.0.3-iqptmd22kibq7j4txrvy4u6knhhn5ge2/lib;/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/c-blosc-1.16.3-dpeknkadhnptttl3eeh7ye6wnwiyiem4/lib;/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/zfp-0.5.0-jfot4w4neh2ekyp4jk4w6cugue55tg5v/lib;/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/sz-1.4.12.3-d34fy2spbopeobfttkipdgyxzczs2jvr/lib;/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/zlib-1.2.11-ypve67owtmigrrbz3b76gg66xvihmhmx/lib;/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/numactl-2.0.12-eokgwk6xob3uh6lfgf7qxf33i4d5r4b4/lib;/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/openblas-0.3.5-7miavkpfij3m35htbo6sabiobj5euoeb/lib;/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/libiconv-1.15-hwi4p6xjukf73hopk6sg34z6v57fbghx/lib;/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/libpciaccess-0.13.5-gcqhf4ksz5ndlmj5aoydcpexg3sizi4s/lib;/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/libxml2-2.9.8-hcsmz5fu7cm7tfbx73f2o2xgo52wgtwt/lib;/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/metis-5.1.0-rn7k363kqpvbznmxb3jkkejwjcgfqpyu/lib;/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-rhel7_4.8.5/gcc-7.3.0-bt47fwrzijla4xzdx4o4au45yljqptsk/lib;/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-rhel7_4.8.5/emacs-26.1-eic6nou3mc6edpozo4i7raeqg734uuv3/lib")
set(CMAKE_C_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
