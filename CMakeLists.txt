# https://cmake.org/cmake/help/v3.6/module/FindPkgConfig.html
# 3.11: add_library() doesn't require dummy source
# 3.12: target_link_libraries() supports OBJECT targets
# 3.20: cmake_path()
cmake_minimum_required(VERSION 3.20)
project(libgtracers
  VERSION 2024.08.1
  DESCRIPTION "libgtracers contains a collection of tracers and associated code for use with the MOM and GOLD ocean models."
  LANGUAGES Fortran
)

include(GNUInstallDirs)

if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel)$")
  message(STATUS "Setting build type to 'Relwithdebinfo' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "Relwithdebinfo"
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|GNU)$")
  message(
    WARNING "Compiler not officially supported: ${CMAKE_Fortran_COMPILER_ID}")
endif()

cmake_path(APPEND CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR} "include")

################################################################################
# Fortran
################################################################################

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # Copied from MOM5/bin/mkmf.template.nci.gfortran
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcray-pointer -fdefault-real-8 -ffree-line-length-none -fno-range-check -Waliasing -Wampersand -Warray-bounds -Wcharacter-truncation -Wconversion -Wline-truncation -Wintrinsics-std -Wsurprising -Wno-tabs -Wunderflow -Wunused-parameter -Wintrinsic-shadow -Wno-align-commons")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -W -fbounds-check")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  # Copied from MOM5/bin/mkmf.template.nci
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-alias -safe-cray-ptr -fpe0 -ftz -assume byterecl -i4 -r8 -traceback -nowarn -check noarg_temp_created -assume nobuffered_io -convert big_endian -grecord-gcc-switches -align all")
  set(CMAKE_Fortran_FLAGS_RELEASE "-g3 -O2 -xCORE-AVX2 -debug all -check none")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g3 -O2 -xCORE-AVX2 -debug all -check none")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g3 -O0 -debug all -check -check noarg_temp_created -check nopointer -warn -warn noerrors -ftrapuv")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  message(FATAL_ERROR "Fortran compiler ${CMAKE_Fortran_COMPILER_ID} is unsupported")

else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

################################################################################

find_package(MPI REQUIRED COMPONENTS Fortran)
find_package(PkgConfig REQUIRED)
pkg_check_modules(NETCDF REQUIRED IMPORTED_TARGET "netcdf-fortran")
pkg_check_modules(FMS REQUIRED IMPORTED_TARGET "FMS")
pkg_check_modules(MOCSY REQUIRED IMPORTED_TARGET "mocsy")

add_library(gtracers STATIC)

target_sources(gtracers PRIVATE
  generic_tracers/FMS_coupler_util.F90
  generic_tracers/FMS_ocmip2_co2calc.F90
  generic_tracers/generic_abiotic.F90
  generic_tracers/generic_age.F90
  generic_tracers/generic_argon.F90
  generic_tracers/generic_BLING.F90
  generic_tracers/generic_blres.F90
  generic_tracers/generic_CFC.F90
  generic_tracers/generic_COBALT.F90
  generic_tracers/generic_ERGOM.F90
  generic_tracers/generic_miniBLING.F90
  generic_tracers/generic_SF6.F90
  generic_tracers/generic_TOPAZ.F90
  generic_tracers/generic_tracer.F90
  generic_tracers/generic_tracer_utils.F90
  generic_tracers/generic_WOMBATlite.F90
)

target_link_libraries(gtracers PUBLIC
			PkgConfig::FMS
			PkgConfig::MOCSY
			PkgConfig::NETCDF)

install(TARGETS gtracers)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ TYPE INCLUDE)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/gfdl-generic-tracers.pc.in
               ${CMAKE_CURRENT_BINARY_DIR}/gfdl-generic-tracers.pc @ONLY)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gfdl-generic-tracers.pc
              DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
              COMPONENT utilities)
