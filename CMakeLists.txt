cmake_minimum_required(VERSION 3.20)
include(FindPackageHandleStandardArgs)
# ---------------------------------------------------------------------
project(nQMCC-FFT VERSION 1.0
        DESCRIPTION ""
        LANGUAGES Fortran)
# ---------------------------------------------------------------------
enable_language(Fortran)
# ---------------------------------------------------------------------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY lib)
set(CMAKE_INSTALL_PREFIX ".") 
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)
# ---------------------------------------------------------------------
add_library(FFT  fft-nqmcc.f)
# ---------------------------------------------------------------------
if(NOT DEFINED FFT_MASTER_PROJECT)
  set(FFT_MASTER_PROJECT OFF)
endif()

if(FFT_MASTER_PROJECT)
  add_executable(fft test-fft.f90)
  target_link_libraries(fft PUBLIC FFT)
endif()

