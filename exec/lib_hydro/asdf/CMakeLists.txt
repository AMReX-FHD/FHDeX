cmake_minimum_required(VERSION 3.17)

enable_language(Fortran)

project(main)
find_package(AMReX)

find_library(
    LIBHYDRO
    NAMES hydro
    PATHS ../cmake/dist/lib
    )


message("LIBHYDRO = ${LIBHYDRO}")

#
# Based on the AMReX Config -- this is HACK, and should be automated
# TODO: this currently assumes macOS (Darwin) -- add linux
#

set(CMAKE_CXX_FLAGS
    "-DAMREX_USE_MPI"
    "-DBL_USE_MPI"
    "-DAMREX_USE_OMP"
    "-DBL_USE_OMP"
    "-DAMREX_SPACEDIM=3"
    "-DBL_SPACEDIM=3"
    "-DBL_FORT_USE_UNDERSCORE"
    "-DAMREX_USE_EB"
    "-DAMREX_USE_F_INTERFACES"
    "-DBL_USE_F_INTERFACES"
    "-DAMREX_PARTICLES"
    "-O3 -DNDEBUG  -fopenmp -Werror=return-type"
    )

if(APPLE)
    set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS}"
        "-DAMREX_Darwin"
        "-DBL_Darwin"
        )
elseif(UNIX)
    set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS}"
        "-DAMREX_Linux"
        "-DBL_Linux"
        )
endif(APPLE)

string(REPLACE ";" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
message(" *** CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")

set(CMAKE_Fortran_FLAGS
    "-DAMREX_USE_MPI"
    "-DBL_USE_MPI"
    "-DAMREX_USE_OMP"
    "-DBL_USE_OMP"
    "-DAMREX_SPACEDIM=3"
    "-DBL_SPACEDIM=3"
    "-DBL_LANG_FORT"
    "-DAMREX_LANG_FORT"
    "-DBL_FORT_USE_UNDERSCORE"
    "-DAMREX_USE_EB"
    "-DAMREX_USE_F_INTERFACES"
    "-DBL_USE_F_INTERFACES"
    "-DAMREX_PARTICLES"
    "-O3 -DNDEBUG -O3  -fopenmp -fimplicit-none -ffixed-line-length-none -ffree-line-length-none"
    )

if(APPLE)
    set(CMAKE_Fortran_FLAGS
        "${CMAKE_Fortran_FLAGS}"
        "-DBL_Darwin"
        )
elseif(UNIX)
    set(CMAKE_Fortran_FLAGS
        "${CMAKE_Fortran_FLAGS}"
        "-DBL_Linux"
        )
endif(APPLE)

string(REPLACE ";" " " CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
message(" *** CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}")


include_directories(
    ../cmake/dist/include
    ../cmake/dist/modules
    )


add_executable(main.ex main.cpp main_driver.cpp)

target_link_libraries(
    main.ex
    AMReX::amrex AMReX::Flags_CXX AMReX::Flags_Fortran AMReX::Flags_FPE
    )

target_link_libraries(
    main.ex
    ${LIBHYDRO}
    )
