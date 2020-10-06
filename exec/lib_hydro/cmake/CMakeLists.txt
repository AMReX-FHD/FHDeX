cmake_minimum_required(VERSION 3.17)

# Define Project
project(lib_hydro)

# Enable Fortan
enable_language(Fortran)

# Include AMReX
find_package(AMReX)

# Enable OpenMP
find_package(OpenMP)
# If OpenMP is available => Add flags compiler + linker
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

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
    "-DAMREX_Darwin"
    "-DBL_Darwin"
    "-DBL_FORT_USE_UNDERSCORE"
    "-DAMREX_USE_EB"
    "-DAMREX_USE_F_INTERFACES"
    "-DBL_USE_F_INTERFACES"
    "-DAMREX_PARTICLES"
    "-O3 -DNDEBUG -Werror=return-type"
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
    "-DAMREX_Darwin"
    "-DBL_LANG_FORT"
    "-DAMREX_LANG_FORT"
    "-DBL_FORT_USE_UNDERSCORE"
    "-DAMREX_USE_EB"
    "-DAMREX_USE_F_INTERFACES"
    "-DBL_USE_F_INTERFACES"
    "-DAMREX_PARTICLES"
    "-O3 -DNDEBUG -O3 -fimplicit-none -ffixed-line-length-none -ffree-line-length-none"
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


# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(
    ../../../src_common
    ../../../src_common/src_F90
    ../../../src_rng
    ../../../src_gmres
    ../../../src_gmres/src_F90
    ../../../src_hydro
    ../../../src_hydro/src_F90
    ../../../src_immersed-boundary
    ../../../src_immersed-boundary/src_F90
    )

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/modules)


add_library(
    hydro SHARED

    ../../../src_common/BCPhysToMath.cpp
    ../../../src_common/ComputeAverages.cpp
    ../../../src_common/ComputeBasicStats.cpp
    ../../../src_common/ComputeDivAndGrad.cpp
    ../../../src_common/ConvertStag.cpp
    ../../../src_common/Debug.cpp
    ../../../src_common/FindCoords.cpp
    ../../../src_common/FindNodalValues.cpp
    ../../../src_common/MultiFabPhysBC.cpp
    ../../../src_common/SqrtMF.cpp
    ../../../src_common/common_functions.cpp
    ../../../src_common/main.cpp
    ../../../src_common/main_driver.cpp

    ../../../src_common/src_F90/common_namelist.F90
    ../../../src_common/src_F90/compute_averages.F90
    ../../../src_common/src_F90/compute_basic_stats.F90
    ../../../src_common/src_F90/debug_mod.F90
    ../../../src_common/src_F90/find_coords.F90
    ../../../src_common/src_F90/find_nodal_values.F90
    ../../../src_common/src_F90/init_funcs.F90
    ../../../src_common/src_F90/species_type.F90
    ../../../src_common/src_F90/time_step.F90

    ../../../src_rng/MultiFabFillRandom.cpp
    ../../../src_rng/bl_random_c.cpp

    ../../../src_rng/src_F90/bl_random.F90
    ../../../src_rng/src_F90/multifab_fill_random.F90
    ../../../src_rng/src_F90/rng_functions.F90

    ../../../src_gmres/ApplyMatrix.cpp
    ../../../src_gmres/GMRES.cpp
    ../../../src_gmres/MacProj.cpp
    ../../../src_gmres/NormInnerProduct.cpp
    ../../../src_gmres/Precon.cpp
    ../../../src_gmres/StagApplyOp.cpp
    ../../../src_gmres/StagMGSolver.cpp
    ../../../src_gmres/Utility.cpp
    ../../../src_gmres/gmres_functions.cpp

    ../../../src_gmres/src_F90/gmres_namelist.F90

    ../../../src_hydro/AddMFluctuations.cpp
    ../../../src_hydro/ConvertMToUmac.cpp
    ../../../src_hydro/HydroBCs.cpp
    ../../../src_hydro/MacProj_hydro.cpp
    ../../../src_hydro/MkAdvMFluxdiv.cpp
    ../../../src_hydro/MkAdvSFluxdiv.cpp
    ../../../src_hydro/StochMomFlux.cpp
    ../../../src_hydro/advance.cpp

    ../../../src_hydro/src_F90/hydro_bounds_mod.F90

    ../../../src_immersed-boundary/IBCore.cpp
    ../../../src_immersed-boundary/IBMarkerContainer.cpp
    ../../../src_immersed-boundary/IBMultiBlobContainer.cpp
    ../../../src_immersed-boundary/IBParticleContainer.cpp
    ../../../src_immersed-boundary/ib_functions.cpp

    ../../../src_immersed-boundary/src_F90/constant_mod.F90
    ../../../src_immersed-boundary/src_F90/ib_fort_utils.F90
    ../../../src_immersed-boundary/src_F90/immbdy_namelist.F90
    ../../../src_immersed-boundary/src_F90/param_mod.F90
    ../../../src_immersed-boundary/src_F90/set_particle_properties.F90
    )

target_link_libraries(
    hydro
    AMReX::amrex AMReX::Flags_CXX AMReX::Flags_Fortran AMReX::Flags_FPE
    )


set(HEADERS
    ../../../src_common/common_functions.H
    ../../../src_common/common_functions_F.H
    ../../../src_common/common_namespace.H
    ../../../src_common/common_namespace_declarations.H
    ../../../src_common/species.H

    ../../../src_rng/bl_random_c.H
    ../../../src_rng/rng_functions.H
    ../../../src_rng/rng_functions_F.H

    ../../../src_gmres/GMRES.H
    ../../../src_gmres/MacProj.H
    ../../../src_gmres/Precon.H
    ../../../src_gmres/StagMGSolver.H
    ../../../src_gmres/gmres_functions.H
    ../../../src_gmres/gmres_functions_F.H
    ../../../src_gmres/gmres_namespace.H
    ../../../src_gmres/gmres_namespace_declarations.H

    ../../../src_hydro/StochMomFlux.H
    ../../../src_hydro/hydro_functions.H
    ../../../src_hydro/hydro_functions_F.H

    ../../../src_immersed-boundary/IBCore.H
    ../../../src_immersed-boundary/IBMarkerContainer.H
    ../../../src_immersed-boundary/IBMarkerContainerBase.H
    ../../../src_immersed-boundary/IBMarkerContainerBaseI.H
    ../../../src_immersed-boundary/IBMultiBlobContainer.H
    ../../../src_immersed-boundary/IBParticleContainer.H
    ../../../src_immersed-boundary/IBParticleInfo.H
    ../../../src_immersed-boundary/MFUtil.H
    ../../../src_immersed-boundary/ib_functions.H
    ../../../src_immersed-boundary/ib_functions_F.H
    ../../../src_immersed-boundary/immbdy_namespace.H
    ../../../src_immersed-boundary/immbdy_namespace_declarations.H
    ../../../src_immersed-boundary/kernel_functions_K.H
    )

message("HEADERS = ${HEADERS}")

set_target_properties(
    hydro PROPERTIES PUBLIC_HEADER
    "${HEADERS}"
    )


install(
    TARGETS hydro
    LIBRARY DESTINATION  lib
    PUBLIC_HEADER DESTINATION  include
    )

# install(
#     DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION include
#     )
