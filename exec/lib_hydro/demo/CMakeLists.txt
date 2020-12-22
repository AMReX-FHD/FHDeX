cmake_minimum_required(VERSION 3.17)

# Define Project
project(main)

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

find_library(
    LIBHYDRO
    NAMES hydro
    HINTS ../cmake/dist/lib
    )

message("LIBHYDRO = ${LIBHYDRO}")

# get_target_property(amrex_compile_defs AMReX::amrex  INTERFACE_COMPILE_DEFINITIONS)
# string(REPLACE       "$"     ""    amrex_compile_defs "${amrex_compile_defs}")
# string(REGEX REPLACE "<.*>;" ""    amrex_compile_defs "${amrex_compile_defs}")
# string(REPLACE       ";"     " -D" amrex_compile_defs "-D${amrex_compile_defs}")

# list(APPEND CMAKE_CXX_FLAGS "${amrex_compile_defs}")
# list(APPEND CMAKE_CXX_FLAGS "-O3 -DNDEBUG -Werror=return-type")
# 
# string(REPLACE ";" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
# message(" *** CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
# 
# list(APPEND CMAKE_Fortran_FLAGS "${amrex_compile_defs}")
# list(APPEND CMAKE_Fortran_FLAGS "-O3 -DNDEBUG -O3 -fimplicit-none -ffixed-line-length-none -ffree-line-length-none")
# 
# string(REPLACE ";" " " CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
# message(" *** CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}")


include_directories(
    ../cmake/dist/include
    ../cmake/dist/modules
    )


add_executable(main.ex main.cpp main_driver.cpp WritePlotFile.cpp)

target_link_libraries(
    main.ex
    AMReX::amrex AMReX::Flags_CXX AMReX::Flags_Fortran AMReX::Flags_FPE
    )

target_link_libraries(
    main.ex
    ${LIBHYDRO}
    )
