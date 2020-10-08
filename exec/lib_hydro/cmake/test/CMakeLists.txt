cmake_minimum_required(VERSION 3.17)

project(hello_world)

enable_language(Fortran)

find_package(AMReX)


# Enable OpenMP
find_package(OpenMP)
# If OpenMP is available => Add flags compiler + linker
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_executable(hello.ex hello.cpp)
target_link_libraries(hello.ex
    AMReX::amrex AMReX::Flags_CXX AMReX::Flags_Fortran AMReX::Flags_FPE)
