cmake_minimum_required(VERSION 3.17)

enable_language(Fortran)

project(hello_world)
find_package(AMReX)
add_executable(hello.ex hello.cpp)
target_link_libraries(hello.ex
    AMReX::amrex AMReX::Flags_CXX AMReX::Flags_Fortran AMReX::Flags_FPE)
