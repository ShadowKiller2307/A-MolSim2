cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# define project name, version
project(PSEMolDyn_GroupC VERSION 0.0.1)

option(ENABLE_OPENMP "Enable OpenMP" OFF)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE
        PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# set Release as the default build type if it is not yet set.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
            "Release"
            CACHE
            STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
            FORCE)

endif(NOT CMAKE_BUILD_TYPE)

# Add O3 for the Release build type
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O3)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-pg)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -pg")
endif()

# modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/)

if(ENABLE_OPENMP)
    message(STATUS "OpenMP is enabled")
    find_package(OpenMP REQUIRED)
    if(OPENMP_FOUND)
        message(STATUS "OpenMP enabled and found")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

    endif()
else ()
    message(STATUS "OpenMP is disabled")
endif ()

# project dependencies
include(boost_program_options)
include(libxsd)
include(spdlog)
include(xercesc)

# tests and documentation
include(gtest) #provides a BUILD_TESTS option
include(doxygen) #provides a BUILD_DOC_DOXYGEN option

# convenience target to format all source files
include(clang-format) #provides a clangformat target

# add subdirectories
add_subdirectory(src project)

add_subdirectory(tests)

add_subdirectory(benchmarks)
#[[
option(SUB "Subdomain Strategy" OFF) # option for the Subdomain parallelization strategy
option(PAR "Particles Strategy" OFF) # option for the Particles parallelization strategy

if (SUB)
    option(Strategy_1 "Subdomain Strategy" ON)
    option(Strategy_2 "Particles Strategy" OFF)
endif ()
if (PAR)
    option(Strategy_2 "Particles Strategy" ON)
    option(Strategy_1 "Subdomain Strategy" OFF)
endif ()
if (Strategy_1)
    message("Strategy 1 - Subdomains parallelization, was chosen!")
endif ()
if (Strategy_2)
    message("Strategy 2 - Particles parallelization, was chosen")
endif ()
#target_compile_options(MolSim PRIVATE -DPARTICLES=ON)
target_compile_options(MolSim PRIVATE $<$<BOOL:${Strategie_1}>:-DSUBDOMAIN=ON>) #-DCELLMAIN=OFF>)
target_compile_options(MolSim PRIVATE $<$<BOOL:${Strategie_2}>:-DPARTICLES=ON>) # -DSUBMAIN=OFF>)
unset(SUB CACHE)
unset(PAR CACHE)
unset(Strategy_1 CACHE)
unset(Strategy_2 CACHE)]]
