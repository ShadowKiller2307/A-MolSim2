include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB_RECURSE PROJECT_SRC 
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

list(REMOVE_ITEM PROJECT_SRC "${CMAKE_CURRENT_SOURCE_DIR}/MolSim.cpp")

add_library(project_lib ${PROJECT_SRC}
        simulation/SimulationParams.cpp)

target_include_directories(project_lib PUBLIC
    ${libxsd_SOURCE_DIR}
    ${spdlog_SOURCE_DIR}/include)

target_link_libraries(project_lib PUBLIC XercesC::XercesC spdlog::spdlog Boost::program_options)

add_executable(MolSim MolSim.cpp
        particles/containers/linkedcells/Subdomains/Subdomain.cpp
        particles/containers/linkedcells/Subdomains/Subdomain.h)
target_link_libraries(MolSim PUBLIC project_lib)

# activate all compiler warnings.
target_compile_options(
    MolSim PRIVATE 
        $<$<CXX_COMPILER_ID:GNU>:-Wall>
        $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>)

#[[target_compile_definitions(MolSim PRIVATE STRATEGY_1) #custom compiler option for strategy 1
target_compile_definitions(MolSim PRIVATE STRATEGY_2) #custom compiler option for strategy 2]]
#compiler options for different parallelization strategies
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

target_compile_options(MolSim PRIVATE $<$<BOOL:${Strategie_1}>:-DSUBDOMAIN=ON>) #-DCELLMAIN=OFF>)
target_compile_options(MolSim PRIVATE $<$<BOOL:${Strategie_2}>:-DPARTICLES=ON>) # -DSUBMAIN=OFF>)
unset(SUB CACHE)
unset(PAR CACHE)
unset(Strategy_1 CACHE)
unset(Strategy_2 CACHE)
