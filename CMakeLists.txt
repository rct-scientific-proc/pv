cmake_minimum_required(VERSION 3.15)
project(phase_visualize C)

# ---- Options --------------------------------------------------------------
option(PV_BUILD_SHARED "Build phase_visualize as a shared library" ON)
option(PV_USE_OPENMP   "Enable OpenMP parallelization"             ON)

# ---- Standards ------------------------------------------------------------
set(CMAKE_C_STANDARD          11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS        OFF)

# Hide all symbols by default; PV_API explicitly exports what we want.
set(CMAKE_C_VISIBILITY_PRESET     hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)

# Place all build artifacts in one place so MATLAB can find the .dll/.so/.dylib.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# ---- Library --------------------------------------------------------------
if(PV_BUILD_SHARED)
    add_library(phase_visualize SHARED phase_visualize.c phase_visualize.h)
    target_compile_definitions(phase_visualize PRIVATE PHASE_VISUALIZE_BUILD_DLL)
else()
    add_library(phase_visualize STATIC phase_visualize.c phase_visualize.h)
endif()

target_include_directories(phase_visualize
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

# Link math library on Unix-like systems.
if(NOT WIN32)
    target_link_libraries(phase_visualize PRIVATE m)
endif()

# ---- OpenMP ---------------------------------------------------------------
if(PV_USE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_C_FOUND)
        target_link_libraries(phase_visualize PUBLIC OpenMP::OpenMP_C)
        message(STATUS "OpenMP enabled")
    else()
        message(WARNING "OpenMP requested but not found - building without it")
    endif()
endif()

# ---- Compiler flags -------------------------------------------------------
if(MSVC)
    target_compile_options(phase_visualize PRIVATE /W3 /O2)
else()
    target_compile_options(phase_visualize PRIVATE -Wall -Wextra -O3)
endif()

# ---- Install --------------------------------------------------------------
include(GNUInstallDirs)
install(TARGETS phase_visualize
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES phase_visualize.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

message(STATUS "phase_visualize will be built as "
        "$<IF:$<BOOL:${PV_BUILD_SHARED}>,SHARED,STATIC>")
