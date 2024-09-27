# Enable FetchContent
include(FetchContent)

# Find Python
find_package(Python COMPONENTS Interpreter Development QUIET)

if(Python_FOUND)
    ## Add pybind11
    # Find existing pybind11 installation
    find_package(pybind11 QUIET)
    # Download pybind11 if an existing installation is not found
    if(NOT ${pybind11_FOUND})
        # Fetch pybind11
        FetchContent_Declare(pybind11
            GIT_REPOSITORY https://github.com/pybind/pybind11.git
            GIT_TAG v2.13.6
        )

        # Make pybind11 available
        FetchContent_MakeAvailable(pybind11)
    endif(NOT ${pybind11_FOUND})

    # Add Eigen
    find_package(Eigen3 3.4 QUIET NO_MODULE)

    # Collect header files
    set(src_headers ${CMAKE_CURRENT_SOURCE_DIR}/include)
    file(GLOB_RECURSE src_hpp CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp)

    # Collect source files
    file(GLOB_RECURSE src_cpp CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

    # Build Python interface
    pybind11_add_module(pythalassa ${src_cpp} ${src_hpp})

    # Include header directory
    target_include_directories(pythalassa PUBLIC ${src_headers})

    # Link PyTHALASSA to CTHALASSA
    target_link_libraries(pythalassa PUBLIC cthalassa)

    # Link PyTHALASSA to Eigen
    target_link_libraries(pythalassa PUBLIC Eigen3::Eigen)
endif(Python_FOUND)
