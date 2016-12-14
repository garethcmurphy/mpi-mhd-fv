cmake_minimum_required(VERSION 3.6.2)
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
find_package(MPI REQUIRED)
#find_package (HDF5 NAMES hdf5 COMPONENTS C static)
find_package(HDF5) # Find non-cmake built HDF5
INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
include_directories(/opt/local/include/)
set(LINK_LIBS ${LINK_LIBS} ${HDF5_C_STATIC_LIBRARY})


SET(H5FLAGS "-DH5_USE_16_API")
add_definitions(${H5FLAGS})


include_directories(${MPI_INCLUDE_PATH})
include_directories(${HDF5_INCLUDE_PATH})
#link_directories(${MPI_CXX_LIBRARIES} /opt/local/lib/openmpi-gcc48/)

add_executable(main main.cpp TwoDim.cpp parflux.cpp
        out.cpp parboundary.cpp parupdate.cpp parflux.cpp parsecond.cpp
        maxspeed.cpp riemann.cpp hlld.cpp sgn.cpp ctop.cpp lf.cpp outhdf5.cpp
        eigenvectors.cpp emf_exchange.cpp maes.cpp blast.cpp
        initialise_uniform.cpp cooling.cpp molcool.cpp tabfind.cpp locate.cpp
        orszagtang.cpp zanni.cpp zannisimple.cpp initialise_jet.cpp rungekutta.cpp minmod.cpp vanleer.cpp
        )

target_include_directories(main PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(main ${MPI_CXX_LIBRARIES} /opt/local/lib/libhdf5.dylib
        /opt/local/lib/libz.dylib)

#if (MPI_COMPILE_FLAGS)
#    set_target_properties(main PROPERTIES
#            COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
#endif ()
#
#if (MPI_LINK_FLAGS)
#    set_target_properties(main PROPERTIES
#            LINK_FLAGS "${MPI_CXX_LINK_FLAGS}" "/opt/local/lib")
#endif ()
