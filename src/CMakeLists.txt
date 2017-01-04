cmake_minimum_required (VERSION 3.7.0)
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
find_package(MPI REQUIRED)
#find_package (HDF5 NAMES hdf5 COMPONENTS C static)
find_package (HDF5) # Find non-cmake built HDF5
INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})
include_directories(/opt/local/include/)
set (LINK_LIBS ${LINK_LIBS} ${HDF5_C_STATIC_LIBRARY})



SET(H5FLAGS "-DH5_USE_16_API")
add_definitions(${H5FLAGS})



include_directories(include)
include_directories(${MPI_INCLUDE_PATH})
include_directories(${HDF5_INCLUDE_PATH})
#link_directories(${MPI_CXX_LIBRARIES} /opt/local/lib/openmpi-gcc48/)

add_executable(main main.cpp TwoDim.cpp                       parflux.cpp
 out.cpp parboundary.cpp parupdate.cpp parflux.cpp parsecond.cpp maxspeed.cpp riemann.cpp hlld.cpp sgn.cpp ctop.cpp lf.cpp outhdf5.cpp eigenvectors.cpp emf_exchange.cpp maes.cpp blast.cpp initialise_uniform.cpp cooling.cpp molcool.cpp tabfind.cpp locate.cpp orszagtang.cpp zanni.cpp zannisimple.cpp initialise_jet.cpp rungekutta.cpp minmod.cpp vanleer.cpp 

PhysConsts.h minmod.h  tnt_fortran_array1d.h TwoDim.h  molcool.h  tnt_fortran_array1d_utils.h boundary.h  orszagtang.h tnt_fortran_array2d.h cooling.h  out.h  tnt_fortran_array2d_utils.h ctop.h  output.h  tnt_fortran_array3d.h eigenvectors.h parboundary.h tnt_fortran_array3d_utils.h emf_exchange.h parflux.h  tnt_i_refvec.h falle.h  parupdate.h  tnt_lapack.h global.h  physics.h  tnt_linalg.h hlld.h  problem.h  tnt_math_utils.h initialise.h riemann.h  tnt_matrix.h initialise_blast.h roe.h  tnt_sparse_matrix.h initialise_jet.h rungekutta.h tnt_sparse_vector.h initialise_maes.h sgn.h  tnt_stopwatch.h initialise_uniform.h tabfind.h  tnt_subscript.h lf.h  tnt.h  tnt_vector.h locate.h  tnt_array1d.h tnt_version.h log.h  tnt_array1d_utils.h vanleer.h main.h  tnt_array2d.h vanleer_fvsplit.h maxspd.h  tnt_array2d_utils.h zanni.h maxspeed.h  tnt_array3d.h zannisimple.h mhd_godunov.h tnt_array3d_utils.h) 
target_include_directories (main PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(main ${MPI_CXX_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(main PROPERTIES
    COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(main PROPERTIES
    LINK_FLAGS "${MPI_CXX_LINK_FLAGS}" "/opt/local/lib")
endif()
