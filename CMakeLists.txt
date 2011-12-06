cmake_minimum_required(VERSION 2.8)

PROJECT(transabs)

include_directories(${PROJECT_SOURCE_DIR}/src)

# BOOST
set(Boost_USE_STATIC_LIBS   ON)
set(Boost_USE_MULTITHREADED OFF)
find_package(Boost COMPONENTS program_options filesystem system REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# HDF5
find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
link_directories(${HDF5_LIBRARY_DIRS})


# Armadillo
FIND_PATH(ARMA_HOME include/armadillo PATHS $ENV{HOME}/sw/x86_64 /opt/hpc)

if(ARMA_HOME)
  include_directories(${ARMA_HOME}/include)
  if (${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR} STREQUAL "Linux-x86_64")
	  link_directories(${ARMA_HOME}/lib64)
  else()
	  link_directories(${ARMA_HOME}/lib)
  endif()
else()
  message(FATAL_ERROR "Armadillo _is_ needed!")
endif()

# MPI
set(MPI_LIBRARIES)
find_package(MPI)
if (MPI_CXX_FOUND OR MPI_FOUND)
	include_directories(${MPI_INCLUDE_PATH})
	add_definitions(-DHAVE_MPI)
	set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
	set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
	message("You should run as: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} PROCS ${MPIEXEC_PREFLAGS} transabs ${MPIEXEC_POSTFLAGS} ARGS")
endif()


ADD_EXECUTABLE(transabs src/Statistics.cpp src/Electron.cpp src/Localization.cpp src/main.cpp src/DirectSum.cpp src/Utils.cpp src/Argon.h src/Atom.h src/Common.h src/DirectSum.h src/Electron.h src/HDF5IO.h src/Laser.h src/Localization.h src/Utils.h src/VectorMap.h)
target_link_libraries(transabs ${Boost_LIBRARIES} ${HDF5_LIBRARIES} ${MPI_LIBRARIES})

#ADD_EXECUTABLE(t t.cc)
