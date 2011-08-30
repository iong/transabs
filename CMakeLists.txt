cmake_minimum_required(VERSION 2.8)

PROJECT(transabs)


ADD_EXECUTABLE(transabs transabs/main.cpp transabs/Utils.cpp)
ADD_EXECUTABLE(t t.cc)
include_directories(${PROJECT_SOURCE_DIR}/transabs)

#BOOST
set(Boost_USE_STATIC_LIBS   ON)
set(Boost_USE_MULTITHREADED ON)
find_package(Boost COMPONENTS program_options REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(transabs ${Boost_LIBRARIES})
endif()

#Eigen3
FIND_PATH(Eigen3_INCLUDE_DIR Eigen/Eigen PATHS $ENV{HOME}/sw/include/eigen3)
if(Eigen3_INCLUDE_DIR)
  include_directories(${Eigen3_INCLUDE_DIR})
else()
  message(FATAL_ERROR "Eigen3 _is_ needed!")
endif()