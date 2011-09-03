cmake_minimum_required(VERSION 2.8)

PROJECT(transabs)

include_directories(${PROJECT_SOURCE_DIR}/transabs)

#BOOST
set(Boost_USE_STATIC_LIBS   ON)
set(Boost_USE_MULTITHREADED ON)
find_package(Boost COMPONENTS program_options REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
else()
  message(FATAL_ERROR "BOOST _is_ needed!")
endif()

#Armadillo
FIND_PATH(ARMA_HOME include/armadillo PATHS $ENV{HOME}/sw/x86_64)
if(ARMA_HOME)
  include_directories(${ARMA_HOME}/include)
  link_directories(${ARMA_HOME}/lib64)
else()
  message(FATAL_ERROR "Armadillo _is_ needed!")
endif()


ADD_EXECUTABLE(transabs transabs/main.cpp transabs/DirectSum.cpp transabs/Utils.cpp)
target_link_libraries(transabs ${Boost_LIBRARIES} armadillo)

ADD_EXECUTABLE(t t.cc)