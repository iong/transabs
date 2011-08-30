# - Try to find Eigen3
# Once done this will define
#  
#  Eigen3_FOUND        - system has OpenGL
#  Eigen3_INCLUDE_DIR  - the GL include directory

#=============================================================================
# Copyright 2011 Ionut Georgescu <george@pks.mpg.de>



FIND_PATH(Eigen3_INCLUDE_DIR Eigen/Eigen PATHS ${HOME}/sw/include/eigen3)

SET( Eigen3_FOUND "NO" )
IF(Eigen3_INCLUDE_DIR)
    SET( Eigen3_FOUND "YES" )
ENDIF(Eigen3_INCLUDE_DIR)

MARK_AS_ADVANCED(
  Eigen3_INCLUDE_DIR
)
