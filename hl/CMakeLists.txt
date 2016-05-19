#/*****************************************************
# **  PIDX Parallel I/O Library                      **
# **  Copyright (c) 2010-2014 University of Utah     **
# **  Scientific Computing and Imaging Institute     **
# **  72 S Central Campus Drive, Room 3750           **
# **  Salt Lake City, UT 84112                       **
# **                                                 **
# **  PIDX is licensed under the Creative Commons    **
# **  Attribution-NonCommercial-NoDerivatives 4.0    **
# **  International License. See LICENSE.md.         **
# **                                                 **
# **  For information about this project see:        **
# **  http://www.cedmav.com/pidx                     **
# **  or contact: pascucci@sci.utah.edu              **
# **  For support: PIDX-support@visus.net            **
# **                                                 **
# *****************************************************/


IF (PIDX_BUILD_HL)

  # ////////////////////////////////////////
  # sources includes and link libraries
  # ////////////////////////////////////////

  INCLUDE_DIRECTORIES(./)
  INCLUDE_DIRECTORIES(../pidx)
  INCLUDE_DIRECTORIES(../)
  INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})

  FILE(GLOB SRCS FALSE ./*.cpp ./*.h)

  SET(PIDX_HL_LIB "pidx_hl")
  PIDX_ADD_LIBRARY(${PIDX_HL_LIB} ${SRCS})

  IF (PIDX_BUILD_METADATA)
    TARGET_LINK_LIBRARIES(${PIDX_HL_LIB} pidx pidx_metadata ${MPI_C_LIBRARIES})
  ELSE()
    TARGET_LINK_LIBRARIES(${PIDX_HL_LIB} pidx ${MPI_C_LIBRARIES})
  ENDIF()

  SET_TARGET_PROPERTIES(${PIDX_HL_LIB} PROPERTIES LINKER_LANGUAGE CXX)

  PIDX_ADD_EXECUTABLE(pidx_hl_test "PIDX_hl_test.cpp")
  TARGET_LINK_LIBRARIES(pidx_hl_test ${PIDX_HL_LIB})
  
ENDIF()
