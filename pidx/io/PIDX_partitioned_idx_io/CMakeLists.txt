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
# ////////////////////////////////////////
# sources
# ////////////////////////////////////////

FILE(GLOB PIDX_PARTITIONED_IDX_IO_SOURCES *.h *.c)


# ////////////////////////////////////////
# includes
# ////////////////////////////////////////

IF (MPI_C_FOUND)
  MESSAGE("Configuring pidx with MPI support")
  INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
  SET(PIDX_PARTITIONED_IDX_IO_LINK_LIBS ${PIDX_PARTITIONED_IDX_IO_LINK_LIBS} ${MPI_C_LIBRARIES})
ENDIF()

IF (ZFP_FOUND)
  MESSAGE("Configuring pidx with ZFP support")
  INCLUDE_DIRECTORIES(${ZFP_INCLUDE_DIR})
  SET(PIDX_PARTITIONED_IDX_IO_LINK_LIBS ${PIDX_PARTITIONED_IDX_IO_LINK_LIBS} ${ZFP_LIBRARIES})
ENDIF()

# ////////////////////////////////////////
# library
# ////////////////////////////////////////

PIDX_ADD_LIBRARY(pidx_partitioned_idx_io "${PIDX_PARTITIONED_IDX_IO_SOURCES}")
TARGET_LINK_LIBRARIES(pidx_partitioned_idx_io pidx_utils pidx_data_handle m pidx_header pidx_rst pidx_block_rst pidx_cmp pidx_hz pidx_agg pidx_file_io pidx_comm pidx_topology ${PIDX_PARTITIONED_IDX_IO_LINK_LIBS})
SET_TARGET_PROPERTIES(pidx_partitioned_idx_io PROPERTIES LINKER_LANGUAGE CXX)
