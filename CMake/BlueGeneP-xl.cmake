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

# the name of the target operating system
set(CMAKE_SYSTEM_NAME BlueGeneP-dynamic)
ADD_DEFINITIONS(-DBLUEGENEP)

# Set search paths to prefer local, admin-installed wrappers for the BG backend compilers
set(BGP_XL_COMPILER_SEARCH_PATHS /usr/local/bin /usr/bin)

# GNU C Compilers
find_program(CMAKE_C_COMPILER       bgxlc    ${BGP_XL_COMPILER_SEARCH_PATHS} /opt/ibmcmp/vac/bg/9.0/bin)
find_program(CMAKE_CXX_COMPILER     bgxlC    ${BGP_XL_COMPILER_SEARCH_PATHS} /opt/ibmcmp/vacpp/bg/9.0/bin)
find_program(CMAKE_Fortran_COMPILER bgxlf90  ${BGP_XL_COMPILER_SEARCH_PATHS} /opt/ibmcmp/xlf/bg/11.1/bin)

# Make sure MPI_COMPILER wrapper matches the gnu compilers.  
# Prefer local machine wrappers to driver wrappers here too.
find_program(MPI_COMPILER NAMES mpixlcxx mpixlc++ mpixlC mpixlc
  PATHS 
  /usr/local/bin
  /usr/bin
  /bgsys/drivers/ppcfloor/comm/bin
  /bgsys/drivers/ppcfloor/comm/default/bin)
