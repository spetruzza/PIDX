/*****************************************************
 **  PIDX Parallel I/O Library                      **
 **  Copyright (c) 2010-2014 University of Utah     **
 **  Scientific Computing and Imaging Institute     **
 **  72 S Central Campus Drive, Room 3750           **
 **  Salt Lake City, UT 84112                       **
 **                                                 **
 **  PIDX is licensed under the Creative Commons    **
 **  Attribution-NonCommercial-NoDerivatives 4.0    **
 **  International License. See LICENSE.md.         **
 **                                                 **
 **  For information about this project see:        **
 **  http://www.cedmav.com/pidx                     **
 **  or contact: pascucci@sci.utah.edu              **
 **  For support: PIDX-support@visus.net            **
 **                                                 **
 *****************************************************/

 /**
 * \file PIDX_comm.h
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * All communicator initialization related
 * functions.
 *
 */
 
#ifndef __PIDX_COMM_H
#define __PIDX_COMM_H 


#ifdef __cplusplus
extern "C" {
#endif

/// PIDX_access is neccessary to manage reading and writing by parallel MPI processes. It must
/// be created using PIDX_create_access() prior to opening or creating a PIDX file using
/// PIDX_file_create() or PIDX_file_open().
struct PIDX_access_struct
{
  int parallel;
  int idx_count[3];
  int sub_div[3];
  int rank_component[3];
  
  //int topology_aware_io;
#if PIDX_HAVE_MPI
  MPI_Comm comm;
#endif
};
typedef struct PIDX_access_struct* PIDX_access;


/// Call this before opening or creating a PIDX_file.
PIDX_return_code PIDX_create_access(PIDX_access* access);


///
PIDX_return_code PIDX_close_access(PIDX_access access);



#if PIDX_HAVE_MPI
///
PIDX_return_code PIDX_set_idx_count(PIDX_access access, int idx_count_x, int idx_count_y, int idx_count_z);


///
PIDX_return_code PIDX_set_mpi_access(PIDX_access access, MPI_Comm comm);


///
PIDX_return_code PIDX_set_process_extent(PIDX_access access, int sub_div_x, int sub_div_y, int sub_div_z);


///
PIDX_return_code PIDX_set_process_rank_decomposition(PIDX_access access, int rank_x, int rank_y, int rank_z);
#endif

#ifdef __cplusplus
} //extern C
#endif

#endif
