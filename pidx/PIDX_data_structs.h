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

#include "Generic_data_structs.h"
#include "PIDX_blocks.h"

#ifndef __PIDX_DATA_STRUCTS_H
#define __PIDX_DATA_STRUCTS_H

enum IO_MODE { PIDX_READ, PIDX_WRITE};

struct PIDX_variable_struct
{
  char* var_name;
  int bits_per_sample;
  int values_per_sample;
  int bits_per_value;
  char* type_name;
  int data_layout;

  int patch_count;
  Ndim_buffer patch[1024];
  HZ_buffer HZ_patch[1024];
  
  block_layout* global_block_layout;
  int *blocks_per_file;
  int existing_file_count;
  int *existing_file_index;
};
typedef struct PIDX_variable_struct* PIDX_variable;

struct idx_file_struct
{
  int current_time_step;
  int variable_count;
  int variable_index_tracker;
  
  PIDX_variable variable[1024];
  
  char* filename;
  int bits_per_block;
  int blocks_per_file;
  int* global_bounds;
  char bitSequence[512];
  char bitPattern[512];
};
typedef struct idx_file_struct* idx_dataset;

struct idx_dataset_derived_metadata_struct
{
  int *file_bitmap;
  int dimension;
  int samples_per_block;
  int maxh;
  int max_file_count;
};
typedef struct idx_dataset_derived_metadata_struct* idx_dataset_derived_metadata;

#endif
