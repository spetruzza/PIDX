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
 
#include "PIDX_agg.h"


struct PIDX_agg_struct 
{
#if PIDX_HAVE_MPI
  MPI_Comm comm;
#endif
  
  //Contains all relevant IDX file info
  //Blocks per file, samples per block, bitmask, box, file name template and more
  idx_dataset idx_ptr;
  
  //Contains all derieved IDX file info
  //number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_derived_ptr;
  
  int start_var_index;
  int end_var_index;
  
  int aggregator_interval;
  
  int ***rank_holder;
  MPI_Win win;
};

PIDX_agg_id PIDX_agg_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived_ptr, int start_var_index, int end_var_index)
{  
  PIDX_agg_id agg_id;

  agg_id = malloc(sizeof (*agg_id));
  memset(agg_id, 0, sizeof (*agg_id));

  agg_id->idx_ptr = (idx_dataset)malloc(sizeof(*(agg_id->idx_ptr)));
  memcpy(agg_id->idx_ptr, idx_meta_data, sizeof(*(agg_id->idx_ptr)));
  
  agg_id->idx_derived_ptr = (idx_dataset_derived_metadata)malloc(sizeof(*(agg_id->idx_derived_ptr)));
  memcpy(agg_id->idx_derived_ptr, idx_derived_ptr, sizeof(*(agg_id->idx_derived_ptr)));
  
  agg_id->start_var_index = start_var_index;
  agg_id->end_var_index = end_var_index;
  
  return agg_id;
}

#if PIDX_HAVE_MPI
int PIDX_agg_set_communicator(PIDX_agg_id agg_id, MPI_Comm comm)
{
  MPI_Comm_dup(comm, &agg_id->comm);
  return 0;
}
#endif

int aggregate_write_read(PIDX_agg_id agg_id, Agg_buffer agg_buffer, int variable_index, unsigned long long hz_start_index, unsigned long long hz_count, unsigned char* hz_buffer, int buffer_offset, int MODE)
{
  int rank = 0, nprocs = 1, itr;
  int bytes_per_datatype;
  int file_no = 0, block_no = 0, negative_block_offset = 0, sample_index = 0, values_per_sample;
  int target_rank = 0;
  //unsigned char *local_addr;
  long long start_agg_index = 0, end_agg_index = 0, target_disp = 0, target_count = 0, hz_start = 0, samples_in_file = 0;
  long long samples_per_file = (long long) agg_id->idx_derived_ptr->samples_per_block * agg_id->idx_ptr->blocks_per_file;

#if PIDX_HAVE_MPI
  MPI_Comm_rank(agg_id->comm, &rank);
  MPI_Comm_size(agg_id->comm, &nprocs);
#endif

  values_per_sample = agg_id->idx_ptr->variable[variable_index]->values_per_sample; //number of samples for variable j

  //starting HZ index for the data buffer at level "level" and for regular box number "box"
  hz_start = hz_start_index;
  
  //file number to which the first element of the buffer belongs to
  file_no = hz_start / samples_per_file;

  //block number for the first element of the buffer
  block_no = hz_start / agg_id->idx_derived_ptr->samples_per_block;

  //number of empty blocks befor block "block_no" in the file "file_no"
  negative_block_offset = find_block_negative_offset(agg_id->idx_ptr->blocks_per_file, block_no, agg_id->idx_ptr->variable[variable_index]->global_block_layout);
  assert(negative_block_offset >= 0);

  //number of samples in file "file_no"
  samples_in_file = agg_id->idx_ptr->variable[variable_index]->blocks_per_file[file_no] * agg_id->idx_derived_ptr->samples_per_block;
  assert(samples_in_file <= samples_per_file);

  //Calculating the hz index of "hz_start" relative to the file to which it belongs also taking into account empty blocks in file
  assert(hz_start >= (samples_per_file * file_no) + (negative_block_offset * agg_id->idx_derived_ptr->samples_per_block));
  target_disp = ((hz_start - ((samples_per_file * file_no) + (negative_block_offset * agg_id->idx_derived_ptr->samples_per_block))) * values_per_sample)
	  %
	  (samples_in_file * values_per_sample);
  assert(target_disp >= 0);

  sample_index = target_disp / (samples_in_file);
  assert(sample_index < agg_id->idx_ptr->variable[variable_index]->values_per_sample);
  
  target_disp = target_disp % (samples_in_file);
  
  target_rank = agg_id->rank_holder[variable_index - agg_id->start_var_index][sample_index][file_no];
  target_count = hz_count * values_per_sample;
  
  bytes_per_datatype = agg_id->idx_ptr->variable[variable_index]->bits_per_value / 8;
  hz_buffer = hz_buffer + buffer_offset * bytes_per_datatype * values_per_sample;
  
  start_agg_index = target_disp / (long long) (samples_in_file);
  end_agg_index = ((target_disp + target_count - 1) / (long long) (samples_in_file));
  assert(start_agg_index >= 0 && end_agg_index >= 0 && end_agg_index >= start_agg_index);
  
  if (start_agg_index != end_agg_index) 
  {
    if(target_rank != rank)
    {
#if PIDX_HAVE_MPI
#ifndef ACTIVE_TARGET
      MPI_Win_lock(MPI_LOCK_SHARED, target_rank, 0 , agg_id->win);
#endif
      if(MODE == PIDX_WRITE)
	MPI_Put(hz_buffer, (samples_in_file - target_disp) * bytes_per_datatype, MPI_BYTE, target_rank, target_disp, (samples_in_file - target_disp) * bytes_per_datatype, MPI_BYTE, agg_id->win);
      else
	MPI_Get(hz_buffer, (samples_in_file - target_disp) * bytes_per_datatype, MPI_BYTE, target_rank, target_disp, (samples_in_file - target_disp) * bytes_per_datatype, MPI_BYTE, agg_id->win);
#ifndef ACTIVE_TARGET
      MPI_Win_unlock(target_rank, agg_id->win);
#endif
#endif
    } 
    else
      if(MODE == PIDX_WRITE)
	memcpy( agg_buffer->buffer + target_disp * bytes_per_datatype, hz_buffer, (samples_in_file - target_disp) * bytes_per_datatype);
      else
	memcpy( hz_buffer, agg_buffer->buffer + target_disp * bytes_per_datatype, (samples_in_file - target_disp) * bytes_per_datatype);
      
      for (itr = 0; itr < end_agg_index - start_agg_index - 1; itr++) 
      {
	if(target_rank != rank)
	{
#if PIDX_HAVE_MPI
#ifndef ACTIVE_TARGET
	  MPI_Win_lock(MPI_LOCK_SHARED, target_rank + agg_id->aggregator_interval, 0, agg_id->win);
#endif
	  if(MODE == PIDX_WRITE)
	    MPI_Put(hz_buffer + ((samples_in_file - target_disp) + (itr * samples_in_file)) * bytes_per_datatype, samples_in_file * bytes_per_datatype, MPI_BYTE, target_rank + agg_id->aggregator_interval, 0, samples_in_file * bytes_per_datatype, MPI_BYTE, agg_id->win);
	  else
	    MPI_Get(hz_buffer + ((samples_in_file - target_disp) + (itr * samples_in_file)) * bytes_per_datatype, samples_in_file * bytes_per_datatype, MPI_BYTE, target_rank + agg_id->aggregator_interval, 0, samples_in_file * bytes_per_datatype, MPI_BYTE, agg_id->win);

#ifndef ACTIVE_TARGET
	  MPI_Win_unlock(target_rank + agg_id->aggregator_interval, agg_id->win);
#endif
#endif
	}
	else
	  if(MODE == PIDX_WRITE)
	    memcpy( agg_buffer->buffer, hz_buffer + ((samples_in_file - target_disp) + (itr * samples_in_file)) * bytes_per_datatype, samples_in_file * bytes_per_datatype);
	  else
	    memcpy( hz_buffer + ((samples_in_file - target_disp) + (itr * samples_in_file)) * bytes_per_datatype, agg_buffer->buffer, samples_in_file * bytes_per_datatype);
      }
      
      if(target_rank + agg_id->aggregator_interval != rank)
      {
#if PIDX_HAVE_MPI
#ifndef ACTIVE_TARGET
	MPI_Win_lock(MPI_LOCK_SHARED, target_rank + agg_id->aggregator_interval, 0, agg_id->win);
#endif
	if(MODE == PIDX_WRITE)
	  MPI_Put(hz_buffer + ((samples_in_file - target_disp) + ((end_agg_index - start_agg_index - 1) * samples_in_file)) * bytes_per_datatype, (target_count - (((end_agg_index - start_agg_index - 1) * (samples_in_file)) + ((samples_in_file) - target_disp))) * bytes_per_datatype, MPI_BYTE, target_rank + agg_id->aggregator_interval, 0, (target_count - ((end_agg_index - start_agg_index) * samples_in_file - target_disp)) * bytes_per_datatype, 
		MPI_BYTE, agg_id->win);
	else
	  MPI_Get(hz_buffer + ((samples_in_file - target_disp) + ((end_agg_index - start_agg_index - 1) * samples_in_file)) * bytes_per_datatype, (target_count - (((end_agg_index - start_agg_index - 1) * (samples_in_file)) + ((samples_in_file) - target_disp))) * bytes_per_datatype, MPI_BYTE, target_rank + agg_id->aggregator_interval, 0, (target_count - ((end_agg_index - start_agg_index) * samples_in_file - target_disp)) * bytes_per_datatype, 
		MPI_BYTE, agg_id->win);
#ifndef ACTIVE_TARGET
	MPI_Win_unlock(target_rank + agg_id->aggregator_interval, agg_id->win);
#endif
#endif
      }
      else
	if(MODE == PIDX_WRITE)
	  memcpy( agg_buffer->buffer, hz_buffer + ((samples_in_file - target_disp) + ((end_agg_index - start_agg_index - 1) * samples_in_file)) * bytes_per_datatype, (target_count - ((end_agg_index - start_agg_index) * samples_in_file - target_disp)) * bytes_per_datatype);    
	else
	  memcpy( hz_buffer + ((samples_in_file - target_disp) + ((end_agg_index - start_agg_index - 1) * samples_in_file)) * bytes_per_datatype, agg_buffer->buffer, (target_count - ((end_agg_index - start_agg_index) * samples_in_file - target_disp)) * bytes_per_datatype);
  }
  else 
  {
#if 0   
    if(rank == 1)
    {
      //printf("[%d] Count = %d x %d x %d\n", rank, hz_count, values_per_sample, bytes_per_datatype);
      int value;
      int u;
      for(u = 0; u < hz_count *values_per_sample; u++)
      {
	memcpy(&value, hz_buffer + u * bytes_per_datatype, bytes_per_datatype); 
	printf("Value at %d = %d\n", u, value);
      }
    }
#endif
    if(target_rank != rank)
    {
#if PIDX_HAVE_MPI
#ifndef ACTIVE_TARGET
      MPI_Win_lock(MPI_LOCK_SHARED, target_rank, 0 , agg_id->win);
#endif
      if(MODE == PIDX_WRITE)
	MPI_Put(hz_buffer, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, target_rank, target_disp, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, agg_id->win);
      else
	MPI_Get(hz_buffer, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, target_rank, target_disp, hz_count * values_per_sample * bytes_per_datatype, MPI_BYTE, agg_id->win);
#ifndef ACTIVE_TARGET
      MPI_Win_unlock(target_rank, agg_id->win);
#endif
#endif
    }
    else
    {
      if(MODE == PIDX_WRITE)
	memcpy( agg_buffer->buffer + target_disp * bytes_per_datatype, hz_buffer, hz_count * values_per_sample * bytes_per_datatype);
      else
	memcpy( hz_buffer, agg_buffer->buffer + target_disp * bytes_per_datatype, hz_count * values_per_sample * bytes_per_datatype);
#if 0
      //printf("Count = %d x %d x %d\n", hz_count, values_per_sample, bytes_per_datatype);
      double value;
      int u;
      for(u = 0; u < hz_count * values_per_sample; u++)
      {
	memcpy(&value, agg_buffer->buffer + (target_disp + u) * bytes_per_datatype, bytes_per_datatype); 
	printf("Value at %d %d = %f\n", target_disp, u, value);
      }
#endif
    }
  }
  
  return PIDX_success;
}

int PIDX_agg_aggregate(PIDX_agg_id agg_id, Agg_buffer agg_buffer) 
{ 
  int i, j, k, var;
  int rank_counter = 0, no_of_aggregators = 0, nprocs = 1, rank = 0;
  
#if PIDX_HAVE_MPI
  MPI_Comm_size(agg_id->comm, &nprocs);
  MPI_Comm_rank(agg_id->comm, &rank);
#endif
  
  for (var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
    no_of_aggregators = no_of_aggregators + agg_id->idx_ptr->variable[var]->values_per_sample * agg_id->idx_ptr->variable[var]->existing_file_count;
  
  agg_id->aggregator_interval = nprocs/no_of_aggregators;
  
  agg_buffer->buffer_size = 0;
  agg_buffer->sample_number = -1;
  agg_buffer->var_number = -1;
  agg_buffer->file_number = -1;

  agg_id->rank_holder = malloc((agg_id->end_var_index - agg_id->start_var_index + 1) * sizeof (int**));
  for (i = agg_id->start_var_index; i <= agg_id->end_var_index; i++) 
  {
    agg_id->rank_holder[i - agg_id->start_var_index] = malloc( agg_id->idx_ptr->variable[i]->values_per_sample  * sizeof (int*));
    for (j = 0; j < agg_id->idx_ptr->variable[i]->values_per_sample; j++)
    {
      agg_id->rank_holder[i - agg_id->start_var_index][j] = malloc(/*agg_id->idx_ptr->variable[i]->existing_file_count*/ agg_id->idx_derived_ptr->max_file_count * sizeof (int));
      memset(agg_id->rank_holder[i - agg_id->start_var_index][j], 0, agg_id->idx_derived_ptr->max_file_count * sizeof (int));
    }
  }
  
  rank_counter = 0;
  for (i = agg_id->start_var_index; i <= agg_id->end_var_index; i++)
  {
    for (j = 0; j < agg_id->idx_ptr->variable[i]->values_per_sample; j++)
    {
      for (k = 0; k < agg_id->idx_ptr->variable[i]->existing_file_count; k++)
      {
	agg_id->rank_holder[i - agg_id->start_var_index][j][agg_id->idx_ptr->variable[i]->existing_file_index[k]] = rank_counter;
	rank_counter = rank_counter + agg_id->aggregator_interval;
	
	if(rank == agg_id->rank_holder[i - agg_id->start_var_index][j][agg_id->idx_ptr->variable[i]->existing_file_index[k]])
	{
	  agg_buffer->file_number = agg_id->idx_ptr->variable[i]->existing_file_index[k];
	  agg_buffer->var_number = i;
	  agg_buffer->sample_number = j;
	  
	  agg_buffer->buffer_size = agg_id->idx_ptr->variable[agg_buffer->var_number]->blocks_per_file[agg_buffer->file_number] * agg_id->idx_derived_ptr->samples_per_block * (agg_id->idx_ptr->variable[agg_buffer->var_number]->bits_per_value/8);
	  agg_buffer->buffer = malloc(agg_buffer->buffer_size);
	  memset(agg_buffer->buffer, 0, agg_buffer->buffer_size);
	  //printf("Aggregator Rank %d Buffer Size %d (Var no: %d) (Sample no: %d) (File no: %d) (%d x %d x %d)\n", rank, agg_buffer->buffer_size, agg_buffer->var_number, agg_buffer->sample_number, agg_buffer->file_number, agg_id->idx_ptr->variable[agg_buffer->var_number]->blocks_per_file[agg_buffer->file_number], agg_id->idx_derived_ptr->samples_per_block, (agg_id->idx_ptr->variable[agg_buffer->var_number]->bits_per_value/8));
	}
      }
    }
  }
  
  return PIDX_success;
}

int PIDX_agg_aggregate_write_read(PIDX_agg_id agg_id, Agg_buffer agg_buffer, int MODE)
{
  int i, p, e1, var;
  int send_index = 0;
  long long index = 0, count = 0, hz_index = 0;
  
#if PIDX_HAVE_MPI
  if(agg_buffer->buffer_size != 0)
    MPI_Win_create(agg_buffer->buffer, agg_buffer->buffer_size, agg_id->idx_ptr->variable[agg_buffer->var_number]->bits_per_value/8, MPI_INFO_NULL, agg_id->comm, &(agg_id->win));
  else
    MPI_Win_create(0, 0, 1, MPI_INFO_NULL, agg_id->comm, &(agg_id->win));    
        
#ifdef ACTIVE_TARGET
  MPI_Win_fence(0, agg_id->win);
#else
  //MPI_Win_free has barrier semantics and therefore adding MPI_Barrier here is unnecessary
#endif
#endif
  
  for (p = 0; p < agg_id->idx_ptr->variable[agg_id->start_var_index]->patch_group_count; p++)
  {
    hz_index = 0, index = 0, count = 0, send_index = 0;
    if(agg_id->idx_ptr->variable[agg_id->start_var_index]->patch_group_ptr[p]->type == 0)
    {
      for (i = 0; i < agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_from; i++) 
	hz_index = hz_index + agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i];
      
      for (i = agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_from; i < agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_to; i++)
      {
	for(e1 = 0; e1 < agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] ; e1++)
	{
	  if(e1 == 0)
	  {
	    index = agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index];
	    send_index = e1;
	    count = 1;
	    
	    if(agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] == 1)
	    {
	      for(var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
	      {
		//printf("[A] Size %lld Offset %lld Send Index %d\n", count, index, send_index);
		aggregate_write_read(agg_id, agg_buffer, var, index, count, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], send_index, MODE);
	      }
	    }
	  }
	  else
	  {
	    if(agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index] - agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index - 1] == 1)
	    {
	      count++;
	      if(e1 == agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] - 1)
	      {
		for(var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
		{
		  //printf("[B] Size %lld Offset %lld Send Index %d\n", count, index, send_index);
		  aggregate_write_read(agg_id, agg_buffer, var, index, count, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], send_index, MODE);
		}
	      }
	    }
	    else
	    {
	      for(var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
	      {
		//printf("[C] Size %lld Offset %lld\n", count, index);
		aggregate_write_read(agg_id, agg_buffer, var, index, count, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], send_index, MODE);
	      }
	      
	      if(e1 == agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->samples_per_level[i] - 1)
	      {
		for(var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
		{
		  //printf("[D] Size %lld Offset %lld\n", count, index);
		  aggregate_write_read(agg_id, agg_buffer, var, agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index], 1, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], e1, MODE);
		}
	      }
	      index = agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->buffer_index[hz_index];
	      count = 1;
	      send_index = e1;
	    }
	  }
	  hz_index++;
	}
      }
    }
    
    else if(agg_id->idx_ptr->variable[agg_id->start_var_index]->patch_group_ptr[p]->type == 1)
    {
      for (i = agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_from; i < agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_to; i++)
      {
	for(var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
	{
	  index = 0;
	  count =  agg_id->idx_ptr->variable[var]->HZ_patch[p]->allign_end_hz[i] - agg_id->idx_ptr->variable[var]->HZ_patch[p]->allign_start_hz[i] + 1;
	  aggregate_write_read(agg_id, agg_buffer, var, agg_id->idx_ptr->variable[var]->HZ_patch[p]->allign_start_hz[i], count, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], 0, MODE);
	}
      }
    }
    
    else if(agg_id->idx_ptr->variable[agg_id->start_var_index]->patch_group_ptr[p]->type == 2)
    {
      int start_block_index, end_block_index, bl;
      for (i = agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_from; i < agg_id->idx_ptr->variable[agg_id->start_var_index]->HZ_patch[p]->HZ_level_to; i++)
      {
	for (var = agg_id->start_var_index; var <= agg_id->end_var_index; var++)
	{
	  start_block_index = agg_id->idx_ptr->variable[var]->HZ_patch[p]->allign_start_hz[i] / agg_id->idx_derived_ptr->samples_per_block;
	  end_block_index = agg_id->idx_ptr->variable[var]->HZ_patch[p]->allign_end_hz[i] / agg_id->idx_derived_ptr->samples_per_block;
	  assert(start_block_index >= 0 && end_block_index >= 0 && start_block_index <= end_block_index);
	  
	  send_index = 0;
	  for (bl = start_block_index; bl <= end_block_index; bl++) 
	  {
	    if (end_block_index == start_block_index) 
	    {
	      index = 0;
	      count = (agg_id->idx_ptr->variable[var]->HZ_patch[p]->allign_end_hz[i] - agg_id->idx_ptr->variable[var]->HZ_patch[p]->allign_start_hz[i] + 1);
	    } 
	    else
	    {
	      if (bl == start_block_index) 
	      {
		index = 0;
		count = ((start_block_index + 1) * agg_id->idx_derived_ptr->samples_per_block) - agg_id->idx_ptr->variable[var]->HZ_patch[p]->allign_start_hz[i];
	      } 
	      else if (bl == end_block_index) 
	      {
		index = (end_block_index * agg_id->idx_derived_ptr->samples_per_block - agg_id->idx_ptr->variable[var]->HZ_patch[p]->allign_start_hz[i]);
		count = agg_id->idx_ptr->variable[var]->HZ_patch[p]->allign_end_hz[i] - ((end_block_index) * agg_id->idx_derived_ptr->samples_per_block) + 1;
	      } 
	      else 
	      {
		index = (bl * agg_id->idx_derived_ptr->samples_per_block - agg_id->idx_ptr->variable[var]->HZ_patch[p]->allign_start_hz[i]);
		count = agg_id->idx_derived_ptr->samples_per_block;
	      }
	    }
	    aggregate_write_read(agg_id, agg_buffer, var, index + agg_id->idx_ptr->variable[var]->HZ_patch[p]->allign_start_hz[i], count, agg_id->idx_ptr->variable[var]->HZ_patch[p]->buffer[i], send_index, MODE);
	    
	    send_index = send_index + count;
	  }
	}
      }
    }
  }

#if PIDX_HAVE_MPI
#ifdef ACTIVE_TARGET
  MPI_Win_fence(0, agg_id->win);		//First Fence
#else
  // MPI_Win_create has barrier semantics and therefore adding MPI_Barrier here is unnecessary
#endif
  MPI_Win_free(&(agg_id->win));
#endif
  
#if 0
  if(agg_id->start_var_index == 0)
  {
    double value;
    for(i = 0; i < agg_buffer->buffer_size / sizeof(double); i++)
    {
      memcpy(&value, agg_buffer->buffer + i * sizeof(double), sizeof(double));
      printf("[%d] Value at %d = %f\n", rank, i, value);
    }
  }
#endif

  return PIDX_success;
}

int PIDX_agg_buf_destroy(Agg_buffer agg_buffer) 
{
  if (agg_buffer->buffer_size != 0) 
  {
    free(agg_buffer->buffer);
    agg_buffer->buffer = 0;
  }
  return PIDX_success;
}

int PIDX_agg_finalize(PIDX_agg_id agg_id) 
{
  int i = 0, j = 0;
  for (i = agg_id->start_var_index; i <= agg_id->end_var_index; i++) 
  {
    for (j = 0; j < agg_id->idx_ptr->variable[i]->values_per_sample; j++)
    {
      free(agg_id->rank_holder[i - agg_id->start_var_index][j]);
      agg_id->rank_holder[i - agg_id->start_var_index][j] = 0;
    }
    free(agg_id->rank_holder[i - agg_id->start_var_index]);
    agg_id->rank_holder[i - agg_id->start_var_index] = 0;
  }
  free(agg_id->rank_holder);
  agg_id->rank_holder = 0;

  free(agg_id->idx_ptr);
  agg_id->idx_ptr = 0;
  
  free(agg_id->idx_derived_ptr);
  agg_id->idx_derived_ptr = 0;
  
#if PIDX_HAVE_MPI
  MPI_Comm_free(&agg_id->comm);
#endif

  free(agg_id);
  agg_id = 0;

  return 0;
}