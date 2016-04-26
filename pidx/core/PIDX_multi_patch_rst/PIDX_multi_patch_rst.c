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
 * \file PIDX_rst.c
 *
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions 
 * declared in PIDX_multi_patch_rst.h
 *
 */

#include "../../PIDX_inc.h"

//Struct for restructuring ID
struct PIDX_multi_patch_rst_struct
{
  //Passed by PIDX API
#if PIDX_HAVE_MPI
  MPI_Comm comm; //Communicator
#endif

  //Contains all relevant IDX file info
  //Blocks per file, samples per block, bitmask, patch, file name template and more
  idx_dataset idx;
  
  //Contains all derieved IDX file info
  //number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_derived;
  
  int init_index;
  int first_index;
  int last_index;
  
  //int if_perform_rst;

  //dimension of the power-two volume imposed patch
  int64_t reg_patch_size[PIDX_MAX_DIMENSIONS];
  int reg_patch_grp_count;
  Ndim_patch_group* reg_patch_grp;
  
  int64_t sim_max_patch_group_count;
  int64_t* sim_multi_patch_r_count;
  int64_t* sim_multi_patch_r_offset;  
  
};

static int maximum_neighbor_count = 256;

#if PIDX_HAVE_MPI
static int intersectNDChunk(Ndim_patch A, Ndim_patch B);
static int getPowerOftwo(int x);
#endif


#if PIDX_HAVE_MPI
/// Function to check if NDimensional data chunks A and B intersects
static int intersectNDChunk(Ndim_patch A, Ndim_patch B)
{
  int d = 0, check_bit = 0;
  for (d = 0; d < /*PIDX_MAX_DIMENSIONS*/3; d++)
    check_bit = check_bit || (A->offset[d] + A->size[d] - 1) < B->offset[d] || (B->offset[d] + B->size[d] - 1) < A->offset[d];
  
  return !(check_bit);
}
#endif


PIDX_multi_patch_rst_id PIDX_multi_patch_rst_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_derived, int first_index, int var_start_index, int var_end_index)
{
  //Creating the restructuring ID
  PIDX_multi_patch_rst_id multi_patch_rst_id;
  multi_patch_rst_id = (PIDX_multi_patch_rst_id)malloc(sizeof (*multi_patch_rst_id));
  memset(multi_patch_rst_id, 0, sizeof (*multi_patch_rst_id));

  multi_patch_rst_id->idx = idx_meta_data;
  multi_patch_rst_id->idx_derived = idx_derived;

  multi_patch_rst_id->init_index = first_index;
  multi_patch_rst_id->first_index = var_start_index;
  multi_patch_rst_id->last_index = var_end_index;

  return (multi_patch_rst_id);
}


#if PIDX_HAVE_MPI
PIDX_return_code PIDX_multi_patch_rst_set_communicator(PIDX_multi_patch_rst_id multi_patch_rst_id, MPI_Comm comm)
{
  if (multi_patch_rst_id == NULL)
    return PIDX_err_id;

  multi_patch_rst_id->comm = comm;

  return PIDX_success;
}
#endif


PIDX_return_code PIDX_multi_patch_rst_meta_data_create(PIDX_multi_patch_rst_id rst_id)
{
  PIDX_variable var0 = rst_id->idx->variable[rst_id->first_index];
  int p = 0, v = 0, j = 0;

#if PIDX_HAVE_MPI
  int r, d, c, nprocs, rank;
  int64_t i, k, l, m, max_vol, patch_count, pc;
  int reg_patch_count, edge_case = 0;
  
  if (rst_id->idx->enable_rst == 0)
    var0->patch_group_count = var0->sim_patch_count;
  else
  {
    MPI_Comm_rank(rst_id->comm, &rank);
    MPI_Comm_size(rst_id->comm, &nprocs);
    
    int start_var_index = rst_id->first_index;
    
    MPI_Allreduce(&rst_id->idx->variable[start_var_index]->sim_patch_count, &rst_id->sim_max_patch_group_count, 1, MPI_INT32_T, MPI_MAX, rst_id->comm);
   
    if (rank == 0)
      printf("loc %d max_patch_group_count %lld\n", rst_id->idx->variable[start_var_index]->sim_patch_count, rst_id->sim_max_patch_group_count);
    
    rst_id->sim_multi_patch_r_count = malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_group_count);
    memset(rst_id->sim_multi_patch_r_count, -1, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_group_count));
    rst_id->sim_multi_patch_r_offset = malloc(sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_group_count);
    memset(rst_id->sim_multi_patch_r_offset, -1, (sizeof (int64_t) * nprocs * PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_group_count));
    
    for(pc=0; pc < rst_id->idx->variable[start_var_index]->sim_patch_count; pc++){
      
      int64_t* tempoff = rst_id->idx->variable[start_var_index]->sim_patch[pc]->offset;
      int64_t* tempsize = rst_id->idx->variable[start_var_index]->sim_patch[pc]->size;
      
      printf("%d:%lld off %lld %lld %lld size %lld %lld %lld\n", rank, pc, tempoff[0],tempoff[1],tempoff[2],tempsize[0], tempsize[1],tempsize[2]);
      
      int64_t index = rank * (PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_group_count) + pc*PIDX_MAX_DIMENSIONS;
      int64_t* curr_patch_offset = &rst_id->sim_multi_patch_r_offset[index];
      int64_t* curr_patch_size = &rst_id->sim_multi_patch_r_count[index];
      
      memcpy(curr_patch_offset, tempoff,sizeof(int64_t) * PIDX_MAX_DIMENSIONS);
      memcpy(curr_patch_size, tempsize,sizeof(int64_t) * PIDX_MAX_DIMENSIONS);
      
      //     printf("%d:%lld off %lld %lld %lld size %lld %lld %lld\n", rank, pc, curr_patch_offset[0],curr_patch_offset[1],curr_patch_offset[2],curr_patch_size[0], curr_patch_size[1],curr_patch_size[2]);
    }
    
    MPI_Allgather(&rst_id->sim_multi_patch_r_count[rank * PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_group_count], PIDX_MAX_DIMENSIONS*rst_id->sim_max_patch_group_count, MPI_LONG_LONG, rst_id->sim_multi_patch_r_count, PIDX_MAX_DIMENSIONS*rst_id->sim_max_patch_group_count, MPI_LONG_LONG, rst_id->comm);
    
    MPI_Allgather(&rst_id->sim_multi_patch_r_offset[rank * PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_group_count], PIDX_MAX_DIMENSIONS*rst_id->sim_max_patch_group_count, MPI_LONG_LONG, rst_id->sim_multi_patch_r_offset, PIDX_MAX_DIMENSIONS*rst_id->sim_max_patch_group_count, MPI_LONG_LONG, rst_id->comm);
    
//      for(int r=0; r<nprocs; r++){
//        for(int pc=0; pc < rst_id->sim_max_patch_group_count; pc++){
//  
//          int64_t index = r * (PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_group_count) + pc*PIDX_MAX_DIMENSIONS;
//       //   int64_t* curr_patch_size = &file->sim_multi_patch_r_count[index];
//          int64_t* curr_patch_offset = &rst_id->sim_multi_patch_r_offset[index];
//  
//          if(curr_patch_offset[0] == -1)
//            printf("patch %d for rank %d doesn't exist\n", pc, r);
//          else if(rank == 2)
//            printf("patch %d for rank %d off %lld %lld %lld\n", pc, r,curr_patch_offset[0],curr_patch_offset[1],curr_patch_offset[2]);
//        }
//        
//      }
    
    //
    
    var0->patch_group_count = 0;
    
    /// STEP 1 : Compute the dimension of the regular patch
    //if (rst_id->idx->reg_patch_size[0] == 0)
    //  set_default_patch_size(rst_id, rst_id->idx_derived->rank_r_count, nprocs);
    //else
    //  memcpy(rst_id->reg_patch_size, rst_id->idx->reg_patch_size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
    
    memcpy(rst_id->reg_patch_size, rst_id->idx->reg_patch_size, sizeof(uint64_t) * PIDX_MAX_DIMENSIONS);
    
    /// extents for the local process(rank)
    Ndim_patch local_proc_patch = (Ndim_patch)malloc(sizeof (*local_proc_patch));
    memset(local_proc_patch, 0, sizeof (*local_proc_patch));
    for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    {
      local_proc_patch->offset[d] = rst_id->idx_derived->rank_r_offset[PIDX_MAX_DIMENSIONS * rank + d];
      local_proc_patch->size[d] = rst_id->idx_derived->rank_r_count[PIDX_MAX_DIMENSIONS * rank + d];
    }
    
    printf("%d: local off %lld %lld %lld local size %lld %lld %lld\n",rank, local_proc_patch->offset[0],local_proc_patch->offset[1],local_proc_patch->offset[2], local_proc_patch->size[0], local_proc_patch->size[1], local_proc_patch->size[2]);
    
    printf("local reg patch size %lld %lld %lld\n", rst_id->reg_patch_size[0], rst_id->reg_patch_size[1], rst_id->reg_patch_size[2]);
    
    int64_t adjusted_bounds[PIDX_MAX_DIMENSIONS];
    memcpy(adjusted_bounds, rst_id->idx->bounds, PIDX_MAX_DIMENSIONS * sizeof(unsigned long long));
    
    for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
    {
      adjusted_bounds[d] = rst_id->idx->bounds[d];
      if (rst_id->idx->bounds[d] % rst_id->idx->chunk_size[d] != 0)
        adjusted_bounds[d] = ((rst_id->idx->bounds[d] / rst_id->idx->chunk_size[d]) + 1) * rst_id->idx->chunk_size[d];
    }
    
    rst_id->reg_patch_grp_count = 0;
    for (i = 0; i < adjusted_bounds[0]; i = i + rst_id->reg_patch_size[0])
      for (j = 0; j < adjusted_bounds[1]; j = j + rst_id->reg_patch_size[1])
        for (k = 0; k < adjusted_bounds[2]; k = k + rst_id->reg_patch_size[2])
          for (l = 0; l < adjusted_bounds[3]; l = l + rst_id->reg_patch_size[3])
            for (m = 0; m < adjusted_bounds[4]; m = m + rst_id->reg_patch_size[4])
            {
              Ndim_patch reg_patch = (Ndim_patch)malloc(sizeof (*reg_patch));
              memset(reg_patch, 0, sizeof (*reg_patch));
              
              //Interior regular patches
              reg_patch->offset[0] = i;
              reg_patch->offset[1] = j;
              reg_patch->offset[2] = k;
              reg_patch->offset[3] = l;
              reg_patch->offset[4] = m;
              reg_patch->size[0] = rst_id->reg_patch_size[0];
              reg_patch->size[1] = rst_id->reg_patch_size[1];
              reg_patch->size[2] = rst_id->reg_patch_size[2];
              reg_patch->size[3] = rst_id->reg_patch_size[3];
              reg_patch->size[4] = rst_id->reg_patch_size[4];
              
              //Edge regular patches
              if ((i + rst_id->reg_patch_size[0]) > adjusted_bounds[0])
                reg_patch->size[0] = adjusted_bounds[0] - i;
              if ((j + rst_id->reg_patch_size[1]) > adjusted_bounds[1])
                reg_patch->size[1] = adjusted_bounds[1] - j;
              if ((k + rst_id->reg_patch_size[2]) > adjusted_bounds[2])
                reg_patch->size[2] = adjusted_bounds[2] - k;
              if ((l + rst_id->reg_patch_size[3]) > adjusted_bounds[3])
                reg_patch->size[3] = adjusted_bounds[3] - l;
              if ((m + rst_id->reg_patch_size[4]) > adjusted_bounds[4])
                reg_patch->size[4] = adjusted_bounds[4] - m;
              
              if (intersectNDChunk(reg_patch, local_proc_patch))
                rst_id->reg_patch_grp_count++;
              
              free(reg_patch);
            }
    
    rst_id->reg_patch_grp = (Ndim_patch_group*)malloc(sizeof(*rst_id->reg_patch_grp) * rst_id->reg_patch_grp_count);
    memset(rst_id->reg_patch_grp, 0, sizeof(*rst_id->reg_patch_grp) * rst_id->reg_patch_grp_count);
    
    //  printf("rst_id->reg_patch_grp_count = %d\n", rst_id->reg_patch_grp_count);
    
    reg_patch_count = 0;
    /// STEP 3 : iterate through extents of all imposed regular patches, and find all the regular patches a process (local_proc_patch) intersects with
    
    for (i = 0; i < adjusted_bounds[0]; i = i + rst_id->reg_patch_size[0])
      for (j = 0; j < adjusted_bounds[1]; j = j + rst_id->reg_patch_size[1])
        for (k = 0; k < adjusted_bounds[2]; k = k + rst_id->reg_patch_size[2])
          for (l = 0; l < adjusted_bounds[3]; l = l + rst_id->reg_patch_size[3])
            for (m = 0; m < adjusted_bounds[4]; m = m + rst_id->reg_patch_size[4])
            {
              Ndim_patch reg_patch = (Ndim_patch)malloc(sizeof (*reg_patch));
              memset(reg_patch, 0, sizeof (*reg_patch));
              
              //Interior regular patches
              reg_patch->offset[0] = i;
              reg_patch->offset[1] = j;
              reg_patch->offset[2] = k;
              reg_patch->offset[3] = l;
              reg_patch->offset[4] = m;
              reg_patch->size[0] = rst_id->reg_patch_size[0];
              reg_patch->size[1] = rst_id->reg_patch_size[1];
              reg_patch->size[2] = rst_id->reg_patch_size[2];
              reg_patch->size[3] = rst_id->reg_patch_size[3];
              reg_patch->size[4] = rst_id->reg_patch_size[4];
              
              //Edge regular patches
              edge_case = 0;
              if ((i + rst_id->reg_patch_size[0]) > adjusted_bounds[0])
              {
                reg_patch->size[0] = adjusted_bounds[0] - i;
                edge_case = 1;
              }
              if ((j + rst_id->reg_patch_size[1]) > adjusted_bounds[1])
              {
                reg_patch->size[1] = adjusted_bounds[1] - j;
                edge_case = 1;
              }
              if ((k + rst_id->reg_patch_size[2]) > adjusted_bounds[2])
              {
                reg_patch->size[2] = adjusted_bounds[2] - k;
                edge_case = 1;
              }
              if ((l + rst_id->reg_patch_size[3]) > adjusted_bounds[3])
              {
                reg_patch->size[3] = adjusted_bounds[3] - l;
                edge_case = 1;
              }
              if ((m + rst_id->reg_patch_size[4]) > adjusted_bounds[4])
              {
                reg_patch->size[4] = adjusted_bounds[4] - m;
                edge_case = 1;
              }
              
              /// STEP 4: If local process intersects with regular patch, then find all other process that intersects with the regular patch.
              if (intersectNDChunk(reg_patch, local_proc_patch))
              {
                //if (rank == 52 && reg_patch->offset[0] == 0 && reg_patch->offset[1] == 768 && reg_patch->offset[2] == 128)
                //printf("[g] reg box %d %d %d : %d %d %d local box %d %d %d : %d %d %d\n", reg_patch->offset[0], reg_patch->offset[1], reg_patch->offset[2], reg_patch->size[0], reg_patch->size[1], reg_patch->size[2], local_proc_patch->offset[0], local_proc_patch->offset[1], local_proc_patch->offset[2], local_proc_patch->size[0], local_proc_patch->size[1], local_proc_patch->size[2]);
                
                rst_id->reg_patch_grp[reg_patch_count] = malloc(sizeof(*(rst_id->reg_patch_grp[reg_patch_count])));
                memset(rst_id->reg_patch_grp[reg_patch_count], 0, sizeof(*(rst_id->reg_patch_grp[reg_patch_count])));
                
                Ndim_patch_group patch_grp = rst_id->reg_patch_grp[reg_patch_count];
                
                patch_grp->source_patch_rank = (int*)malloc(sizeof(int) * maximum_neighbor_count);
                patch_grp->patch = malloc(sizeof(*patch_grp->patch) * maximum_neighbor_count);
                patch_grp->reg_patch = malloc(sizeof(*patch_grp->reg_patch));
                memset(patch_grp->source_patch_rank, 0, sizeof(int) * maximum_neighbor_count);
                memset(patch_grp->patch, 0, sizeof(*patch_grp->patch) * maximum_neighbor_count);
                memset(patch_grp->reg_patch, 0, sizeof(*patch_grp->reg_patch));
                
                patch_count = 0;
                patch_grp->count = 0;
                if(edge_case == 0)
                  patch_grp->type = 1;
                else
                  patch_grp->type = 2;
                
                //Iterate through all processes
                for (r = 0; r < nprocs; r++)
                {
                  for(pc = 0; pc < rst_id->sim_max_patch_group_count; pc++)
                  {
                    //Extent of process with rank r
                    Ndim_patch curr_patch = malloc(sizeof (*curr_patch));
                    memset(curr_patch, 0, sizeof (*curr_patch));
                    
                    for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
                    {
                      int64_t index = r * (PIDX_MAX_DIMENSIONS * rst_id->sim_max_patch_group_count) + pc*PIDX_MAX_DIMENSIONS;
                      
                      curr_patch->offset[d] = rst_id->sim_multi_patch_r_offset[index+d];
                      curr_patch->size[d] = rst_id->sim_multi_patch_r_count[index+d];
   
//                      curr_patch->offset[d] = rst_id->idx_derived->rank_r_offset[PIDX_MAX_DIMENSIONS * r + d];
//                      curr_patch->size[d] = rst_id->idx_derived->rank_r_count[PIDX_MAX_DIMENSIONS * r + d];
                    }
                    
                    //If process with rank r intersects with the regular patch, then calculate the offset, count and volume of the intersecting volume
                    //if (rank == 52 && reg_patch->offset[0] == 0 && reg_patch->offset[1] == 768 && reg_patch->offset[2] == 128)
                    printf("[l %d] reg box %lld %lld %lld : %lld %lld %lld local box %lld %lld %lld : %lld %lld %lld\n", r, reg_patch->offset[0], reg_patch->offset[1], reg_patch->offset[2], reg_patch->size[0], reg_patch->size[1], reg_patch->size[2], curr_patch->offset[0], curr_patch->offset[1], curr_patch->offset[2], curr_patch->size[0], curr_patch->size[1], curr_patch->size[2]);
                    
                    if (intersectNDChunk(reg_patch, curr_patch))
                    {
                      //if (rank == 52 && reg_patch->offset[0] == 0 && reg_patch->offset[1] == 768 && reg_patch->offset[2] == 128)
                      printf("[li] reg box %lld %lld %lld : %lld %lld %lld local box %lld %lld %lld : %lld %lld %lld\n", reg_patch->offset[0], reg_patch->offset[1], reg_patch->offset[2], reg_patch->size[0], reg_patch->size[1], reg_patch->size[2], curr_patch->offset[0], curr_patch->offset[1], curr_patch->offset[2], curr_patch->size[0], curr_patch->size[1], curr_patch->size[2]);
                      
                      patch_grp->patch[patch_count] = malloc(sizeof(*(patch_grp->patch[patch_count])));
                      memset(patch_grp->patch[patch_count], 0, sizeof(*(patch_grp->patch[patch_count])));
                      
                      for (d = 0; d < PIDX_MAX_DIMENSIONS; d++)
                      {
                        //STEP 5 : offset and count of intersecting chunk of process with rank r and regular patch
                        if (curr_patch->offset[d] <= reg_patch->offset[d] && (curr_patch->offset[d] + curr_patch->size[d] - 1) <= (reg_patch->offset[d] + reg_patch->size[d] - 1))
                        {
                          patch_grp->patch[patch_count]->offset[d] = reg_patch->offset[d];
                          patch_grp->patch[patch_count]->size[d] = (curr_patch->offset[d] + curr_patch->size[d] - 1) - reg_patch->offset[d] + 1;
                        }
                        else if (reg_patch->offset[d] <= curr_patch->offset[d] && (curr_patch->offset[d] + curr_patch->size[d] - 1) >= (reg_patch->offset[d] + reg_patch->size[d] - 1))
                        {
                          patch_grp->patch[patch_count]->offset[d] = curr_patch->offset[d];
                          patch_grp->patch[patch_count]->size[d] = (reg_patch->offset[d] + reg_patch->size[d] - 1) - curr_patch->offset[d] + 1;
                        }
                        else if (( reg_patch->offset[d] + reg_patch->size[d] - 1) <= (curr_patch->offset[d] + curr_patch->size[d] - 1) && reg_patch->offset[d] >= curr_patch->offset[d])
                        {
                          patch_grp->patch[patch_count]->offset[d] = reg_patch->offset[d];
                          patch_grp->patch[patch_count]->size[d] = reg_patch->size[d];
                        }
                        else if (( curr_patch->offset[d] + curr_patch->size[d] - 1) <= (reg_patch->offset[d] + reg_patch->size[d] - 1) && curr_patch->offset[d] >= reg_patch->offset[d])
                        {
                          patch_grp->patch[patch_count]->offset[d] = curr_patch->offset[d];
                          patch_grp->patch[patch_count]->size[d] = curr_patch->size[d];
                        }
                        
                        //offset and count of intersecting regular patch
                        patch_grp->reg_patch->offset[d] = reg_patch->offset[d];
                        patch_grp->reg_patch->size[d] = reg_patch->size[d];
                      }
                      //if (rank == 52 && reg_patch->offset[0] == 0 && reg_patch->offset[1] == 768 && reg_patch->offset[2] == 128)
                       //printf("[%d] oc : %lld %lld %lld :: %lld %lld %lld\n", r, patch_grp->patch[patch_count]->offset[0], patch_grp->patch[patch_count]->offset[1], patch_grp->patch[patch_count]->offset[2], patch_grp->patch[patch_count]->size[0], patch_grp->patch[patch_count]->size[1], patch_grp->patch[patch_count]->size[2]);
                      
                      patch_grp->source_patch_rank[patch_count] = r;
                      patch_count++;
                      
                      if (patch_count >= maximum_neighbor_count)
                      {
                        maximum_neighbor_count = maximum_neighbor_count * 2;
                        
                        int *temp_buffer2 = realloc(patch_grp->source_patch_rank, maximum_neighbor_count * sizeof(int));
                        if (temp_buffer2 == NULL)
                        {
                          fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
                          return PIDX_err_rst;
                        }
                        else
                          patch_grp->source_patch_rank = temp_buffer2;
                        
                        Ndim_patch *temp_buffer3 = realloc(patch_grp->patch, maximum_neighbor_count * sizeof(*patch_grp->patch));
                        if (temp_buffer3 == NULL)
                        {
                          fprintf(stderr, "[%s] [%d] realloc() failed.\n", __FILE__, __LINE__);
                          return PIDX_err_rst;
                        }
                        else
                          patch_grp->patch = temp_buffer3;
                        
                        if (rank == 0)
                          printf("[ERROR] maximum_neighbor_count needs to be increased\n");
                        return PIDX_err_rst;
                      }
                      
                      patch_grp->count = patch_count;
                    }
                    free(curr_patch);
                  }
                }
                
                patch_grp->max_patch_rank = patch_grp->source_patch_rank[0];
                max_vol = 1;
                for(d = 0; d < PIDX_MAX_DIMENSIONS; d++)
                  max_vol = max_vol * patch_grp->patch[0]->size[d];
                int64_t c_vol = 1;
                for(c = 1; c < patch_grp->count ; c++)
                {
                  c_vol = 1;
                  for(d = 0; d < PIDX_MAX_DIMENSIONS; d++)
                    c_vol = c_vol * patch_grp->patch[c]->size[d];
                  if(c_vol > max_vol)
                  {
                    max_vol = c_vol;
                    patch_grp->max_patch_rank = patch_grp->source_patch_rank[c];
                  }
                }
                
                printf("max_patch_rank %d\n", patch_grp->max_patch_rank);
                if(rank == patch_grp->max_patch_rank)
                  var0->patch_group_count = var0->patch_group_count + 1;
                printf("%d\n", var0->patch_group_count);
                reg_patch_count++;
              }
              free(reg_patch);
            }
    
    free(local_proc_patch);
    //free(rank_r_offset);
    //free(rank_r_count);
    
    //return num_output_buffers;
  }
#else
  rst_id->idx->enable_rst = 0;
  var0->patch_group_count = var0->sim_patch_count;
#endif
  
  for (v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    PIDX_variable var = rst_id->idx->variable[v];
    var->patch_group_count = var0->patch_group_count;
    
    var->patch_group_count = rst_id->idx->variable[rst_id->first_index]->patch_group_count;
    
    var->rst_patch_group = malloc(var->patch_group_count * sizeof(*(var->rst_patch_group)));
    memset(var->rst_patch_group, 0, var->patch_group_count * sizeof(*(var->rst_patch_group)));
    for (p = 0; p < var->patch_group_count; p++)
    {
      var->rst_patch_group[p] = malloc(sizeof(*(var->rst_patch_group[p])));
      memset(var->rst_patch_group[p], 0, sizeof(*(var->rst_patch_group[p])));
    }
  }
  
  j = 0;
  v = 0;
  p = 0;
  if(rst_id->idx->enable_rst == 1)
  {
#if PIDX_HAVE_MPI
    int rank = 0, cnt = 0, i = 0;
    MPI_Comm_rank(rst_id->comm, &rank);
    for (v = rst_id->first_index; v <= rst_id->last_index; v++)
    {
      PIDX_variable var = rst_id->idx->variable[v];
      cnt = 0;
      for (i = 0; i < rst_id->reg_patch_grp_count; i++)
      {
        if (rank == rst_id->reg_patch_grp[i]->max_patch_rank)
        {
          Ndim_patch_group patch_group = var->rst_patch_group[cnt];
          patch_group->count = rst_id->reg_patch_grp[i]->count;
          patch_group->type = rst_id->reg_patch_grp[i]->type;
          patch_group->patch = malloc(sizeof(*(patch_group->patch)) * rst_id->reg_patch_grp[i]->count);
          memset(patch_group->patch, 0, sizeof(*(patch_group->patch)) * rst_id->reg_patch_grp[i]->count);
          
          patch_group->reg_patch = malloc(sizeof(*(patch_group->reg_patch)));
          memset(patch_group->reg_patch, 0, sizeof(*(patch_group->reg_patch)));
          
          for(j = 0; j < rst_id->reg_patch_grp[i]->count; j++)
          {
            patch_group->patch[j] = malloc(sizeof(*(patch_group->patch[j])));
            memset(patch_group->patch[j], 0, sizeof(*(patch_group->patch[j])));
            
            memcpy(patch_group->patch[j]->offset, rst_id->reg_patch_grp[i]->patch[j]->offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
            memcpy(patch_group->patch[j]->size, rst_id->reg_patch_grp[i]->patch[j]->size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
          }
          memcpy(patch_group->reg_patch->offset, rst_id->reg_patch_grp[i]->reg_patch->offset, sizeof(int64_t) * PIDX_MAX_DIMENSIONS);
          memcpy(patch_group->reg_patch->size, rst_id->reg_patch_grp[i]->reg_patch->size, sizeof(int64_t) * PIDX_MAX_DIMENSIONS);
          cnt++;
        }
      }
      if (cnt != var->patch_group_count)
        return PIDX_err_rst;
    }
#endif
  }
  else
  {
    for (v = rst_id->first_index; v <= rst_id->last_index; v++)
    {
      PIDX_variable var = rst_id->idx->variable[v];
      for (p = 0; p < var->patch_group_count; p++)
      {
        Ndim_patch_group patch_group = var->rst_patch_group[p];
        patch_group->count = 1;
        patch_group->type = 0;
        patch_group->patch = malloc(sizeof(*(patch_group->patch)) * patch_group->count);
        memset(patch_group->patch, 0, sizeof(*(patch_group->patch)) * patch_group->count);
        
        patch_group->reg_patch = malloc(sizeof(*(patch_group->reg_patch)));
        memset(patch_group->reg_patch, 0, sizeof(*(patch_group->reg_patch)));
        
        for(j = 0; j < patch_group->count; j++)
        {
          patch_group->patch[j] = malloc(sizeof(*(patch_group->patch[j])));
          memset(patch_group->patch[j], 0, sizeof(*(patch_group->patch[j])));
          
          memcpy(patch_group->patch[j]->offset, var->sim_patch[p]->offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
          memcpy(patch_group->patch[j]->size, var->sim_patch[p]->size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
        }
        memcpy(patch_group->reg_patch->offset, var->sim_patch[p]->offset, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
        memcpy(patch_group->reg_patch->size, var->sim_patch[p]->size, PIDX_MAX_DIMENSIONS * sizeof(int64_t));
      }
    }
  }
  
  return PIDX_success;
}


PIDX_return_code PIDX_multi_patch_rst_buf_create(PIDX_multi_patch_rst_id multi_patch_rst_id)
{
  return PIDX_err_not_implemented;
}


PIDX_return_code PIDX_multi_patch_rst_write(PIDX_multi_patch_rst_id multi_patch_rst_id)
{
  return PIDX_err_not_implemented;
}


PIDX_return_code PIDX_multi_patch_rst_read(PIDX_multi_patch_rst_id multi_patch_rst_id)
{
  return PIDX_err_not_implemented;
}



PIDX_return_code PIDX_multi_patch_rst_buf_destroy(PIDX_multi_patch_rst_id multi_patch_rst_id)
{
  return PIDX_err_not_implemented;
}


PIDX_return_code PIDX_multi_patch_rst_buf_aggregate_read(PIDX_multi_patch_rst_id multi_patch_rst_id)
{
  return PIDX_err_not_implemented;
}

PIDX_return_code PIDX_multi_patch_rst_aggregate_buf_destroy(PIDX_multi_patch_rst_id multi_patch_rst_id)
{
  return PIDX_err_not_implemented;
}



PIDX_return_code PIDX_multi_patch_rst_buf_aggregate_write(PIDX_multi_patch_rst_id multi_patch_rst_id)
{
  return PIDX_err_not_implemented;
}


PIDX_return_code PIDX_multi_patch_rst_meta_data_destroy(PIDX_multi_patch_rst_id multi_patch_rst_id)
{
  return PIDX_err_not_implemented;
}

PIDX_return_code PIDX_multi_patch_rst_finalize(PIDX_multi_patch_rst_id multi_patch_rst_id)
{
  free(multi_patch_rst_id);
  multi_patch_rst_id = 0;

  return PIDX_success;
}
