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

#include "../../PIDX_inc.h"
#define MAX_TEMPLATE_DEPTH 6

static uint32_t* headers;
static int populate_meta_data(PIDX_header_io_id header_io_id, PIDX_block_layout block_layout, int file_number, char* bin_file, int mode);

struct PIDX_header_io_struct 
{
#if PIDX_HAVE_MPI
  MPI_Comm comm;
#endif
    
  //Contains all relevant IDX file info
  //Blocks per file, samples per block, bitmask, box, file name template and more
  idx_dataset idx;
  
  //Contains all derieved IDX file info
  //number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_d;

  int enable_raw_dump;

  int first_index;
  int last_index;
  char filename_template[1024];
};

PIDX_header_io_id PIDX_header_io_init(idx_dataset idx_meta_data,
                                      idx_dataset_derived_metadata idx_d,
                                      int first_index, int last_index )
{
  PIDX_header_io_id header_io_id;

  //Creating the IO ID
  header_io_id = (PIDX_header_io_id)malloc(sizeof (*header_io_id));
  if (!header_io_id) 
    memset(header_io_id, 0, sizeof (*header_io_id));
  
  header_io_id->idx = idx_meta_data;
  header_io_id->idx_d = idx_d;

  header_io_id->first_index = first_index;
  header_io_id->last_index = last_index;

  header_io_id->enable_raw_dump = 0;

  if (first_index == 0)
  {
    int total_header_size;

    total_header_size = (10 + (10 * header_io_id->idx->blocks_per_file)) * sizeof (uint32_t) * header_io_id->idx->variable_count;
    header_io_id->idx_d->start_fs_block = total_header_size / header_io_id->idx_d->fs_block_size;
    if (total_header_size % header_io_id->idx_d->fs_block_size)
      header_io_id->idx_d->start_fs_block++;

    headers = (uint32_t*)malloc(total_header_size);
    memset(headers, 0, total_header_size);
  }

  return header_io_id;
}

#if PIDX_HAVE_MPI
PIDX_return_code PIDX_header_io_set_communicator(PIDX_header_io_id header_io, MPI_Comm comm)
{
  if (header_io == NULL)
    return PIDX_err_id;

  header_io->comm = comm;

  return PIDX_success;
}
#endif

PIDX_return_code PIDX_header_io_enable_raw_dump(PIDX_header_io_id header_io)
{
  if (header_io == NULL)
    return PIDX_err_id;

  header_io->enable_raw_dump = 1;

  return PIDX_success;
}


PIDX_return_code PIDX_header_io_write_idx (PIDX_header_io_id header_io, char* data_set_path, int current_time_step)
{
  int l = 0, rank = 0, N;
  FILE* idx_file_p;
  char dirname[1024], basename[1024];
  
#if PIDX_HAVE_MPI
  if (header_io->idx_d->parallel_mode == 1)
    MPI_Comm_rank(header_io->comm, &rank);
#endif
  
  int nbits_blocknumber = (header_io->idx_d->maxh - header_io->idx->bits_per_block - 1);
  VisusSplitFilename(data_set_path, dirname, basename);

  //remove suffix
  for (N = strlen(basename) - 1; N >= 0; N--) 
  {
    int ch = basename[N];
    basename[N] = 0;
    if (ch == '.') break;
  }

#if 0
  //if i put . as the first character, if I move files VisusOpen can do path remapping
  sprintf(pidx->filename_template, "./%s", basename);
#endif
  //pidx does not do path remapping 
  strcpy(header_io->filename_template, data_set_path);
  for (N = strlen(header_io->filename_template) - 1; N >= 0; N--) 
  {
    int ch = header_io->filename_template[N];
    header_io->filename_template[N] = 0;
    if (ch == '.') break;
  }

  //can happen if I have only only one block
  if (nbits_blocknumber == 0) 
    strcat(header_io->filename_template, "/%01x.bin");
   
  else 
  {
    //approximate to 4 bits
    if (nbits_blocknumber % 4) 
    {
      nbits_blocknumber += (4 - (nbits_blocknumber % 4));
      //assert(!(nbits_blocknumber % 4));
    }
    if (nbits_blocknumber <= 8) 
      strcat(header_io->filename_template, "/%02x.bin"); //no directories, 256 files
    else if (nbits_blocknumber <= 12) 
      strcat(header_io->filename_template, "/%03x.bin"); //no directories, 4096 files
    else if (nbits_blocknumber <= 16) 
      strcat(header_io->filename_template, "/%04x.bin"); //no directories, 65536  files
    else 
    {
      while (nbits_blocknumber > 16) 
      {
        strcat(header_io->filename_template, "/%02x"); //256 subdirectories
        nbits_blocknumber -= 8;
      }
      strcat(header_io->filename_template, "/%04x.bin"); //max 65536  files
      nbits_blocknumber -= 16;
      //assert(nbits_blocknumber <= 0);
    }
  }
  
  if (strncmp(".idx", &data_set_path[strlen(data_set_path) - 4], 4) != 0) 
  {
    fprintf(stderr, "[%s] [%d] Bad file name extension.\n", __FILE__, __LINE__);
    return 1;
  }
    
  if (rank == 0)
  {
    //fprintf(stderr, "writing IDX file...\n", __FILE__, __LINE__);

    if (header_io->idx->compression_type != PIDX_NO_COMPRESSION)
    {
      char visus_data_path[PATH_MAX];
      char filename_skeleton[PATH_MAX];
      strncpy(filename_skeleton, data_set_path, strlen(data_set_path) - 4);
      filename_skeleton[strlen(filename_skeleton) - 4] = '\0';
      sprintf(visus_data_path, "%s_visus.idx", filename_skeleton);

      idx_file_p = fopen(visus_data_path, "w");
      if (!idx_file_p)
      {
        fprintf(stderr, " [%s] [%d] idx_dir is corrupt.\n", __FILE__, __LINE__);
        return -1;
      }

      fprintf(idx_file_p, "(version)\n6\n");
      fprintf(idx_file_p, "(logic_to_physic)\n%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", header_io->idx->transform[0], header_io->idx->transform[1], header_io->idx->transform[2], header_io->idx->transform[3], header_io->idx->transform[4], header_io->idx->transform[5], header_io->idx->transform[6], header_io->idx->transform[7], header_io->idx->transform[8], header_io->idx->transform[9], header_io->idx->transform[10], header_io->idx->transform[11], header_io->idx->transform[12], header_io->idx->transform[13], header_io->idx->transform[14], header_io->idx->transform[15]);


      if (header_io->idx->compression_type != PIDX_NO_COMPRESSION)
      {
        fprintf(idx_file_p, "(compression type)\n%d\n", header_io->idx->compression_type);
        fprintf(idx_file_p, "(box)\n0 %lld 0 %lld 0 %lld 0 %lld 0 %lld\n", (long long)(header_io->idx->chunked_bounds[0] - 1), (long long)(header_io->idx->chunked_bounds[1] - 1), (long long)(header_io->idx->chunked_bounds[2] - 1), (long long)(header_io->idx->chunked_bounds[3] - 1), (long long)(header_io->idx->chunked_bounds[4] - 1));

        fprintf(idx_file_p, "(original box)\n0 %lld 0 %lld 0 %lld 0 %lld 0 %lld\n", (long long)(header_io->idx->bounds[0] - 1), (long long)(header_io->idx->bounds[1] - 1), (long long)(header_io->idx->bounds[2] - 1), (long long)(header_io->idx->bounds[3] - 1), (long long)(header_io->idx->bounds[4] - 1));

        fprintf(idx_file_p, "(compression bit rate)\n%d\n", header_io->idx->compression_bit_rate);

        fprintf(idx_file_p, "(fields)\n");
        for (l = 0; l < header_io->last_index; l++)
        {
          fprintf(idx_file_p, "%s %d*float64", header_io->idx->variable[l]->var_name, header_io->idx->compression_bit_rate);
          if (l != header_io->last_index - 1)
            fprintf(idx_file_p, " + \n");
        }
      }

      fprintf(idx_file_p, "\n(bits)\n%s\n", header_io->idx->bitSequence);
      fprintf(idx_file_p, "(bitsperblock)\n%d\n(blocksperfile)\n%d\n", header_io->idx->bits_per_block, header_io->idx->blocks_per_file);
      fprintf(idx_file_p, "(filename_template)\n./%s\n", header_io->filename_template);

      fprintf(idx_file_p, "(time)\n0 %d time%%09d/"/*note: uintah starts at timestep 1, but we shouldn't assume...*/, header_io->idx->current_time_step);
      fclose(idx_file_p);
    }

    idx_file_p = fopen(data_set_path, "w");
    if (!idx_file_p) 
    {
      fprintf(stderr, " [%s] [%d] idx_dir is corrupt.\n", __FILE__, __LINE__);
      return -1;
    }

    fprintf(idx_file_p, "(version)\n6\n");
    fprintf(idx_file_p, "(logic_to_physic)\n%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", header_io->idx->transform[0], header_io->idx->transform[1], header_io->idx->transform[2], header_io->idx->transform[3], header_io->idx->transform[4], header_io->idx->transform[5], header_io->idx->transform[6], header_io->idx->transform[7], header_io->idx->transform[8], header_io->idx->transform[9], header_io->idx->transform[10], header_io->idx->transform[11], header_io->idx->transform[12], header_io->idx->transform[13], header_io->idx->transform[14], header_io->idx->transform[15]);

    fprintf(idx_file_p, "(box)\n0 %lld 0 %lld 0 %lld 0 %lld 0 %lld\n", (long long)(header_io->idx->bounds[0] - 1), (long long)(header_io->idx->bounds[1] - 1), (long long)(header_io->idx->bounds[2] - 1), (long long)(header_io->idx->bounds[3] - 1), (long long)(header_io->idx->bounds[4] - 1));
    
    if (header_io->enable_raw_dump == 1)
      fprintf(idx_file_p, "(raw_dump)\n%lld %lld %lld\n", (long long)header_io->idx->reg_patch_size[0], (long long)header_io->idx->reg_patch_size[1], (long long)header_io->idx->reg_patch_size[2]);

    if (header_io->idx->compression_type != PIDX_NO_COMPRESSION)
    {
      fprintf(idx_file_p, "(compression type)\n%d\n", header_io->idx->compression_type);
      fprintf(idx_file_p, "(compressed box)\n%lld %lld %lld %lld %lld\n", (long long)(header_io->idx->chunk_size[0]), (long long)(header_io->idx->chunk_size[1]), (long long)(header_io->idx->chunk_size[2]), (long long)(header_io->idx->chunk_size[3]), (long long)(header_io->idx->chunk_size[4]));
      fprintf(idx_file_p, "(compression bit rate)\n%d\n", header_io->idx->compression_bit_rate);
    }
    fprintf(idx_file_p, "(fields)\n");

    for (l = 0; l < header_io->last_index; l++)
    {
      fprintf(idx_file_p, "%s %s", header_io->idx->variable[l]->var_name, header_io->idx->variable[l]->type_name);
      if (l != header_io->last_index - 1)
        fprintf(idx_file_p, " + \n");
    }

    if (header_io->enable_raw_dump != 1)
    {
      fprintf(idx_file_p, "\n(bits)\n%s\n", header_io->idx->bitSequence);
      fprintf(idx_file_p, "(bitsperblock)\n%d\n(blocksperfile)\n%d\n", header_io->idx->bits_per_block, header_io->idx->blocks_per_file);

      //if (header_io->idx_d->res_to != 0)
      //  fprintf(idx_file_p, "(resolution)\n%d\n", header_io->idx_d->res_to);

      fprintf(idx_file_p, "(filename_template)\n./%s\n", header_io->filename_template);
    }

    fprintf(idx_file_p, "\n(time)\n0 %d time%%09d/"/*note: uintah starts at timestep 1, but we shouldn't assume...*/, header_io->idx->current_time_step);
    fclose(idx_file_p);
  }
  
  return 0;
}

///
int PIDX_header_io_file_create(PIDX_header_io_id header_io_id, PIDX_block_layout block_layout)
{
  int i = 0, rank = 0, j, ret;
  char bin_file[PATH_MAX];
  char last_path[PATH_MAX] = {0};
  char this_path[PATH_MAX] = {0};
  char tmp_path[PATH_MAX] = {0};
  char* pos;

#if PIDX_HAVE_MPI
  int nprocs = 1;
  if (header_io_id->idx_d->parallel_mode == 1)
  {
    MPI_Comm_rank(header_io_id->comm, &rank);
    MPI_Comm_size(header_io_id->comm, &nprocs);
  }
#endif


  for (i = 0; i < header_io_id->idx_d->max_file_count; i++)
  {
#if PIDX_HAVE_MPI
    if (i % nprocs == rank && block_layout->file_bitmap[i] == 1)
#else
      if (rank == 0 && block_layout->file_bitmap[i] == 1)
#endif
      {
        ret = generate_file_name(header_io_id->idx->blocks_per_file, header_io_id->idx->filename_template, i, bin_file, PATH_MAX);
        if (ret == 1)
        {
          fprintf(stderr, "[%s] [%d] generate_file_name() failed.\n", __FILE__, __LINE__);
          return 1;
        }

        //TODO: the logic for creating the subdirectory hierarchy could
        //be made to be more efficient than this. This implementation
        //walks up the tree attempting to mkdir() each segment every
        //time we switch to a new directory when creating binary files.

        // see if we need to make parent directory
        strcpy(this_path, bin_file);
        if ((pos = rindex(this_path, '/')))
        {
          pos[1] = '\0';
          if (!strcmp(this_path, last_path) == 0)
          {
            //this file is in a previous directory than the last
            //one; we need to make sure that it exists and create
            //it if not.
            strcpy(last_path, this_path);
            memset(tmp_path, 0, PATH_MAX * sizeof (char));
            //walk up path and mkdir each segment
            for (j = 0; j < (int)strlen(this_path); j++)
            {
              if (j > 0 && this_path[j] == '/')
              {
                //printf("path = %s %s [T %s]\n", tmp_path, bin_file, header_io_id->idx->filename_template);
                ret = mkdir(tmp_path, S_IRWXU | S_IRWXG | S_IRWXO);
                if (ret != 0 && errno != EEXIST)
                {
                  perror("mkdir");
                  fprintf(stderr, "Error: failed to mkdir %s\n", tmp_path);
                  return 1;
                }
              }
              tmp_path[j] = this_path[j];
            }
          }
        }

#if PIDX_HAVE_MPI
        if (header_io_id->idx_d->parallel_mode == 1)
        {
          MPI_File fh = 0;
          MPI_File_open(MPI_COMM_SELF, bin_file, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
          MPI_File_close(&fh);
        }
        else
        {
          int fp = 0;
          fp = open(bin_file, O_CREAT, 0664);
          close(fp);
        }
#else
        int fp = 0;
        fp = open(bin_file, O_CREAT, 0664);
        close(fp);
#endif
      }
  }
  
#if PIDX_HAVE_MPI
  if (header_io_id->idx_d->parallel_mode == 1)
    MPI_Barrier(header_io_id->comm);
#endif
  
  return PIDX_success;
}

PIDX_return_code PIDX_header_io_file_write(PIDX_header_io_id header_io_id, PIDX_block_layout block_layout, int mode)
{
  int i = 0, rank = 0, ret;
  char bin_file[PATH_MAX];
  
#if PIDX_HAVE_MPI
  int nprocs = 1;
  if (header_io_id->idx_d->parallel_mode == 1)
  {
    MPI_Comm_rank(header_io_id->comm, &rank);
    MPI_Comm_size(header_io_id->comm, &nprocs);
  }
#endif

  for (i = 0; i < header_io_id->idx_d->max_file_count; i++)
  {
#if PIDX_HAVE_MPI
    if (i % nprocs == rank && block_layout->file_bitmap[i] == 1)
#else
    if (rank == 0 && block_layout->file_bitmap[i] == 1)
#endif
    {
      ret = generate_file_name(header_io_id->idx->blocks_per_file, header_io_id->idx->filename_template, i, bin_file, PATH_MAX);
      if (ret == 1)
      {
        fprintf(stderr, "[%s] [%d] generate_file_name() failed.\n", __FILE__, __LINE__);
        return 1;
      }

      ret = populate_meta_data(header_io_id, block_layout, i, bin_file, mode);
      if (ret != PIDX_success) return PIDX_err_header;
    }
  }
  
  return PIDX_success;
}

static int populate_meta_data(PIDX_header_io_id header_io_id, PIDX_block_layout block_layout, int file_number, char* bin_file, int mode)
{
  int block_negative_offset = 0;
  int i = 0, j = 0, k = 0, all_scalars = 0;
  off_t data_offset = 0, base_offset = 0;
  //int total_file_size = 0;

  int64_t total_chunk_size = (header_io_id->idx->chunk_size[0] * header_io_id->idx->chunk_size[1] * header_io_id->idx->chunk_size[2] * header_io_id->idx->chunk_size[3] * header_io_id->idx->chunk_size[4]);// / (64 / header_io_id->idx->compression_bit_rate);


  /*
  for (i = 0; i < header_io_id->idx->variable_count; i++)
  {
    if (header_io_id->idx->variable[i]->values_per_sample != 1)
    {
      all_scalars = 0;
      break;
    }
  }
  */

  int total_header_size = (10 + (10 * header_io_id->idx->blocks_per_file)) * sizeof (uint32_t) * header_io_id->idx->variable_count;
  //memset(headers, 0, total_header_size);

  //int nprocs;
  //MPI_Comm_size(header_io_id->comm, &nprocs);

  for (i = 0; i < header_io_id->idx->blocks_per_file; i++)
  {
    //if (nprocs == 1)
    //  printf("[B] [%d %d] : %d\n", header_io_id->idx->variable[header_io_id->first_index]->global_block_layout->resolution_from, header_io_id->idx->variable[header_io_id->first_index]->global_block_layout->resolution_to, i);
    if (PIDX_blocks_is_block_present((i + (header_io_id->idx->blocks_per_file * file_number)), header_io_id->idx->variable[header_io_id->first_index]->global_block_layout))
    {
        //if (nprocs == 1)
        //  printf("[A] %d\n", i);
      block_negative_offset = PIDX_blocks_find_negative_offset(header_io_id->idx->blocks_per_file, (i + (header_io_id->idx->blocks_per_file * file_number)), header_io_id->idx->variable[header_io_id->first_index]->global_block_layout);


      for (j = header_io_id->first_index; j < header_io_id->last_index; j++)
      //for (j = 0; j < header_io_id->idx->variable_count; j++)
      {
        base_offset = 0;
        if (all_scalars == 0)
        {
          for (k = 0; k < j; k++)
            base_offset = base_offset + ((block_layout->block_count_per_file[file_number]) * (header_io_id->idx->variable[k]->bits_per_value / 8) * total_chunk_size * header_io_id->idx_d->samples_per_block * header_io_id->idx->variable[k]->values_per_sample) / (header_io_id->idx->compression_factor);
            //base_offset = base_offset + ((header_io_id->idx->variable[header_io_id->first_index]->block_count_per_file[file_number]) * (header_io_id->idx->variable[k]->bits_per_value / 8) * total_chunk_size * header_io_id->idx_d->samples_per_block * header_io_id->idx->variable[k]->values_per_sample) / (header_io_id->idx->variable[k]->bits_per_value / header_io_id->idx->compression_bit_rate);
        }
        else
          base_offset =  j * (block_layout->block_count_per_file[file_number]) * (header_io_id->idx->variable[header_io_id->first_index]->bits_per_value / 8) * total_chunk_size * header_io_id->idx_d->samples_per_block * header_io_id->idx->variable[header_io_id->first_index]->values_per_sample / (header_io_id->idx->compression_factor);
          //base_offset =  j * (header_io_id->idx->variable[header_io_id->first_index]->block_count_per_file[file_number]) * (header_io_id->idx->variable[header_io_id->first_index]->bits_per_value / 8) * total_chunk_size * header_io_id->idx_d->samples_per_block * header_io_id->idx->variable[header_io_id->first_index]->values_per_sample / (header_io_id->idx->variable[j]->bits_per_value / header_io_id->idx->compression_bit_rate);

        data_offset = (((i) - block_negative_offset) * header_io_id->idx_d->samples_per_block) * (header_io_id->idx->variable[j]->bits_per_value / 8) * total_chunk_size * header_io_id->idx->variable[j]->values_per_sample  / (header_io_id->idx->compression_factor);

        //printf("BLOCK %d = %d + %d + %d (%d %d)\n", i, (int)base_offset, (int)data_offset, (int)(header_io_id->idx_d->start_fs_block * header_io_id->idx_d->fs_block_size), header_io_id->idx_d->start_fs_block, header_io_id->idx_d->fs_block_size);
        data_offset = base_offset + data_offset + header_io_id->idx_d->start_fs_block * header_io_id->idx_d->fs_block_size;

        //TODO
        //printf("%d %d: %d (%d * %d * %d * %d)\n", i, j, header_io_id->idx_d->samples_per_block * (header_io_id->idx->variable[j]->bits_per_value / 8) * total_chunk_size * header_io_id->idx->variable[j]->values_per_sample, header_io_id->idx_d->samples_per_block, (header_io_id->idx->variable[j]->bits_per_value / 8), total_chunk_size, header_io_id->idx->variable[j]->values_per_sample);
        //if (file_number == 2)
        //if (nprocs == 1)
        //printf("[%d] offset : count = %lld %lld\n", i, (unsigned long long)data_offset, (unsigned long long)(header_io_id->idx_d->samples_per_block * (header_io_id->idx->variable[j]->bits_per_value / 8) * total_chunk_size * header_io_id->idx->variable[j]->values_per_sample));
        headers[12 + ((i + (header_io_id->idx->blocks_per_file * j))*10 )] = htonl(data_offset);
        headers[14 + ((i + (header_io_id->idx->blocks_per_file * j))*10)] = htonl(header_io_id->idx_d->samples_per_block * (header_io_id->idx->variable[j]->bits_per_value / 8) * total_chunk_size * header_io_id->idx->variable[j]->values_per_sample / (header_io_id->idx->compression_factor));

        //total_file_size = data_offset + header_io_id->idx_d->samples_per_block * (header_io_id->idx->variable[j]->bits_per_value / 8) * total_chunk_size * header_io_id->idx->variable[j]->values_per_sample / (header_io_id->idx->compression_factor);
      }
    }
  }
#if 1
  if (mode == 1)
  {

#if PIDX_HAVE_MPI
    if (header_io_id->idx_d->parallel_mode == 1)
    {
      MPI_File fh;
      MPI_Status status;
      int ret = 0;
      ret = MPI_File_open(MPI_COMM_SELF, bin_file, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
      if (ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed on %s\n", __FILE__, __LINE__, bin_file);
        return PIDX_err_io;
      }

      /*
      if(total_file_size != 0)
      {
        ret = MPI_File_set_size(fh, total_file_size);
        if(ret != MPI_SUCCESS)
        {
          fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed.\n", __FILE__, __LINE__);
          return PIDX_err_io;
        }
      }
      */

      ret = MPI_File_write_at(fh, 0, headers, total_header_size, MPI_BYTE, &status);
      if (ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_write_at() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }

      ret = MPI_File_close(&fh);
      if (ret != MPI_SUCCESS)
      {
        fprintf(stderr, "[%s] [%d] MPI_File_open() failed on %s\n", __FILE__, __LINE__, bin_file);
        return PIDX_err_io;
      }
    }
    else
    {
      int fh = 0;
      fh = open(bin_file, O_WRONLY, 0664);
      //ftruncate(fh, total_file_size);
      ssize_t ret_size = pwrite(fh, headers, total_header_size, 0);
      if (ret_size != total_header_size)
      {
        fprintf(stderr, "[%s] [%d] pwrite() failed on %s\n", __FILE__, __LINE__, bin_file);
        return PIDX_err_io;
      }

      close(fh);
    }
#else
    int fh = 0;
    fh = open(bin_file, O_WRONLY, 0664);
    //ftruncate(fh, total_file_size);
    ssize_t ret_size = pwrite(fh, headers, total_header_size, 0);
    if (ret_size != total_header_size)
    {
      fprintf(stderr, "[%s] [%d] pwrite() failed on %s\n", __FILE__, __LINE__, bin_file);
      return PIDX_err_io;
    }
    close(fh);
#endif

    //int total_header_size = (10 + (10 * header_io_id->idx->blocks_per_file)) * sizeof (uint32_t) * header_io_id->idx->variable_count;
  }
#endif

  return PIDX_success;
}

PIDX_return_code PIDX_header_io_finalize(PIDX_header_io_id header_io_id)
{

  if (header_io_id->idx->variable_count == header_io_id->last_index)
    free(headers);

  free(header_io_id);
  header_io_id = 0;

  return PIDX_success;
}
