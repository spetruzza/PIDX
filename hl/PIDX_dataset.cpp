/***************************************************
 ** ViSUS Visualization Project                    **
 ** Copyright (c) 2010 University of Utah          **
 ** Scientific Computing and Imaging Institute     **
 ** 72 S Central Campus Drive, Room 3750           **
 ** Salt Lake City, UT 84112                       **
 **                                                **
 ** For information about this project see:        **
 ** http://www.pascucci.org/visus/                 **
 **                                                **
 **      or contact: pascucci@sci.utah.edu         **
 **                                                **
 ****************************************************/

#include "PIDX_dataset.h"

//static MPI_Comm NEW_COMM_WORLD; 

void terminate()
{
#if PIDX_HAVE_MPI
  MPI_Abort(MPI_COMM_WORLD, -1);
#else
  exit(-1);
#endif
}

void terminate_with_error_msg(const char *format, ...)
{
  va_list arg_ptr;
  va_start(arg_ptr, format);
  vfprintf(stderr, format, arg_ptr);
  va_end(arg_ptr);
  terminate();
}

PIDX_Dataset::PIDX_Dataset(int* global_size_ptr, 
	double* phy_dim_ptr, int* local_offset_ptr, int* local_size_ptr)
{
//	if(MPI_Comm_dup(MPI_COMM_WORLD, &NEW_COMM_WORLD)!= MPI_SUCCESS)
//    terminate_with_error_msg("ERROR: MPI_Comm_dup error\n");

  if (MPI_Comm_size(MPI_COMM_WORLD, &process_count) != MPI_SUCCESS)
    terminate_with_error_msg("ERROR: MPI_Comm_size error\n");
  if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS)
    terminate_with_error_msg("ERROR: MPI_Comm_rank error\n");

  variable = NULL;

  if(global_size_ptr == NULL)
  {
  	fprintf(stderr, "NULL global size\n");
  	assert(false);
  	return;
  }

  if(phy_dim_ptr != NULL)
    memcpy(&phy_dim[0], phy_dim_ptr, 3*sizeof(double));
  else
  {
  	phy_dim[0] = (double) global_size_ptr[0];
  	phy_dim[1] = (double) global_size_ptr[1];
  	phy_dim[2] = (double) global_size_ptr[2];
  }
  if(local_offset_ptr == NULL)
  {
  	local_offset_ptr = (int*)malloc(sizeof(int)*3);
  	memset(local_offset_ptr, 0, sizeof(int)*3);
  }

  if(local_size_ptr == NULL)
  {
  	local_size_ptr = global_size_ptr;
  }

  PIDX_set_point_5D(global_size, (int64_t)global_size_ptr[0], (int64_t)global_size_ptr[1], (int64_t)global_size_ptr[2], 1, 1);
  PIDX_set_point_5D(local_offset, (int64_t)local_offset_ptr[0], (int64_t)local_offset_ptr[1],(int64_t)local_offset_ptr[2], 0, 0);
  PIDX_set_point_5D(local_size, (int64_t)local_size_ptr[0], (int64_t)local_size_ptr[1], (int64_t)local_size_ptr[2], 1, 1);

}


void PIDX_Dataset::open(std::string name, PIDX_flags flags)
{

  filename = name; 
	write_mode = false;
  variable_count = 0;
  
  int ret = PIDX_create_access(&access);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_create_access");
#if PIDX_HAVE_MPI
  ret = PIDX_set_mpi_access(access, MPI_COMM_WORLD);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_mpi_access");
#endif

  char output_file_name[512];
  sprintf(output_file_name, "%s.idx", name.c_str());

  if( flags == PIDX_MODE_CREATE || flags == PIDX_MODE_WRONLY )
  {  	
	  ret = PIDX_file_create(output_file_name, flags, access, &file);
	  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_file_create");
	  ret = PIDX_set_dims(file, global_size);
	  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_dims");
	  ret = PIDX_set_current_time_step(file, 0);
	  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_current_time_step"); 

 		//PIDX_set_block_count(file, 256);
  	//PIDX_set_block_size(file, 15);

	  write_mode = true; 
  }
  else
  {
  	if(rank == 0)
  		printf("opening %s\n", output_file_name);
	  //  PIDX mandatory calls
	  ret = PIDX_file_open(output_file_name, flags, access, &file);
	  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_file_create");

	  ret = PIDX_get_dims(file, global_size);
	  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_dims");

	  ret = PIDX_get_variable_count(file, &variable_count);
	  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_variable_count");

	  ret = PIDX_get_first_tstep(file, &first_tstep);
	  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_get_first_tstep");
	  ret = PIDX_get_last_tstep(file, &last_tstep);
	  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_get_last_tstep");

	  int ntsteps = last_tstep-first_tstep;
	  if(rank == 0)
	    printf("found %d variables ntsteps %d ft %d lt %d\n", variable_count, ntsteps, first_tstep, last_tstep);

	  variable = (PIDX_variable*)malloc(sizeof(*variable) * variable_count);
	  memset(variable, 0, sizeof(*variable) * variable_count);

	  ret = PIDX_set_current_time_step(file, 0);
	  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_current_time_step");
	 
	  int var = 0;
	  for (var = 0; var < variable_count; var++)
	  {
	    ret = PIDX_get_next_variable(file, &variable[var]);
	    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_get_next_variable");

	    std::string varname = variable[var]->var_name;
	    var_map[varname] = var;

	   // printf("adding var %s to index %d\n", varname.c_str(), var);
	    ret = PIDX_read_next_variable(file, variable[var]);
	    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_read_next_variable");
	  }

	  ret = PIDX_reset_variable_counter(file);
	  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_reset_variable_counter");
	}

#if PIDX_HAVE_METADATA
  if(rank==0)
  {
    std::string metadata_filename = "./"+filename+"/"+filename+".xml";
    struct stat temp;
    if( stat (metadata_filename.c_str(), &temp) == 0 ) {
      PIDX_metadata_load(&metadata, metadata_filename.c_str());
    }
    else
    {
      PIDX_metadata_create(&metadata, metadata_filename.c_str());
      PIDX_metadata_add_simple_box(metadata, global_size, phy_dim);
      PIDX_metadata_save(metadata);
    }
  }
#endif

}

void PIDX_Dataset::write(std::string var_name, const void* buf, PIDX_data_type dtype)
{
  PIDX_variable variable;
  int ret = PIDX_set_variable_count(file, ++variable_count);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_variable_count");

  char name[512];
  sprintf(name, "%s", var_name.c_str());
  
  ret = PIDX_set_current_time_step(file, curr_tstep);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_current_time_step");
   
  int bits_per_sample = 0;
  ret = PIDX_default_bits_per_datatype(dtype, &bits_per_sample);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_default_bytes_per_datatype");

  ret = PIDX_variable_create(name,  bits_per_sample, dtype, &variable);
  if (ret != PIDX_success)  terminate_with_error_msg("A PIDX_variable_create_pressure");

  ret = PIDX_variable_write_data_layout(variable, local_offset, local_size, buf, PIDX_row_major);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_data_layout");

  ret = PIDX_append_and_write_variable(file, variable);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_append_and_write_variable");

  write_mode = true;

}

void PIDX_Dataset::close()
{
	int ret = PIDX_close(file);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_close");

  ret = PIDX_close_access(access);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_close_access");

#if PIDX_HAVE_METADATA 
  if(rank==0 && write_mode) // save timestep information (if not present)
  { 
    double curr_t = 0;
    ret = PIDX_metadata_get_timestep(metadata, curr_tstep, curr_t);
    
    if(ret != PIDX_success)
    {
      PIDX_metadata_add_timestep(metadata, curr_tstep, curr_simtime);
      PIDX_metadata_save(metadata);
    }    
  }
#endif

  if(variable != NULL)
  { 
    free(variable);
    variable = NULL;
  }
}

void PIDX_Dataset::setCurrentTime(int time_index, double simtime){
  curr_tstep = time_index;
  curr_simtime = simtime;
}

int PIDX_Dataset::getTimeIndex(double simtime){

  int time_step_count = -1;

#if PIDX_HAVE_METADATA
  int ret = 0;

  if(rank == 0)
  {
    //printf("search in %d to %d\n", first_tstep, last_tstep);
    for(int ts=first_tstep; ts < last_tstep; ts++)
    {
      double curr_t = 0;
      ret = PIDX_metadata_get_timestep(metadata, ts, curr_t);
      
      if (ret != PIDX_success)  terminate_with_error_msg("PIDX_metadata_get_timestep");
          
      if(curr_t == simtime)
      {
        time_step_count = ts;
        break;
      }
    }
    if(time_step_count == -1)
    {
      fprintf(stderr,"Timestep not found!\n");
      assert(false);
      return -1;
    }
	}	

  MPI_Bcast(&time_step_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
#else
  fprintf(stderr, "No PIDX metadata module installed\n");
  assert(false);
#endif

  return time_step_count;
}

void PIDX_Dataset::read(std::string var_name, void* buf){
	std::map<std::string, int>::iterator varit = var_map.find(var_name);

	if(varit == var_map.end())
	{
		fprintf(stderr, "Variable %s not found\n", var_name.c_str());
	}
	else
		read(varit->second, buf);
}

void PIDX_Dataset::read(int var_index, void* buf)
{	 
  int ret = PIDX_set_current_time_step(file, curr_tstep);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_current_time_step");
 
  ret = PIDX_get_next_variable(file, &variable[var_index]);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_get_next_variable");

  ret = PIDX_variable_read_data_layout(variable[var_index], local_offset, local_size, buf, PIDX_row_major);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_read_data_layout");

  ret = PIDX_read_next_variable(file, variable[var_index]);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_read_next_variable");

  write_mode = false;
}
