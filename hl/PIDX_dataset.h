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

#ifndef _pidx_dataset_h
#define _pidx_dataset_h

#include <stdint.h>
#include <string>
#include "PIDX.h"
#include "metadata/PIDX_metadata.h"

class PIDX_Dataset{
public:
 PIDX_Dataset(int* global_size_ptr, double* phy_dim_ptr=NULL, 
 				int* local_offset_ptr=NULL, int* local_size_ptr=NULL);

 void open(std::string name);
 void write(std::string var_name, const void* buf, PIDX_data_type dtype, double simtime=0.0);
 void read(std::string var_name, void* buf, double simtime);
 void read(int var_index, void* buf, double simtime);
 void close();

private:

 int process_count, rank;
 std::string filename;


 PIDX_file file;            
 PIDX_access access;
 PIDX_point global_size, local_offset, local_size;
 PIDX_metadata metadata;
 PIDX_variable* variable;   // variable descriptor
 int *v_per_sample;
 int *bits_per_sample;
 char **type_name;
 int first_tstep, last_tstep;
 int variable_count;
 double phy_dim[3];
 int ntsteps;
 int curr_tstep;
 double last_simtime;

 std::map<std::string, int> var_map;

};
  
#endif
