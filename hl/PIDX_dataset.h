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

#if PIDX_HAVE_METADATA
#include "metadata/PIDX_metadata.h"
#endif

class PIDX_Dataset
{
public:

 PIDX_Dataset(int* global_size_ptr, double* phy_dim_ptr=NULL, 
 				int* local_offset_ptr=NULL, int* local_size_ptr=NULL);

 void open(std::string name, PIDX_flags flags);
 void write(std::string var_name, const void* buf, PIDX_data_type dtype);
 void read(std::string var_name, void* buf);
 void read(int var_index, void* buf);
 void setCurrentTime(int time_index, double simtime);
 int getTimeIndex(double simtime);
 void close();

private:

 int process_count, rank;
 std::string filename;

 PIDX_file file;            
 PIDX_access access;
 PIDX_point global_size, local_offset, local_size;
 PIDX_variable* variable;   // variable descriptor

 #if PIDX_HAVE_METADATA
 PIDX_metadata metadata;
 #endif

 int first_tstep, last_tstep;
 int variable_count;
 int curr_tstep;
 double curr_simtime;
 double phy_dim[3];
 bool write_mode;
 
 std::map<std::string, int> var_map;

};
  
#endif
