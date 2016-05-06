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

#ifndef _pidx_metadata_h
#define _pidx_metadata_h

#include <stdint.h>
#include "PIDX_error_codes.h"

class TiXmlDocument;
  
struct PIDX_metadata_struct{
  char filename[1024];
  TiXmlDocument* doc;
};
typedef struct PIDX_metadata_struct* PIDX_metadata;
  
PIDX_return_code PIDX_metadata_load(PIDX_metadata* metadata, const char* _filename);

PIDX_return_code PIDX_metadata_create(PIDX_metadata* metadata, const char* _filename);

PIDX_return_code PIDX_metadata_save(PIDX_metadata metadata);

PIDX_return_code PIDX_metadata_add_timestep(PIDX_metadata metadata, int index, double value);

PIDX_return_code PIDX_metadata_get_timestep(PIDX_metadata metadata, int index, double& value);

PIDX_return_code PIDX_metadata_add_simple_box(PIDX_metadata metadata, int64_t* log_size, double* phy_size);

PIDX_return_code PIDX_metadata_destroy(PIDX_metadata metadata);

#endif
