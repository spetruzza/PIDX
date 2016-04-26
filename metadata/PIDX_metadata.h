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

class TiXmlDocument;
  
struct PIDX_metadata_struct{
  char filename[1024];
  TiXmlDocument* doc;
};
typedef struct PIDX_metadata_struct* PIDX_metadata;
  
int PIDX_metadata_load(PIDX_metadata* metadata, const char* _filename);

int PIDX_metadata_create(PIDX_metadata* metadata, const char* _filename);

int PIDX_metadata_save(PIDX_metadata metadata);

int PIDX_metadata_add_timestep(PIDX_metadata metadata, int index, float value);

int PIDX_metadata_get_timestep(PIDX_metadata metadata, int index, float& value);

int PIDX_metadata_destroy(PIDX_metadata metadata);

#endif
