#include <cstdio>
#include "PIDX_metadata.h"



int main(){
  PIDX_metadata metadata;
  
  PIDX_metadata_create(&metadata, "data.xml");
  PIDX_metadata_add_timestep(metadata, 0, 0.000000000000);
  PIDX_metadata_add_timestep(metadata, 2, 0.000000001);
  PIDX_metadata_add_timestep(metadata, 3, 2.998);
  
  double findv;
  PIDX_metadata_get_timestep(metadata, 2, findv);
  printf("2 value %g\n", findv);
  PIDX_metadata_save(metadata);
  
  PIDX_metadata_load(&metadata, "data.xml");
  PIDX_metadata_add_timestep(metadata, 4, 3.598);
  PIDX_metadata_save(metadata);

  long int log_size[3] = {300, 450, 680};
  double phy_size[3] = {12.5, 45.330, 6.80};
  PIDX_metadata_add_simple_box(metadata, log_size, phy_size);
  PIDX_metadata_save(metadata);
  
  return 0;
}

/*
int main(){
  
  MetadataIDX& meta = *MetadataIDX::getInstance();
  
  meta.create("data.xml");
  meta.addTimeStep(1, 0.998);
  meta.addTimeStep(2, 1.998);
  meta.addTimeStep(3, 2.998);
  
  float findv;
  meta.getTimeStep(3, findv);
  printf("3 value %f\n", findv);
  meta.save();
  
  meta.load("data.xml");
  meta.addTimeStep(4, 3.598);
  meta.save();
  
  return 0;
}
*/
