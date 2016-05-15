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

int main(int argc,  char** argv){

  MPI_Init(&argc,&argv);

  int global_size_ptr[3] = {64,64,64};

  PIDX_Dataset* dataset = new PIDX_Dataset(global_size_ptr);
  double* Pa = (double*)malloc(sizeof(double)*64*64*64);
  
  for(int i=0; i<64*64*64; i++)
    Pa[i] = (double)i;

  int ts = 0;
  double simtime = 0.0;

  // write
  dataset->open("test", PIDX_MODE_CREATE);
  dataset->setCurrentTime(ts, simtime);
  dataset->write("Pressure", Pa, FLOAT64);
  
  dataset->close();

  MPI_Abort(MPI_COMM_WORLD, -1);

}
