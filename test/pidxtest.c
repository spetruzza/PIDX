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

#include "pidxtest.h"
#include "testdefs.h"

#include <getopt.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>


/*usage*/
static void usage(enum Kind kind) 
{
  printf("Usage for pidxtest -k %s:\n\n",kindToStr(kind));

  switch (kind)
  {
  //case READER:          usage_reader();    break;
  case WRITER:          usage_writer();    break;
  case ONE_VAR_WRITER:  usage_one_var_writer();    break;
  case MULTI_VAR_WRITER:  usage_multi_var_writer();    break;
  case DEFAULT:
  default:              usage_one_var_writer();
  }
}

/*main*/
int main(int argc, char **argv) 
{
  int ret, nprocs, rank;   /* process count and rank */
  struct Args args;        /* initialize args */

  /*MPI initialization*/
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  /* Rank 0 parses the command line arguments */
  if (rank == 0) 
  {
    ret = parse_args(&args, argc, argv);
    if (ret < 0) 
    {
      usage(args.kind);
      print_error("ret error", __FILE__, __LINE__);
    }
    if (args.count_local[0] == 0 || args.count_local[1] == 0 || args.count_local[2] == 0) 
    {
      usage(args.kind);
      print_error("Local Dimension cannot be 0", __FILE__, __LINE__);
    } 
    else 
    {
      if ((args.extents[0] / args.count_local[0]) * (args.extents[1] / args.count_local[1]) * (args.extents[2] / args.count_local[2]) != nprocs) 
      {
	usage(args.kind);
	print_error("Wrong Number of Processes\n", __FILE__, __LINE__);
      }
    }
  }
  
  MPI_Bcast(&args.kind, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  /* run the specified test */
  
  switch (args.kind)
  {
    //case READER:          return test_reader(args, rank, nprocs);
    case WRITER:          
      test_writer(args, rank, nprocs);
      break;
    case ONE_VAR_WRITER:  
      test_one_var_writer(args, rank, nprocs);
      break;
    case MULTI_VAR_WRITER:  
      test_multi_var_writer(args, rank, nprocs);
      break;
    default:
      test_one_var_writer(args, rank, nprocs);
  }
  
  MPI_Finalize();
  return 0;
}

/*parse_args*/
int parse_args(struct Args *args, int argc, char **argv) 
{
  char flags[] = "k:g:l:f:t:v:";
  int one_opt = 0;
  char strkind[128];

  if (!args)
    return -1;
  
  memset(args, 0, sizeof(struct Args)); 
  
  while ((one_opt = getopt(argc, argv, flags)) != EOF) 
  {
    /* postpone error checking for after while loop */
    switch (one_opt) 
    {
      case('k'):
	  sprintf(strkind, "%s", optarg);
	  args->kind = strToKind(strkind);
	  break;
      case('g'):
	  sscanf(optarg, "%dx%dx%d", &args->extents[0], &args->extents[1], &args->extents[2]);
	  break;
      case('l'):
	  sscanf(optarg, "%dx%dx%d", &args->count_local[0], &args->count_local[1], &args->count_local[2]);
	  break;
      case('f'):
	  sprintf(args->output_file_template, "%s", optarg);
	  break;
      case('t'):
	  sscanf(optarg, "%d", &args->time_step);
	  break;
      case('v'):
	  sscanf(optarg, "%d", &args->variable_count);
	  break;
      case('?'):
	  return (-1);
    }
  }
  
  printf("COUNT = %d\n", args->variable_count);
  
  /* need positive dimensions */
  if (args->extents[0] < 1 || args->extents[1] < 1 || args->extents[2] < 1 || args->count_local[0] < 1 || args->count_local[1] < 1 || args->count_local[2] < 1) 
  {
    printf("Error: bad dimension specification.\n");
    return (-1);
  }

  /* need global dimension to be larger than the local */
  if (args->extents[0] < args->count_local[0] || args->extents[1] < args->count_local[1] || args->extents[2] < args->count_local[2]) 
  {
    printf("Error: global dimensions and local dimensions aren't evenly divisible\n");
    return (-1);
  }

  args->extents[3] = 1;
  args->extents[4] = 1;

  return (0);
}

/*print_error*/
int print_error(char *error_message, char* file, int line) 
{
  fprintf(stderr, "File [%s] Line [%d] Error [%s]\n", error_message, line, file);
  MPI_Abort(MPI_COMM_WORLD, -1);
  return 0;
}

/*kindToStr*/
char* kindToStr(enum Kind k)
{
  switch (k)
  {
    //case READER:          return "reader";
    case WRITER:           return "writer";
    case ONE_VAR_WRITER:   return "one_var";
    case MULTI_VAR_WRITER:   return "multi_var";
    case DEFAULT:
    default:               return "default";
  }
}

/*strToKind*/
enum Kind strToKind(const char *str)
{
//   if (strcmp(str,"reader")    == 0) return READER;
  if (strcmp(str,"writer")    == 0) return WRITER;
  if (strcmp(str,"one_var")    == 0) return ONE_VAR_WRITER;
  if (strcmp(str,"multi_var")    == 0) return MULTI_VAR_WRITER;
  else                              return DEFAULT;
}