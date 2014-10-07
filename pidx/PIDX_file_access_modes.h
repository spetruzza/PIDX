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
#ifndef __PIDX_FILE_ACCESS_MODES_H
#define __PIDX_FILE_ACCESS_MODES_H

extern const unsigned PIDX_file_excl;             // Error creating a file that already exists.
extern const unsigned PIDX_file_trunc;           // Create the file if it does not exist.
extern const unsigned PIDX_file_rdwr;           // Read only.
extern const unsigned PIDX_file_rdonly;             // Reading and writing.
extern const unsigned PIDX_file_wronly;           // Write only. 

#endif