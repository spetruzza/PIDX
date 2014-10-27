/*****************************************************
 **  PIDX Parallel I/O Library                      **
 **  Copyright (c) 2010-2014 University of Utah     **
 **  Scientific Computing and Imaging Institute     **
 **  72 S Central Campus Drive, Room 3750           **
 **  Salt Lake City, UT 84112                       ***
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
 /**
 * \file PIDX_data_layout.h
 *
 * \mainpage
 *
 * \author Sidharth Kumar
 * \author Cameron Christensen
 * \author Giorgio Scorzelli
 * \author Valerio Pascucci
 * \date   10/09/14
 *
 * PIDX is an I/O library that enables HPC applications to write distributed 
 * multi-dimensional data directly into a hierarchical multi-resolution 
 * data format (IDX) with minimal overhead. 
 *
 */
 
#ifndef __PIDX_DATA_LAYOUT_H
#define __PIDX_DATA_LAYOUT_H

typedef unsigned int PIDX_data_layout;

/// Data layout options 
/// \param PIDX_row_major row major layout
/// \param PIDX_column_major column major layout
///
extern PIDX_data_layout PIDX_row_major;
extern PIDX_data_layout PIDX_column_major;

#endif
