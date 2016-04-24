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
 
#ifndef __PIDX_UTILS_H
#define __PIDX_UTILS_H


#define Min2ab(a,b)      (((a)<=(b))?(a):(b))
#define Max2ab(a,b)      (((a)> (b))?(a):(b))
#define OffsetFor(_D_,_From_,_Off_) for ((_D_)=(_From_);(_D_)<(PIDX_MAX_DIMENSIONS+(_Off_));(_D_)++)
#define For(_D_) for ((_D_)=0;(_D_)<PIDX_MAX_DIMENSIONS;(_D_)++)
#define PGET(_Point_,_Coordinate_) ((&((_Point_).x))[(_Coordinate_)])

typedef struct {int x,y,z,u,v;} PointND;

unsigned int getNumBits ( unsigned int v );

uint64_t getPowerOf2(int x);

unsigned int getLevelFromBlock (int64_t block, int bits_per_block);

unsigned int getLeveL (uint64_t index);

int isValidBox(int** box);

void Deinterleave(const char* bitmask, int maxh, uint64_t zaddress, int* point);

uint64_t ZBitmask(const char* bitmask,int maxh);

uint64_t ZStart(const char* bitmask,int maxh,int BlockingH);

uint64_t ZEnd(const char* bitmask,int maxh,int BlockingH);

void ZDelta(const char* bitmask, int maxh, int BlockingH, int* point);

void GetBoxIntersection(int** a, int** b, int** c);

int** AlignEx(int** box, int* p0, int* delta);

void revstr(char* str);

void GuessBitmaskPattern(char* _bits, PointND dims);

void Align(int maxh, int H, const char* bitmask, int** userBox, int** a_offset, int** a_count, int** nsamples);

int RegExBitmaskBit(const char* bitmask_pattern,int N);

int64_t xyz_to_HZ(const char* bitmask, int maxh, PointND xyz);

void Hz_to_xyz(const char* bitmask,  int maxh, int64_t hzaddress, int64_t* xyz);

int VisusSplitFilename(const char* filename,char* dirname,char* basename);

#endif
