/**********************************************************************
 * PACKAGE  : geometrix
 * COPYRIGHT: (C) 2015, Bernard Mourrain, Inria
 **********************************************************************/
#pragma once

enum regularity_t
{
    OUTSIDE=0,
    BOUNDARY=1,
    INSIDE=2,              //     10
    BOUNDARY_REGULAR1=5,   //    101
    BOUNDARY_REGULAR2=9,   //   1001
    BOUNDARY_REGULAR3=13,  //   1101
    BOUNDARY_CENTER=17,    //  10001
    BOUNDARY_SINGULAR=33,  // 100001
    UNKNOWN=-1
};

#define ISBOUNDARY(x) (x>0 && ((x&3)==1))
#define ISINSIDE(x)   (x>0 && ((x&3)==2))
