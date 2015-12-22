/**********************************************************************
 * PACKAGE  : geometrix
 * COPYRIGHT: (C) 2015, Bernard Mourrain, Inria
 **********************************************************************/
#pragma once

<<<<<<< HEAD
//#define OUTSIDE  0
//#define BOUNDARY 1
//#define BOUNDARY_XREGULAR 5   //    101
//#define BOUNDARY_YREGULAR 9   //   1001
//#define BOUNDARY_ZREGULAR 13  //   1101
//#define BOUNDARY_CENTER   17  //  10001
//#define BOUNDARY_SINGULAR 33  // 100001
//#define INSIDE   2            //     10
//#define UNKNOWN  -1

enum regularity_t {
    OUTSIDE=0,
    BOUNDARY=1,
    BOUNDARY_XREGULAR=5,   //    101
    BOUNDARY_YREGULAR=9,   //   1001
    BOUNDARY_ZREGULAR=13,  //   1101
    BOUNDARY_CENTER=17,    //  10001
    BOUNDARY_SINGULAR=33,  // 100001
    INSIDE=2,              //     10
=======
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
>>>>>>> 2984f6fa500a4945a67475f63558b2328403ba17
    UNKNOWN=-1
};

#define ISBOUNDARY(x) (x>0 && ((x&3)==1))
#define ISINSIDE(x)   (x>0 && ((x&3)==2))
