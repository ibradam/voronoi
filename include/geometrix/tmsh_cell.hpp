/**********************************************************************
 * PACKAGE  : geometrix
 * COPYRIGHT: (C) 2015, Bernard Mourrain, Inria
 **********************************************************************/
#pragma once

#include <geometrix/tuple.hpp>

//====================================================================
namespace mmx {
//--------------------------------------------------------------------

template<int N> struct tmsh_cell;

template<>
struct tmsh_cell<2> {
    typedef tuple<2,int> Tuple;

    tmsh_cell(void): m_reg(0) {
        for(unsigned i=0;i<4;i++) corners[i]=i;
    }

    tmsh_cell(const tmsh_cell<2>& cel): m_reg(cel.m_reg) {
       for(unsigned i=0;i<4;i++) corners[i]=cel[i];
    }

    tmsh_cell(int idx[4]): m_reg(0) {
        for(unsigned i=0;i<4;i++) corners[i]=idx[i];
    }

    int  operator[](int i) const { return corners[i]; }
    int& operator[](int i)       { return corners[i]; }

    int  idx(int i) const { return corners[i]; }
    int& idx(int i)       { return corners[i]; }

    void regularity(int r)       { m_reg = r; }
    int  regularity(void)  const { return m_reg; }

    void add_center(int i)       { m_center.push_back(i); }
    int  center(int i)     const { return (m_center[i]); }
    bool has_center(void)  const { return (m_center.size() > 0); }

    const std::vector<int>& boundary(void) const  { return m_boundary;}
    std::vector<int>&       boundary(void)        { return m_boundary;}
    int                     boundary(int k) const { return m_boundary[k]; }

    void add_boundary_point(int i, int mu) {
        if(std::find(m_boundary.begin(), m_boundary.end(),i) == m_boundary.end())
            m_boundary.push_back(i);
        if(mu>1) m_boundary.push_back(i);
    }

    static unsigned dim;

public:
    int  corners[4];

    static tuple<2,int>  Face[2][2];
    static tuple<2,int>  FaceEdge[1];
    static unsigned      nb_fedge;


private:
    int              m_reg;
    std::vector<int> m_boundary;
    std::vector<int> m_center;
};

unsigned tmsh_cell<2>::dim = 3;
unsigned tmsh_cell<2>::nb_fedge = 1;

tuple<2,int> tmsh_cell<2>::FaceEdge[1] = {
    tuple<2,int>(0,1)
};

tuple<2,int> tmsh_cell<2>::Face[2][2]={
    {tuple<2,int>(0,2), tuple<2,int>(1,3)},
    {tuple<2,int>(0,1), tuple<2,int>(2,3)}
};

//--------------------------------------------------------------------

//--------------------------------------------------------------------
template<>
/**
 * @brief The tmsh_cell<3> struct
 * The directions, vertices and edges are numbered as follows:
 *
 * \begincode
 *    v=2           v=1
 *     |            /
 *     |           /
 *     |   6 ----6---- 7
 *     |  /|     /    /|
 *     | 7 |    /    5 |
 *     |/  11  /    /  10
 *     4 ----4---- 5   |
 *     |   | /     |   |
 *     |   2 ----2-|-- 3
 *     8  /        9  /
 *     | 3         | 1
 *     |/          |/
 *     0 ----0---- 1 --------- v = 0
 * \endcode
 */
struct tmsh_cell<3> {
    typedef tuple<4,int> Tuple;

    tmsh_cell(void) {
        for(unsigned i=0;i<8;i++) corners[i]=i;
    }

    tmsh_cell(int idx[8]) {
        for(unsigned i=0;i<8;i++) corners[i]=idx[i];
    }

    int  operator[](int i) const { return corners[i]; }
    int& operator[](int i)       { return corners[i]; }

    int  idx(int i)        const { return corners[i]; }
    int& idx(int i)              { return corners[i]; }

    int  lowercorner(int v, int s) const { return corners[Face[v][s][0]]; }
    int  uppercorner(int v, int s) const { return corners[Face[v][s][3]]; }

    int corners[8];

    static unsigned dim;
    static tuple<2,int> Edge[12];
    static int EdgeDir[12];
    static tuple<4,int> Face[3][2];
    static tuple<2,int> FaceEdge[4];
    static unsigned nb_fedge;
};

unsigned tmsh_cell<3>::dim = 3;

tuple<2,int> tmsh_cell<3>::Edge[12] = {
    tuple<2,int>(0,1), tuple<2,int>(1,3), tuple<2,int>(2,3), tuple<2,int>(0,2),
    tuple<2,int>(4,5), tuple<2,int>(5,7), tuple<2,int>(6,7), tuple<2,int>(4,6),
    tuple<2,int>(0,4), tuple<2,int>(1,5), tuple<2,int>(3,7), tuple<2,int>(2,6)
};

int tmsh_cell<3>::EdgeDir[12] = {
    0,1,0,1,
    0,1,0,1,
    2,2,2,2
};

tuple<4,int> tmsh_cell<3>::Face[3][2]={
    {tuple<4,int>(0,2,4,6), tuple<4,int>(1,3,5,7)},
    {tuple<4,int>(0,1,4,5), tuple<4,int>(2,3,6,7)},
    {tuple<4,int>(0,1,2,3), tuple<4,int>(4,5,6,7)}
};

unsigned tmsh_cell<3>::nb_fedge = 4;

tuple<2,int> tmsh_cell<3>::FaceEdge[4] = {
    tuple<2,int>(0,1), tuple<2,int>(2,3), tuple<2,int>(0,2), tuple<2,int>(1,3)
};

//--------------------------------------------------------------------
} /* namespace mmx */
//====================================================================
#undef TMPL
