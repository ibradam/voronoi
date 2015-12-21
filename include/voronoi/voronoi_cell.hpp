/**********************************************************************
 * PACKAGE  : voronoi
 * COPYRIGHT: (C) 2015, Ibrahim Adamou, Bernard Mourrain, Inria
 **********************************************************************/
#pragma once

#include <geometrix/regularity.hpp>
#include <geometrix/tmsh.hpp>

//====================================================================
struct voronoi_cell : mmx::tmsh_cell<3> {

    typedef double              C;
    typedef voronoi_cell  Cell;

public:

    voronoi_cell(): tmsh_cell<3>() {}
    voronoi_cell(const voronoi_cell& cl): tmsh_cell<3>(cl) {}

    void subdivide(int v, voronoi_cell &left, voronoi_cell &right);

    void set_distance(int i, double d) { m_distance[i]=d; }

    double m_distance[8];
    std::vector<int> m_active;

};

//--------------------------------------------------------------------
void voronoi_cell::subdivide(int v, voronoi_cell & left, voronoi_cell & right) {

}
//====================================================================
