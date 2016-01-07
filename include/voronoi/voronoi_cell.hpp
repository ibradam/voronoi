/**********************************************************************
 * PACKAGE  : voronoi
 * COPYRIGHT: (C) 2015, Ibrahim Adamou, Bernard Mourrain, Inria
 **********************************************************************/
#pragma once

#include <geometrix/regularity.hpp>
#include <geometrix/tmsh_cell.hpp>

//====================================================================
/*!
 * \brief The voronoi_cell struct
 */
struct voronoi_cell : mmx::tmsh_cell<3> {

    typedef double              C;
    typedef voronoi_cell  Cell;

public:

    voronoi_cell(): tmsh_cell<3>() {}
    voronoi_cell(const voronoi_cell& cl): tmsh_cell<3>(cl) {}

    void add_active(int i) {
        //mdebug()<<"add active"<<i<< "size"<<m_active.size();
        if(std::find(m_active.begin(), m_active.end(),i) == m_active.end())
            m_active.push_back(i);
        //mdebug()<<">add active"<<i<< "size"<<m_active.size();
    }

    unsigned nba () const { return m_active.size(); }

    std::vector<int> m_active;

};
//====================================================================
