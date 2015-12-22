/**********************************************************************
 * PACKAGE  : voronoi
 * COPYRIGHT: (C) 2015, Ibrahim Adamou, Bernard Mourrain, Inria
 **********************************************************************/
#pragma once

#include <geometrix/mesher.hpp>
#include <geometrix/polygonizer3d.hpp>
#include <voronoi/voronoi_cell.hpp>
#include <voronoi/voronoi_controler.hpp>
#include <voronoi/distfield.hpp>
#include <voronoi/hline.h>

//====================================================================
template<class C>
struct voronoi_hline {
    typedef  mmx::mesh<C>                               Mesh;
    typedef  mmx::tmsh_vertex<3,double>                 Vertex;
<<<<<<< HEAD
    typedef  voronoi_cell                          Cell;
    typedef  distfield<hline>                      Distfield;
=======
    typedef  voronoi_cell                               Cell;
    typedef  distfield<hline>                           Distfield;
>>>>>>> 2984f6fa500a4945a67475f63558b2328403ba17
    typedef  voronoi_controler<Cell, Vertex, Distfield> Controler;
    typedef  mmx::polygonizer3d<Controler>              Polygonizer;
};

//====================================================================
