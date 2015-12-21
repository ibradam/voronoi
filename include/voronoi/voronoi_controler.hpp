/**********************************************************************
 * PACKAGE  : voronoi
 * COPYRIGHT: (C) 2015, Ibrahim Adamou, Bernard Mourrain, Inria
 **********************************************************************/
#pragma once

#include <geometrix/regularity.hpp>
#include <geometrix/tmsh.hpp>

//====================================================================
#define TMPL template<class CELL, class VERTEX, class DISTFIELD>
#define CTRL voronoi_controler<CELL,VERTEX,DISTFIELD>
//--------------------------------------------------------------------
TMPL
struct voronoi_controler: mmx::tmsh<CELL,VERTEX>
{

    typedef double    C;
    typedef CELL      Cell;
    typedef VERTEX    Vertex;

    voronoi_controler() {}

    Cell* init_cell(DISTFIELD* fl,
                    double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);

    regularity_t regularity(Cell* cl);

    void process_regular (Cell* cl);
    void process_singular(Cell* cl);

    void subdivide(int v, Cell* cl, Cell*& left, Cell*& right) {
        left  = new Cell(*cl);
        right = new Cell(*cl);
        this->mmx::tmsh<CELL,VERTEX>::subdivide(v,*cl,*left,*right);
        cl->subdivide(v,*left,*right);
    }

    void tag_corner(Cell* cl);

    void boundary_point(Cell* cl);

    void interior_point(Cell* cl);

    std::vector<Cell*> m_regular;
    std::vector<Cell*> m_singular;

    DISTFIELD*         f;

};

//--------------------------------------------------------------------
TMPL
CELL* CTRL::init_cell(DISTFIELD* fld,
                      double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{

    f=fld;
    Vertex p[8] = {
        Vertex(xmin,ymin,zmin),
        Vertex(xmax,ymin,zmin),
        Vertex(xmin,ymax,zmin),
        Vertex(xmax,ymax,zmin),
        Vertex(xmin,ymin,zmax),
        Vertex(xmax,ymin,zmax),
        Vertex(xmin,ymax,zmax),
        Vertex(xmax,ymax,zmax)
    };

    this->mmx::tmsh<CELL,VERTEX>::init(p);

    return new CELL;
}

//--------------------------------------------------------------------
TMPL
regularity_t CTRL::regularity(Cell *cl) {

}

//--------------------------------------------------------------------
TMPL
void CTRL::tag_corner(CELL *cl) {
  int n, t;
  double d;
  for(unsigned i=0; i<8;i++) {
    n = cl->idx(i);
    d = f->distance2(this->vertex(n)[0],this->vertex(n)[1],this->vertex(n)[2],t);
    this->vertex(n).tag(t);
    cl->set_distance(i,d);
  }
}

//--------------------------------------------------------------------
TMPL
void CTRL::boundary_point(Cell* cl) {

}

//--------------------------------------------------------------------
TMPL
void CTRL::interior_point(Cell* cl) {

}

//--------------------------------------------------------------------
TMPL
void CTRL::process_regular(Cell* cl) {
    m_regular.push_back(cl);
}

//--------------------------------------------------------------------
TMPL
void CTRL::process_singular(Cell* cl) {
    m_singular.push_back(cl);
}

//====================================================================
#undef TMPL
#undef CTRL
