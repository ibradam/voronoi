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
<<<<<<< HEAD
=======
/*!
 * \brief The voronoi_controler struct
 *
 * controls the subdvision process by the following methods:
 *  - init_cell: create the initial cell of the subdivision;
 *  - regularity: determine the regularity of a cell;
 *  - subdivide: split the cell;
 *  - process_regular: process the regular cells;
 *  - process_singular: process the singular cells;
 */
>>>>>>> 2984f6fa500a4945a67475f63558b2328403ba17
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
<<<<<<< HEAD

=======
>>>>>>> 2984f6fa500a4945a67475f63558b2328403ba17
    void interior_point(Cell* cl);

    std::vector<Cell*> m_regular;
    std::vector<Cell*> m_singular;

    DISTFIELD*         f;

};

//--------------------------------------------------------------------
TMPL
<<<<<<< HEAD
=======
/*!
 * \brief init_cell
 * \param fld
 * \param xmin
 * \param xmax
 * \param ymin
 * \param ymax
 * \param zmin
 * \param zmax
 * \return a cell
 */
>>>>>>> 2984f6fa500a4945a67475f63558b2328403ba17
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
<<<<<<< HEAD
regularity_t CTRL::regularity(Cell *cl) {
=======
/*!
 * \brief regularity
 * \param cl: cell
 * \return the regularity of the cell
 * Its regularity is
 *   - INSIDE: one active site
 *   - BOUNDARY: 2,3 active sites
 *   _ UNKNOWN: otherwise
 */
regularity_t CTRL::regularity(Cell *cl) {
    int a = cl->m_active_site.size();
    if(a==1)
        return INSIDE;
    else if(a==2)
        return BOUNDARY_REGULAR2;
    else if(a==3)
        return BOUNDARY_REGULAR3;
    return UNKNOWN;
>>>>>>> 2984f6fa500a4945a67475f63558b2328403ba17

}

//--------------------------------------------------------------------
TMPL
<<<<<<< HEAD
=======
/*!
 * \brief tag_corner
 * \param cl: cell
 *
 *  - tag the corner of the cell with the index of the closest site.
 *  - store the distance to the closest site.
 *
 */
>>>>>>> 2984f6fa500a4945a67475f63558b2328403ba17
void CTRL::tag_corner(CELL *cl) {
  int n, t;
  double d;
  for(unsigned i=0; i<8;i++) {
    n = cl->idx(i);
<<<<<<< HEAD
    d = f->distance2(this->vertex(n)[0],this->vertex(n)[1],this->vertex(n)[2],t);
=======
    d = f->distance2(this->vertex(n)[0], this->vertex(n)[1], this->vertex(n)[2], t);
>>>>>>> 2984f6fa500a4945a67475f63558b2328403ba17
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
