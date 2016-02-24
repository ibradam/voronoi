/**********************************************************************
 * PACKAGE  : voronoi
 * COPYRIGHT: (C) 2015, Ibrahim Adamou, Bernard Mourrain, Inria
 **********************************************************************/
#pragma once

#include <geometrix/regularity.hpp>
#include <geometrix/tmsh_controler3.hpp>

//====================================================================
#define TMPL template<class CELL, class VERTEX, class DISTFIELD>
#define CTRL voronoi_controler<CELL,VERTEX,DISTFIELD>
//--------------------------------------------------------------------
TMPL
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
struct voronoi_controler: mmx::tmsh_controler<3,CELL,VERTEX>
{

    typedef double    C;
    typedef CELL      Cell;
    typedef VERTEX    Vertex;

    voronoi_controler() {}

    Cell* init_cell(DISTFIELD* fl,
                    double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);

    regularity_t regularity(Cell* cl);

    void vertex_tag(int n, int t) { this->vertex(n).tag(t); }

    double distance(int n, int& t);

    void process_regular (Cell* cl);
    void process_singular(Cell* cl);

    void subdivide(int v, Cell* cl, Cell*& left, Cell*& right);
    void tag_corner(Cell* cl);
    void boundary_point(Cell* cl);
    void interior_point(Cell* cl);

    std::vector<Cell*> m_regular;
    std::vector<Cell*> m_singular;

    DISTFIELD*         f;

};

//--------------------------------------------------------------------
TMPL
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

    this->mmx::tmsh_controler<3,CELL,VERTEX>::init(p);
    CELL *cl = new CELL;
    this->tag_corner(cl);

    return cl;
}
//--------------------------------------------------------------------
TMPL
double CTRL::distance(int n, int& t) {
   return f->distance2(this->vertex(n)[0], this->vertex(n)[1], this->vertex(n)[2], t);
}

//--------------------------------------------------------------------
TMPL
void CTRL::subdivide(int v, Cell* cl, Cell*& left, Cell*& right) {

    left  = new Cell(*cl);
    right = new Cell(*cl);
    this->mmx::tmsh_controler<3,CELL,VERTEX>::subdivide(v,*cl,*left,*right);

    const typename Cell::Tuple& f0 = Cell::Face[v][0];
    const typename Cell::Tuple& f1 = Cell::Face[v][1];

    // add the active sites of the corners to left and right
    for(unsigned k=0;k<Cell::Tuple::size;k++) {
        left ->add_active(this->vertex(cl->idx(f0[k])).tag());
        right->add_active(this->vertex(cl->idx(f1[k])).tag());
    }

    int t;
    for(unsigned k=0;k<Cell::Tuple::size;k++) {
        // get the closest site
        this->distance(left->idx(f1[k]), t);

        //mdebug()<<">> closest"<< t;
        // tag the vertex with the index (+1) of the closest site
        this->vertex_tag(left->idx(f1[k]),t);

        left->add_active(t);
        right->add_active(t);
    }

    /* std::cout<<"split";
    for(unsigned k=0;k<cl->nba();k++) std::cout<<" "<<cl->m_active[k]; std::cout<<std::endl;
    std::cout<<"left";
    for(unsigned k=0;k<left->nba();k++) std::cout<<" "<<left->m_active[k]; std::cout<<std::endl;
    std::cout<<"right";
    for(unsigned k=0;k<right->nba();k++) std::cout<<" "<<right->m_active[k]; std::cout<<std::endl;
    */
}

//--------------------------------------------------------------------
TMPL
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
    int a = cl->nba();

    f->site(cl->m_active[0]);

    //mdebug()<<"regularity"<<a;

    if(a<2)
        return INSIDE;
    else
        return UNKNOWN;
    /*
     if(a==2)
        return BOUNDARY_REGULAR2;
    else if(a==3)
        return BOUNDARY_REGULAR3;
    return UNKNOWN;
    */
}


//--------------------------------------------------------------------
TMPL
/*!
 * \brief tag_corner
 * \param cl: cell
 *
 *  - tag the corner of the cell with the index of the closest site.
 *  - store the distance to the closest site.
 *
 */
void CTRL::tag_corner(CELL *cl) {
  int n, t;
  for(unsigned i=0; i<8;i++) {
    n = cl->idx(i);
    f->distance2(this->vertex(n)[0], this->vertex(n)[1], this->vertex(n)[2], t);
    this->vertex_tag(n,t);
    cl->add_active(t);
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
