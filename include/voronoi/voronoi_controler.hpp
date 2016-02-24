/**********************************************************************
 * PACKAGE  : voronoi
 * COPYRIGHT: (C) 2015, Ibrahim Adamou, Bernard Mourrain, Inria
 **********************************************************************/
#pragma once

#include <voronoi/mediatrice.h>
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

    //void vertex_tag(int n, int t) { this->vertex(n).tag(t); }

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
/*!
 * \brief CTRL::distance
 * \param n
 * \param t
 * \return
 */
double CTRL::distance(int n, int& t) {
    return f->distance2(this->vertex(n)[0], this->vertex(n)[1], this->vertex(n)[2], t);
}

//--------------------------------------------------------------------
TMPL
/*!
 * \brief CTRL::subdivide
 * \param v
 * \param cl
 * \param left
 * \param right
 */
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
        this->vertex(left->idx(f1[k])).tag(t);

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
 *   - BOUNDARY: 2,3,4 active sites
 *   _ UNKNOWN: otherwise
 */
regularity_t CTRL::regularity(Cell *cl) {
    int a = cl->nba();

    //mdebug()<<"regularity"<<a;

    if(a==1)
        return INSIDE;
    else  if(a==2)
        return BOUNDARY;
    else
        return UNKNOWN;

    //-----------------------------------------------------------
    // To be used later

    mmx::tmsh_vertex < 3, double> F0 = this-> vertex(cl-> idx(0));
    mmx::tmsh_vertex < 3, double> F1 = this-> vertex(cl-> idx(1));
    mmx::tmsh_vertex < 3, double> F2 = this-> vertex(cl-> idx(2));
    mmx::tmsh_vertex < 3, double> F3 = this-> vertex(cl-> idx(3));
    mmx::tmsh_vertex < 3, double> F4 = this-> vertex(cl-> idx(4));
    mmx::tmsh_vertex < 3, double> F5 = this-> vertex(cl-> idx(5));
    mmx::tmsh_vertex < 3, double> F6 = this-> vertex(cl-> idx(6));
    mmx::tmsh_vertex < 3, double> F7 = this-> vertex(cl-> idx(7));
    double m1,m2,m3,M1,M2,M3;
    m1=F0[0];m2=F0[1];m3=F0[2];
    M1=F7[0];M2=F7[1];M3=F7[2];

    if(a==1)
        return INSIDE;

    if(a==2)
    {
        for(unsigned v=0; v<3; v++ ) {
            for(unsigned s=0; s<2; s++) {
                unsigned c=0;
                if (this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][0])).tag() !=
                        this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][1])).tag() ) c++;
                if (this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][1])).tag() !=
                        this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][3])).tag() ) c++;
                if (this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][3])).tag() !=
                        this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][2])).tag() ) c++;
                if (this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][2])).tag() !=
                        this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][0])).tag() ) c++;
                if (c>2) return UNKNOWN;
            }
        }
        return BOUNDARY;

    }

    else  if(a==3)
    {
        for(unsigned v=0; v<3; v++ ) {
            for(unsigned s=0; s<2; s++) {
                unsigned c=0;
                if (this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][0])).tag() !=
                        this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][1])).tag() ) c++;
                if (this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][1])).tag() !=
                        this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][3])).tag() ) c++;
                if (this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][3])).tag() !=
                        this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][2])).tag() ) c++;
                if (this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][2])).tag() !=
                        this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][0])).tag() ) c++;
                if (c>3) return UNKNOWN;
            }
        }
        return BOUNDARY;
    }

    else  if(a==4)
    {
         int info=0 ;

        mmx::point<double> q=equidist( f->site(cl->m_active[0]),f->site(cl->m_active[1]), f->site(cl->m_active[2]),f->site(cl->m_active[3]),info);

        for(unsigned v=0; v<3; v++ ) {
            for(unsigned s=0; s<2; s++) {
                unsigned c=0;
                if (this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][0])).tag() !=
                        this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][1])).tag() ) c++;
                if (this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][1])).tag() !=
                        this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][3])).tag() ) c++;
                if (this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][3])).tag() !=
                        this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][2])).tag() ) c++;
                if (this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][2])).tag() !=
                        this->vertex(cl->idx(mmx::tmsh_cell<3>::Face[v][s][0])).tag() ) c++;
                if (c>4 || (q[0]>m1 && q[0]<M1 && q[1]>m2 && q[1]<M2 && q[2]>m3 && q[2]<M3) ) return UNKNOWN;
            }
        }
        return BOUNDARY;
    }

    if(a>5) return UNKNOWN;

    return UNKNOWN;
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
        this->vertex(n).tag(t);
        cl->add_active(t);
    }
}

//--------------------------------------------------------------------
TMPL
void CTRL::boundary_point(Cell* cl) {

    int n;
    mmx::point<double> A(0,0,0), B(0,0,0), m(0,0,0);

    for(int i=0; i<12; i++) {
        int n0 = cl->idx(mmx::tmsh_cell<3>::Edge[i][0]),
            n1 = cl->idx(mmx::tmsh_cell<3>::Edge[i][1]);
        int v = mmx::tmsh_cell<3>::EdgeDir[i];
        if( this->vertex(n0).tag() != this->vertex(n1).tag()
                && this->vertex(n0).tag() != -1
                && this->vertex(n1).tag() != -1) {
            //mdebug()<<"distinct"<<n0<<n1;
            for(unsigned k=0;k<3;k++) {
                A[k] = this->vertex(n0)[k];
                B[k] = this->vertex(n1)[k];
            }

            mdebug()<<"equidist"<< this->vertex(n0).tag() << this->vertex(n1).tag()
                   << ":" << f->site(this->vertex(n0).tag())<<","<< f->site(this->vertex(n1).tag())
                   << "\nin segment"<< A<<"-"<<B;
            m = equidist(f->site(this->vertex(n0).tag()), f->site(this->vertex(n1).tag()), A, B);
            mdebug() <<" --> " << m;
            n = this->insert_vertex( Vertex(m[0],m[1],m[2]), n0, n1, v);
            mdebug() <<"inserted point" <<n<< m<<"\n";
            this->vertex(n).tag(-1);
        }
    }

}
//--------------------------------------------------------------------
TMPL
void CTRL::interior_point(Cell* cl) {

}

//--------------------------------------------------------------------
TMPL
void CTRL::process_regular(Cell* cl) {
    tag_corner(cl);
    boundary_point(cl);
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
