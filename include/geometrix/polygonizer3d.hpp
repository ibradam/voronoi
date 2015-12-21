/**********************************************************************
 * PACKAGE  : geometrix
 * COPYRIGHT: (C) 2015, Bernard Mourrain, Inria
 **********************************************************************/
#pragma once 

#include <vector>
#include <utility>
#include <geometrix/node.hpp>

#define TMPL  template<class CONTROLER>
#define SELF  polygonizer3d<CONTROLER>

//====================================================================
namespace mmx {
//--------------------------------------------------------------------

template<class CONTROLER >
class polygonizer3d {
public:

    typedef CONTROLER                 Controler;
    typedef typename CONTROLER::Cell  Cell;
    typedef std::vector<Cell*>        Input;
    typedef std::pair<int,int>        Edge;
    typedef std::vector<int>          Face;

    polygonizer3d(Controler* ctrl);

    ~polygonizer3d(void) ;

    void    set_regular  (Input* lc) { m_regular = lc; }
    void    set_singular (Input* lc) { m_singular = lc; }

    Input*  regular_cells (void) { return m_regular; }
    Input*  singular_cells(void) { return m_singular; }

    Controler* controler() { return d; }

    struct order {
        order(Controler* d, unsigned v): m_d(d), m_v(v) {}
        bool operator() (int i,int j) { return (m_d->vertex(i)[m_v]<m_d->vertex(j)[m_v]);}
        Controler* m_d;
        unsigned m_v;
    };

    void run(void);

    /* Remove all vertices, edges and faces from the polygonizer3d. */
    virtual void clear();

    void add_edge (int i0, int i1) { m_edge.push_back(Edge(i0,i1)); }
    void add_face (const Face& f)  { m_face.push_back(f); }

    void run_regular(Cell *cl);
    void points_on_face(int i0, int i1, int v, std::vector<int>& l);

    template<class MESH> void get_output(MESH* out);

private:
    Controler*          d;
    Input*              m_regular;
    Input*              m_singular;

public:
    std::vector< Edge > m_edge;
    std::vector< Face > m_face;

};
//--------------------------------------------------------------------
TMPL SELF::polygonizer3d(Controler* ctrl): d(ctrl)
{ 
    //m_output = new Output;
}
//--------------------------------------------------------------------
TMPL SELF::~polygonizer3d(void) {
    //    delete m_output;
}

//--------------------------------------------------------------------
TMPL void SELF::clear() {

}
//--------------------------------------------------------------------
/** @brief Extract boundary points of a face from the min and max points. 

  Search the points of type 1 on the edges between vertex i0 and vertex e.

  @param b Index of the lower point.
  @param e Index of the upper point.
  @param v Normal direction to the face.
  @param l Array where the indices of the boundary points are stored.

  @note { may go several times through the same point (to be improved). }

**/
TMPL void SELF::points_on_face(int b, int e, int v, std::vector<int>& l) {
    for(int w=0; w<3; w++)
        if(w != v) {
            int n = d->vertex(b).next(w);

            if( n != -1 && d->vertex(n)[w] <= d->vertex(e)[w] ) {

                if(d->vertex(n).tag() == 1 && std::find(l.begin(),l.end(),n) == l.end()) {
                    l.push_back(n);
                }

                this->points_on_face(n,e,v,l);

            }
        }
}

//--------------------------------------------------------------------
/** @brief Polygonization of a cell.
**/
TMPL void SELF::run_regular(Cell *cl) {
    // For each face:
    //   Extract the boundary points.
    //   Sort them by "x" reg order.
    //   Connect the adjacent points not extreme in y direction.
    //   Add the edges (pair of points ordered lexicographically) to the cell.
    // Extract the connected components to get the faces of the cell.

    int b, e;
    std::vector<int> l;
    for(unsigned v=0;v<3;v++) {
        b = cl->lowercorner(v,0);
        e = cl->uppercorner(v,0);
        points_on_face(b, e, v, l);

        std::sort(l.begin(), l.end(), order(d,v));
        
        b = cl->lowercorner(v,1);
        e = cl->uppercorner(v,1);
        points_on_face(b, e, v, l);

    }

}
//--------------------------------------------------------------------
TMPL void SELF::run() {
    foreach(Cell* cl, *this->regular_cells()) {
        run_regular(cl);
    }
}
//--------------------------------------------------------------------
TMPL
template<class MESH>
void SELF::get_output(MESH* out) {
    typedef typename MESH::Point Point;
    typedef typename MESH::Face  Face;

    for(unsigned i=0; i<d->nbv() ; i++)
        out->add_vertex(new Point(d->vertex(i)[0], d->vertex(i)[1]));
    for(unsigned i=0; i<m_edge.size() ; i++)
        out->add_edge(m_edge[i].first, m_edge[i].second);
    for(unsigned i=0; i<m_face.size() ; i++) {
        Face * f = new Face;
        for(unsigned j=0; j< m_face[i].size(); j++) {
            f->insert(m_face[i][j]);
        }
        out->add_face(f);
    }
}
//--------------------------------------------------------------------
} /* namespace mmx */
//====================================================================
#undef TMPL1
#undef TMPL
#undef SELF 

