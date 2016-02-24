/**********************************************************************
 * PACKAGE  : geometrix
 * COPYRIGHT: (C) 2015, Bernard Mourrain, Inria
 **********************************************************************/
#pragma once

#include <vector>
#include <utility>
#include <geometrix/mdebug.hpp>
#include <geometrix/node.hpp>
#include <geometrix/tmsh_cell.hpp>

#define TMPL  template<class CONTROLER>
#define SELF  polygonizer3d_mc<CONTROLER>

//====================================================================
namespace mmx {

extern int marching_cube_tri_table[256][16] ;

//--------------------------------------------------------------------
template<class CONTROLER >
class polygonizer3d_mc {
public:

    typedef CONTROLER                 Controler;
    typedef typename CONTROLER::Cell  Cell;
    typedef std::vector<Cell*>        Input;
    typedef std::pair<int,int>        Edge;
    typedef std::vector<int>          Face;

    polygonizer3d_mc(void);
    polygonizer3d_mc(Controler* ctrl);

    ~polygonizer3d_mc(void) ;

    Input*  regular_cells (void) { return m_regular; }
    Input*  singular_cells(void) { return m_singular; }

    Controler* controler() { return d; }

    void run(void);

    /* Remove all vertices, edges and faces from the polygonizer3d. */
    virtual void clear();

    void set_regular  (Input* lc) { m_regular = lc; }
    void set_singular (Input* lc) { m_singular = lc; }

    void add_edge (int i0, int i1) { m_edge.push_back(Edge(i0,i1)); }
    void add_face (const Face& f)  { m_face.push_back(f); }

    int  mc_index(Cell *cl);
    int  mc_edge_point(Cell *cl, int n0, int n1, int v);
    void mc_edge_points(std::vector<int>& edges, Cell *cl, int cubeindex);

    void face_regular(Cell *cl);

    template<class MESH> void get_mesh(MESH* out);

private:
    Controler*          d;
    Input*              m_regular;
    Input*              m_singular;

public:
    std::vector< Edge > m_edge;
    std::vector< Face > m_face;

};
//--------------------------------------------------------------------
TMPL SELF::polygonizer3d_mc(void): d(new Controler), m_regular(NULL), m_singular(NULL)
{

}

//--------------------------------------------------------------------
TMPL SELF::polygonizer3d_mc(Controler* ctrl): d(ctrl)
{
    //m_output = new Output;
}
//--------------------------------------------------------------------
TMPL SELF::~polygonizer3d_mc(void)
{
    //    delete m_output;
}

//--------------------------------------------------------------------
TMPL void SELF::clear() {

}

//--------------------------------------------------------------------
TMPL
/*!
 * \brief SELF::mc_index
 * \param cl: cell
 * \return the marching cube index of the cell cl.
 */
int SELF::mc_index(Cell *cl) {
    int l = -1;
    for(unsigned k=0; k<8; k++) {
        //mdebug()<<"  "<<k<< "tag "<<this->controler()->vertex(cl->idx(k)).tag();
        l = std::max(l, this->controler()->vertex(cl->idx(k)).tag());
    }

    //mdebug()<<"max tag"<<l;
    int CubeIndex=0;

    /*int s=1;
    for (unsigned i=0; i< 8; i++) {
        if (this->controler()->vertex(cl->idx(i)).tag() < l) CubeIndex |= s;
        s*=2;
    }
    */

    if (this->controler()->vertex(cl->idx(0)).tag() < l) CubeIndex |= 1;
    if (this->controler()->vertex(cl->idx(1)).tag() < l) CubeIndex |= 2;
    if (this->controler()->vertex(cl->idx(3)).tag() < l) CubeIndex |= 4;
    if (this->controler()->vertex(cl->idx(2)).tag() < l) CubeIndex |= 8;
    if (this->controler()->vertex(cl->idx(4)).tag() < l) CubeIndex |= 16;
    if (this->controler()->vertex(cl->idx(5)).tag() < l) CubeIndex |= 32;
    if (this->controler()->vertex(cl->idx(7)).tag() < l) CubeIndex |= 64;
    if (this->controler()->vertex(cl->idx(6)).tag() < l) CubeIndex |= 128;

    return CubeIndex;
}
//--------------------------------------------------------------------
TMPL
/*!
 * \brief SELF::mc_edge_point
 * \param cl: cell
 * \param n0: index of first point
 * \param n1: index of last point
 * \param v: direction
 *
 * \return the index of the first point with tag -1 between n0 and n1.
 *
 * If it does not exists, returns n1;
 */
int SELF::mc_edge_point(Cell *cl, int c0, int c1, int v) {

    int n = cl->idx(c0), n1 = cl->idx(c1), w = v;

    //    mdebug()<<c0<<"->"<<c1<<" v:"<<v<< "tag:"<<this->controler()->vertex(n).tag()
    //           <<this->controler()->vertex(n1).tag();
    //    mdebug()<<" -->"<<n;
    while(n != n1 && n != -1 && this->controler()->vertex(n).tag() != -1) {
        // mdebug()<<" -->"<<n<< "tag:"<<this->controler()->vertex(n).tag();
        n = this->controler()->vertex(n).next(v);
    }

    if(n == n1 && this->controler()->vertex(n).tag() != -1)
        mdebug()<<"mc_edge_point not found in"<<cl->idx(c0)<<"--"<<n1;
    //mdebug()<<" ==>"<<n;
    return n;
}

//--------------------------------------------------------------------
TMPL
/*!
 * \brief SELF::mc_edge_points
 * \param edges: array of point index on edges
 * \param cl: cell
 * \param cubeindex: marching cube index of the cell
 *
 * Store in edges the indices of the vertices with tag -1 on the edges.
 */
void SELF::mc_edge_points(std::vector<int>& edges, Cell *cl, int cubeindex) {

    //mdebug()<<"mc_index"<<cubeindex;
    if(cubeindex==0 || cubeindex==255)
    //if (marching_cube_edge_table[cubeindex] == 0)
        return;

    int c0, c1, v;
    for(unsigned i=0;i<12; i++) {
        c0 = tmsh_cell<3>::Edge[i][0];
        c1 = tmsh_cell<3>::Edge[i][1];
        v  = tmsh_cell<3>::EdgeDir[i];
        if(this->controler()->vertex(cl->idx(c0)).tag() != this->controler()->vertex(cl->idx(c1)).tag() ) {
            edges[i] = mc_edge_point(cl,c0,c1,v);
            //mdebug()<<"edge"<<i<<" "<<c0<<c1 << "index "<<edges[i];
        }
    }

}

//--------------------------------------------------------------------
/** @brief Polygonization of a cell.
**/
TMPL void SELF::face_regular(Cell *cl) {

    int cubeindex = mc_index(cl);
    std::vector<int> edges(12);
    mc_edge_points(edges, cl, cubeindex);
    int f;
    for (int i=0; marching_cube_tri_table[cubeindex][i]!=-1; i+=3) {
        Face face;
        //        mdebug()<<"triangle"<< cubeindex<< ": "
        //          <<edges[marching_cube_tri_table[cubeindex][i]]
        //          <<edges[marching_cube_tri_table[cubeindex][i+1]]
        //          <<edges[marching_cube_tri_table[cubeindex][i+2]];

        face.push_back(edges[marching_cube_tri_table[cubeindex][i]]);
        face.push_back(edges[marching_cube_tri_table[cubeindex][i+1]]);
        face.push_back(edges[marching_cube_tri_table[cubeindex][i+2]]);
        this->add_face(face);
    }
}

//--------------------------------------------------------------------
TMPL void SELF::run() {
    foreach(Cell* cl, *this->regular_cells()) {
        face_regular(cl);
    }
}
//--------------------------------------------------------------------
TMPL
template<class MESH>
        void SELF::get_mesh(MESH* out) {
    typedef typename MESH::Point Point;
    typedef typename MESH::Face  Face;

    for(unsigned i=0; i<d->nbv() ; i++)
        out->add_vertex(new Point(d->vertex(i)[0], d->vertex(i)[1],d->vertex(i)[2]));
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
#undef TMPL
#undef SELF

