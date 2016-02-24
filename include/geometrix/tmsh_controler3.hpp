/**********************************************************************
 * PACKAGE  : geometrix
 * COPYRIGHT: (C) 2015, Bernard Mourrain, Inria
 **********************************************************************/
#pragma once

#include <vector>
#include <iostream>

#include <geometrix/mdebug.hpp>
#include <geometrix/tmsh_cell.hpp>
#include <geometrix/tmsh_vertex.hpp>

#define TMPL template <class CELL, class VERTEX>
#define SELF tmsh_controler<3,CELL,VERTEX>

//====================================================================
namespace mmx {
//--------------------------------------------------------------------
template <int N, class CELL, class VERTEX> struct tmsh_controler;

template <class CELL, class VERTEX>
struct tmsh_controler<3,CELL,VERTEX> {

    typedef VERTEX Vertex;
    typedef CELL   Cell;

    tmsh_controler(): m_eps(1e-6) {}

    const Vertex&  vertex(int i) const { return m_vertices[i]; }
    Vertex&        vertex(int i)       { return m_vertices[i]; }

    void init(VERTEX cl[]);

    double  min(int v, const CELL& c);
    double  max(int v, const CELL& c);

    double xmin(const CELL& c);
    double ymin(const CELL& c);
    double zmin(const CELL& c);

    double xmax(const CELL& c);
    double ymax(const CELL& c);
    double zmax(const CELL& c);

    double xsize(const CELL& c);
    double ysize(const CELL& c);
    double zsize(const CELL& c);
    double  size(const CELL& c);
    double  size(CELL* c) { return this->size(*c); }

    int  neighbor(int idx, int& w);
    int  neighbor(int idx, int &x0, int &x1, int &w);

    int  split_direction(const Cell& cl);
    //int  split_direction3d(const Cell& cl);

    int  insert_vertex(const Vertex &p);
    int  insert_vertex(const Vertex& p, int i0, int i1 , int v);
    int  insert_vertex(const Vertex& p, int i0, int i1);
    int  insert_middle(int i0, int i1);

    void insert_edge  (int i0, int i1, int v);
    void insert_edge  (int i0, int i1);

    int  insert_cell  (const Cell& c);

    void subdivide(int v, const CELL& c0, CELL& left, CELL& right);

    unsigned nbv() const { return m_vertices.size(); }
    unsigned nbc() const { return m_cells.size(); }

    bool is_identifiable(double a, double b) const {
        return (a-b<m_eps && a-b>-m_eps);
    }

    void set_tag_corner(CELL* cl, int t);

private:
    double              m_eps;
    std::vector<VERTEX> m_vertices;
    std::vector<CELL>   m_cells;
};


TMPL
void SELF::init(VERTEX p[]) {
    for(unsigned k=0;k< 2*CELL::Tuple::size; k++)
        this->insert_vertex(p[k]);
    for(unsigned v=0;v< VERTEX::dim;v++) {
        for(unsigned k=0;k< CELL::Tuple::size; k++)
            this->insert_edge(CELL::Face[v][0][k], CELL::Face[v][1][k], v);
    }
}

TMPL
double SELF::min(int v, const CELL &cl) {
    return (this->vertex(cl[0])[v]);
}

TMPL
double SELF::max(int v, const CELL &cl) {
    return (this->vertex(cl[2*v+1])[v]);
}

TMPL
double SELF::xmin(const CELL &cl) {
    return (this->vertex(cl[0])[0]);
}

TMPL
double SELF::ymin(const CELL &cl) {
    return (this->vertex(cl[0])[1]);
}

TMPL
double SELF::zmin(const CELL &cl) {
    return (this->vertex(cl[0])[2]);
}

TMPL
double SELF::xmax(const CELL &cl) {
    return (this->vertex(cl[1])[0]);
}

TMPL
double SELF::ymax(const CELL &cl) {
    return (this->vertex(cl[2])[1]);
}

TMPL
double SELF::zmax(const CELL &cl) {
    return (this->vertex(cl[4])[2]);
}

TMPL
double SELF::xsize(const CELL &cl) {
    return (this->vertex(cl[1])[0]-this->vertex(cl[0])[0]);
}

TMPL
double SELF::ysize(const CELL &cl) {
    return (this->vertex(cl[2])[1]-this->vertex(cl[0])[1]);
}

TMPL
double SELF::zsize(const CELL &cl) {
    return (this->vertex(cl[4])[2]-this->vertex(cl[0])[2]);
}

TMPL
double SELF::size(const CELL& c) {
    return std::max(this->xsize(c),this->ysize(c));
}

/** @brief neighbor on the cell boundary of side w
 *
 * Compute the next vertex index in direction w and update w if
 * there is an element in the direction w+1 mod 4;
**/
TMPL
int SELF::neighbor(int idx, int& w) {
    int x0 = w%2, x1=(w+1)%2;
    return neighbor(idx, x0,x1, w);
}

/**
 * @brief neighbor on the face of dimensions x0 and x1
 */
TMPL
int SELF::neighbor(int idx, int &x0, int &x1, int &w) {
    if(idx == -1) return -1;
    int d0 = w/2;

    int n = (d0==0?this->vertex(idx).next(x0):this->vertex(idx).previous(x0));
    if (n != -1) {
        int d1 = ((w+1)%4)/2;
        int n1 = (d1==0?this->vertex(n).next(x1):this->vertex(n).previous(x1));

        if(n1 != -1) {
            w++;
            w%=4;
            int tmp = x0;
            x0 = x1;
            x1 = tmp;
        }
    }
    return n;
}

TMPL
int SELF::split_direction(const CELL& cl)
{
    double s[3] = { this->xsize(cl), this->ysize(cl), this->zsize(cl) };
    int v = 0;
    if (s[1]>s[0]) v = 1;
    if (s[2]>s[v]) v = 2;
    return v;
}

TMPL
void SELF::subdivide(int v, const CELL &c0, CELL &left, CELL& right) {

    Vertex p;

    int i0, i1;
    const typename CELL::Tuple& f0 = CELL::Face[v][0];
    const typename CELL::Tuple& f1 = CELL::Face[v][1];
    int nwi[4];
    for(unsigned k=0;k<CELL::Tuple::size;k++) {

        i0 = c0[f0[k]];
        i1 = c0[f1[k]];

        set_middle(p, this->vertex(i0), this->vertex(i1));

        nwi[k] = this->insert_vertex(p,i0,i1,v);

        left [f0[k]]=c0[f0[k]];
        right[f1[k]]=c0[f1[k]];

        left [f1[k]]=nwi[k];
        right[f0[k]]=nwi[k];

    }

    for(unsigned k=0;k<CELL::nb_fedge;k++)
        insert_edge(nwi[CELL::FaceEdge[k][0]],nwi[CELL::FaceEdge[k][1]]);

}

TMPL
int SELF::insert_vertex(const VERTEX &p) {
    int i=m_vertices.size();
    m_vertices.push_back(p);
    return i;
}

TMPL
int SELF::insert_vertex(const VERTEX &p, int i0, int i1) {
    int v=(vertex(i0)[0]!=vertex(i1)[0]?0:(vertex(i0)[1]!=vertex(i1)[1]?1:2));
    return insert_vertex(p, i0, i1, v);
}

TMPL
int SELF::insert_vertex(const VERTEX &p, int i0, int i1, int v) {
    // mdebug()<<"insert_vertex v:"<<v << " "<<i0<<i1<<"   "<<p[0]<<p[1]
    //  <<" ["<< vertex(i0)[0]<<vertex(i0)[1]
    //  <<"--"<< vertex(i1)[0]<<vertex(i1)[1]
    //  <<"]";

    if(this->is_identifiable(p[v],vertex(i0)[v]))
        return i0;
    if(this->is_identifiable(p[v],vertex(i1)[v]))
        return i1;

    int k, k0=i0, k1 = vertex(i0).next(v);
    // mdebug()<<"insert loop"<<k0<<k1;
    while(vertex(k1).next(v)>=0 && p[v]>vertex(k1)[v]+m_eps) {
        k0 = k1;
        k  = k1;
        k1 = vertex(k).next(v);
        // mdebug()<<"insert loop"<<k0<<k1;
    }

    int j=(this->is_identifiable(p[v],vertex(k1)[v])?k1:insert_vertex(p));

    int j1 = (j==k1?vertex(k1).next(v):k1);
    //mdebug()<<"insert"<<i0<<i1<<"  "<<k0<<k1<<"  "<<j<<j1;
    m_vertices[k0].m_neighbor[2*v]   = j;
    m_vertices[j1].m_neighbor[2*v+1] = j;
    m_vertices [j].m_neighbor[2*v+1] = k0;
    m_vertices [j].m_neighbor[2*v]   = j1;

    // mdebug()<<"insert idx"<<j
    // <<vertex(j).ngbr(0,0)<< vertex(j).ngbr(0,1)
    // <<vertex(j).ngbr(1,0)<< vertex(j).ngbr(1,1);
    return j;
}

TMPL
int SELF::insert_middle(int i0, int i1) {
    Vertex p;
    for(unsigned i=0;i<Vertex::dim;i++) {
        p[i] = (this->vertex(i0)[i]+this->vertex(i1)[i])/2;
    }
    return insert_vertex(p,i0,i1);
}


TMPL
void SELF::insert_edge(int i0, int i1, int v) {
    m_vertices[i0].m_neighbor[2*v]=i1;
    m_vertices[i1].m_neighbor[2*v+1]=i0;
}

TMPL
void SELF::insert_edge(int i0, int i1) {
    int v=(vertex(i0)[0]!=vertex(i1)[0]?0:(vertex(i0)[1]!=vertex(i1)[1]?1:2));
    insert_edge(i0, i1, v);
}

TMPL
int SELF::insert_cell(const CELL &c) {
    int i=m_cells.size();
    m_cells.push_back(c);
    return i;
}

TMPL
void SELF::set_tag_corner(CELL* cl, int t) {
    for(unsigned k=0; k<4; k++) {
        this->vertex(cl->vertex_idx(k)).tag(t);
    }
}

template<class OSTREAM, class CELL, class VERTEX>
OSTREAM& operator<<(OSTREAM& os, const SELF& msh) {
    os<<"<mesh size=\"0.1\" color=\"0 0 255\">\n";
    unsigned c=0;
    for(unsigned i=0;i<msh.nbv();i++) {
        for(unsigned j=0; j<SELF::Vertex::dim;j++)
            if(msh.vertex(i).m_neighbor[2*j]>=0)
                c++;
    }
    os<<"<count>"<< msh.nbv()<<" "<<c<<" 0</count>\n";
    os<<"<points>\n";
    for(unsigned i=0;i<msh.nbv();i++) {
        os<< msh.vertex(i)[0]<<" "<<msh.vertex(i)[1]<<" "<<msh.vertex(i)[2]<<"\n";
    }
    os<<"</points>\n";

    os<<"<edges>\n";
    for(unsigned i=0;i<msh.nbv();i++) {
        for(unsigned j=0; j<SELF::Vertex::dim;j++)
            if(msh.vertex(i).m_neighbor[2*j]>=0)
                os<<"2 "<<i<<" "<< msh.vertex(i).m_neighbor[2*j]<<"\n";
    }
    os<<"</edges>\n";
    os<<"</mesh>\n";
    return os;
}
//--------------------------------------------------------------------
} /* namespace mmx */
//====================================================================
#undef SELF
#undef TMPL

