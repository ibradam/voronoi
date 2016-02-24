/**********************************************************************
 * PACKAGE  : geometrix 
 * COPYRIGHT: (C) 2015, Bernard Mourrain, Inria
 **********************************************************************/
#pragma once
 
#define TMPL template <class C>

//====================================================================
namespace mmx {
//--------------------------------------------------------------------

template <int N, class C> struct tmsh_vertex;

TMPL
struct tmsh_vertex<2,C> {

    tmsh_vertex(void): m_id(-1), m_tag(0) {
        m_coord[0]=0; m_coord[1]=0; m_coord[2]=0;
        for(unsigned i=0;i<4;i++) m_neighbor[i]=-1;
    }

    tmsh_vertex(const C& X, const C& Y): m_id(-1), m_tag(0) {
        m_coord[0]=X; m_coord[1]=Y; m_coord[2]=0;
        for(unsigned i=0;i<4;i++) m_neighbor[i]=-1;
    }

    int        next (int v) const { return m_neighbor[2*v]; }
    int    previous (int v) const { return m_neighbor[2*v+1]; }
    int    neighbor (int& v) const;    

    C    operator[] (int v) const { return m_coord[v]; }
    C&   operator[] (int v)       { return m_coord[v]; }

    int  id(void)    const { return m_id; }
    int& id(void)          { return m_id; }

    int  tag(void)   const { return m_tag; }
    void tag(int r)        { m_tag = r; }

    bool has_no_tag (void) const { return m_tag==0; }

    static unsigned dim;

    C   m_coord[3];
    int m_neighbor[4];

    // The index of the point
    int m_id;
    // The tag of the point
    int m_tag;
};

TMPL unsigned tmsh_vertex<2,C>::dim = 2;

TMPL
struct tmsh_vertex<3,C> {

    tmsh_vertex(void): m_id(-1), m_tag(0) {
        m_coord[0]=0; m_coord[1]=0; m_coord[2]=0;
        for(unsigned i=0;i<6;i++) m_neighbor[i]=-1;
    };

    tmsh_vertex(const C& X, const C& Y, const C& Z): m_id(-1), m_tag(0) {
        m_coord[0]=X; m_coord[1]=Y; m_coord[2]=Z;
        for(unsigned i=0;i<6;i++) m_neighbor[i]=-1;
    };

    int        next (int v) const { return m_neighbor[2*v]; }
    int    previous (int v) const { return m_neighbor[2*v+1]; }

    C    operator[] (int v) const { return m_coord[v]; }
    C&   operator[] (int v)       { return m_coord[v]; }

    int  tag(void)   const { return m_tag; }
    void tag(int r)        { m_tag = r; }

    static unsigned dim;
    C   m_coord[3];
    int m_neighbor[6];

    // The index of the point
    int m_id;

    // The tag of the point
    int m_tag;

};

TMPL unsigned tmsh_vertex<3,C>::dim = 3;

template <class M, class P>
void set_middle(M& m, const P& a, const P& b) {
    for(unsigned i=0;i<M::dim;i++) m[i]=(a[i]+b[i])/2.;
}
//--------------------------------------------------------------------
} /* namespace mmx */
//====================================================================
#undef TMPL
