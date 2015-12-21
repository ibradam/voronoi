/**********************************************************************
 * PACKAGE  : geometrix
 * COPYRIGHT: (C) 2015, Bernard Mourrain, Inria
 **********************************************************************/
#pragma once

#include <geometrix/mdebug.hpp>
#include <geometrix/foreach.hpp>
#include <geometrix/point.hpp>
#include <geometrix/edge.hpp>
#include <geometrix/face.hpp>
#include <geometrix/color.hpp>
#include <map>

#define TMPL   template<class C,class V>
#define SELF   mesh<C,V>

//====================================================================
namespace mmx {
//--------------------------------------------------------------------

template<class C>
struct mesh_def {
    typedef point<C>     Point;
    typedef edge         Edge;
    typedef face         Face;
    typedef color        Color;
};

//--------------------------------------------------------------------
template<class C, class V=mesh_def<C> >
class mesh {
public:

    typedef typename V::Point    Point;
    typedef typename V::Edge     Edge;
    typedef typename V::Face     Face;
    typedef typename std::vector<Point*>::iterator          PointIterator;
    typedef typename std::vector<Point*>::const_iterator    PointConstIterator;
    typedef typename V::Color    Color;

    mesh();
    ~mesh(void) ;

    /// Number of vertices.
    unsigned nbv() const { return m_vertices.size(); }
    /// Number of edges.
    unsigned nbe() const { return m_edges.size(); }
    /// Number of faces.
    unsigned nbf() const { return m_faces.size(); }

    /// Number of colors.
    unsigned nbc() const { return m_colors.size(); }

    /// Remove all vertices, edges and faces.
    virtual void clear();


    inline std::vector<Point *> vertices(void) const { return m_vertices ; }
    inline std::vector<Edge *>  edges(void)    const { return m_edges ; }
    inline std::vector<Face *>  faces(void)    const { return m_faces ; }

    virtual       Point* vertex(int i)       { return m_vertices[i] ; }
    virtual const Point* vertex(int i) const { return m_vertices[i] ; }

    virtual       Edge*  edge  (int i)       { return m_edges[i] ; }
    virtual const Edge*  edge  (int i) const { return m_edges[i] ; }

    virtual       Face*  face  (int i)       { return m_faces[i] ; }
    virtual const Face*  face  (int i) const { return m_faces[i] ; }

    Color&             color(int i)          {return m_colors[i];}
    Color              color(int i)    const {return m_colors[i];}


    PointConstIterator begin() const { return this->m_vertices.begin() ; }
    PointIterator      begin()       { return this->m_vertices.begin() ; }
    PointConstIterator end()   const { return this->m_vertices.end(); }
    PointIterator      end()         { return this->m_vertices.end(); }

    void add_vertex(Point*);
    void add_vertex(const C& x, const C& y, const C& z) {
        this->add_vertex(new Point(x,y,z));
    }

    void push_back_vertex(Point* p) { this->add_vertex(p); }
    void push_back_vertex(const C& x, const C& y, const C& z) {
        this->add_vertex(new Point(x,y,z));
    }

    template<class POINT>
    void push_back_vertex(const POINT& p) {
        this->add_vertex(new Point(p.x(),p.y(),p.z()));
    }

    void push_back_color(const Color& c)        { m_colors.push_back(c); }
    void add_color      (const Color& c)        { m_colors.push_back(c); }
    void set_color      (int i, const Color& c) { m_colors[i]=c; }
    Color get_color(int i) const                { return m_colors[i]; }

    int index_of(Point* p) const    ;
    void set_index(Point* p, int i) ;

    virtual void add_edge(Point*, Point*);
    virtual void add_edge(Edge *);
    virtual void push_back_edge(int b, int e) {
        this->add_edge(b,e);
    }
    virtual void add_edge(int b, int e) {
        m_edges.push_back(new Edge(b,e));
    }

    virtual void add_face(Face *);
    virtual void add_face(Point*, Point*, Point*);
    virtual void add_face(Point*, Point*, Point*, Point*);

    virtual void push_back_face(int i0, int i1, int i2) {
        this->add_face(this->vertex(i0), this->vertex(i1), this->vertex(i2));
    }

    virtual void add_face(int i0, int i1, int i2) {
        Face* f= new Face;
        f->insert(i0);
        f->insert(i1);
        f->insert(i2);
        m_faces.push_back(f);
        //this->add_face(this->vertex(i0), this->vertex(i1), this->vertex(i2));
    }


    //private:

    std::vector<Point *>          m_vertices ;
    std::vector<Edge  *>          m_edges ;
    std::vector<Face  *>          m_faces ;
    std::vector<Color>            m_colors ;
    std::map<Point*,int>          m_indices ;
    int                           m_nbv ;
};
//--------------------------------------------------------------------
TMPL SELF::mesh(void) : m_nbv(0) {
}

TMPL SELF::~mesh(void) {
}

TMPL void
SELF::add_vertex(Point *p) {

    if (this->index_of(p)<0) {
        this->set_index(p,m_nbv);
        m_vertices.push_back(p);
        m_nbv++;
    }

}

TMPL int
SELF::index_of(Point * p) const {
    //    return p->index();
    typename std::map<Point*,int>::const_iterator it = m_indices.find( p );
    if ( it == m_indices.end() ) {
        return -1;
    }
    else {
        return it->second;
    }
}

TMPL void
SELF::set_index(Point * p, int i) {
    //    p->set_index(i); return;
    m_indices[p]=i;
}


TMPL void
SELF::add_edge(Point * p1, Point * p2) {
    //std::cout<<"Insert edge"<<std::endl;
    if (this->index_of(p1)<0) this->add_vertex(p1);
    if (this->index_of(p2)<0) this->add_vertex(p2);
    m_edges.push_back(new Edge(this->index_of(p1),this->index_of(p2)));
}


TMPL void
SELF::add_edge(Edge* e) {
    m_edges.push_back(e);
}

TMPL void
SELF::add_face(Face * f) {
    //std::cout<< " add face "<< f->size() << std::endl;
    m_faces.push_back(f);
}

TMPL void
SELF::add_face(Point * p1, Point* p2, Point* p3) {
    Face* f= new Face;
    if (this->index_of(p1)<0) this->add_vertex(p1);
    if (this->index_of(p2)<0) this->add_vertex(p2);
    if (this->index_of(p3)<0) this->add_vertex(p3);
    f->insert(this->index_of(p1));
    f->insert(this->index_of(p2));
    f->insert(this->index_of(p3));
    m_faces.push_back(f);
}

TMPL void
SELF::add_face(Point * p1, Point* p2, Point* p3, Point* p4) {
    Face* f= new Face;
    if (this->index_of(p1)<0) this->add_vertex(p1);
    if (this->index_of(p2)<0) this->add_vertex(p2);
    if (this->index_of(p3)<0) this->add_vertex(p3);
    if (this->index_of(p4)<0) this->add_vertex(p4);
    f->insert(this->index_of(p1));
    f->insert(this->index_of(p2));
    f->insert(this->index_of(p3));
    f->insert(this->index_of(p4));
    m_faces.push_back(f);
}

TMPL void
SELF::clear() {
    // Should we clean up these points?
    m_faces.resize( 0 );
    m_edges.resize( 0 );
    m_vertices.resize( 0 );
}

template<class GRAPHIC, class BOUNDINGBOX>
void insert_bounding_box2d(GRAPHIC* g, BOUNDINGBOX *bx, bool cross = false) {
    typedef typename GRAPHIC::Point Point;
    Point
            *p0= new Point(bx->xmin(),bx->ymin()),
            *p1= new Point(bx->xmin(),bx->ymax()),
            *p2= new Point(bx->xmax(),bx->ymax()),
            *p3= new Point(bx->xmax(),bx->ymin());

    g->add_vertex(p0);g->add_vertex(p1); g->add_edge(p0,p1);
    g->add_vertex(p1);g->add_vertex(p2); g->add_edge(p1,p2);
    g->add_vertex(p2);g->add_vertex(p3); g->add_edge(p2,p3);
    g->add_vertex(p3);g->add_vertex(p0); g->add_edge(p3,p0);

    if(cross) {
        g->add_edge(p0,p2);
        g->add_edge(p1,p3);
    }
}

template<class GRAPHIC, class BOUNDINGBOX>
void insert_boundingbox(GRAPHIC* g, BOUNDINGBOX *bx, bool cross = false) {
    typedef typename GRAPHIC::Point Point;
    Point
            *p0= new Point(bx->xmin(),bx->ymin(),bx->zmin()),
            *p1= new Point(bx->xmin(),bx->ymax(),bx->zmin()),
            *p2= new Point(bx->xmax(),bx->ymax(),bx->zmin()),
            *p3= new Point(bx->xmax(),bx->ymin(),bx->zmin());

    g->add_vertex(p0);g->add_vertex(p1); g->add_edge(p0,p1);
    g->add_vertex(p1);g->add_vertex(p2); g->add_edge(p1,p2);
    g->add_vertex(p2);g->add_vertex(p3); g->add_edge(p2,p3);
    g->add_vertex(p3);g->add_vertex(p0); g->add_edge(p3,p0);

    Point
            *q0= new Point(bx->xmin(),bx->ymin(),bx->zmax()),
            *q1= new Point(bx->xmin(),bx->ymax(),bx->zmax()),
            *q2= new Point(bx->xmax(),bx->ymax(),bx->zmax()),
            *q3= new Point(bx->xmax(),bx->ymin(),bx->zmax());

    g->add_vertex(q0);g->add_vertex(q1); g->add_edge(q0,q1);
    g->add_vertex(q1);g->add_vertex(q2); g->add_edge(q1,q2);
    g->add_vertex(q2);g->add_vertex(q3); g->add_edge(q2,q3);
    g->add_vertex(q3);g->add_vertex(q0); g->add_edge(q3,q0);

    g->add_vertex(p0);g->add_vertex(q0);g->add_edge(p0,q0);
    g->add_vertex(p1);g->add_vertex(q1);g->add_edge(p1,q1);
    g->add_vertex(p2);g->add_vertex(q2);g->add_edge(p2,q2);
    g->add_vertex(p3);g->add_vertex(q3);g->add_edge(p3,q3);

    if(cross) {
        Point
                *r0= new Point(bx->xmin(),bx->ymin(),bx->zmax()),
                *r1= new Point(bx->xmin(),bx->ymax(),bx->zmax()),
                *r2= new Point(bx->xmax(),bx->ymax(),bx->zmax()),
                *r3= new Point(bx->xmax(),bx->ymin(),bx->zmax());

        g->add_vertex(r0);g->add_vertex(r2);g->add_edge(r0,r2);
        g->add_vertex(r1);g->add_vertex(r3);g->add_edge(r1,r3);
    }
}


//--------------------------------------------------------------------
template<class VIEWER, class C, class V, class COLOR>
VIEWER& axl_print(VIEWER& out, const SELF& tp, COLOR c= COLOR()) {

    typedef typename SELF::Point  Point;
    typedef typename SELF::Edge   Edge;
    typedef typename SELF::Face   Face;
    typedef COLOR  Color;
    //mdebug()<<"print mesh"<<tp.nbv()<<tp.nbe()<<tp.nbf();
    out<<"<mesh";
    if(c != Color()) out<< " color=\""<<(unsigned)(c.r*255)<<" "<<(unsigned)(c.g*255)<<" "<<(unsigned)(c.b*255)<<"\"";
    out<<     " size=\"0.1\" >\n";
    out<<"<count>"<<tp.nbv()<<" "<< tp.nbe()<<" "<< tp.nbf() <<"</count>\n";

    if (tp.nbv()>0){
        out<<"<points>\n";
        int c=0;
        foreach(Point* p, tp.vertices()) {
            out <<" "<<p->x()<<" "<<p->y()<<" "<<p->z()<<"\n";
            //mdebug()<<c<< p->x()<<p->y()<<p->z(); c++;
        }
        out<<"</points>\n";
    }
    //mdebug()<<"Points";

    if (tp.nbc()==tp.nbv()){
        out<<"<colors>\n";
        for(unsigned i=0; i< tp.nbc(); i++) {
            out <<" "<<tp.color(i).red()<<" "<<tp.color(i).green()<<" "<<tp.color(i).blue() <<"\n";
        }
        out<<"</colors>\n";
    }
    //std::cout<<"Colors"<<std::endl;

    if (tp.nbe()>0){
        out<<"<edges>\n";
        foreach(Edge* e, tp.edges()) {
            out <<"2 "<< e->source()<<" "<< e->destination() <<"\n";
        }
        out<<"</edges>\n";
    }
    //std::cout<<"Edges"<<std::endl;

    if(tp.nbf()>0) {
        out<<"<faces>\n";
        foreach(Face* f, tp.faces()) {
            //      out<<"3 ";
            out<<f->size();
            foreach(int p, f->points())
            {
                out <<" "<<p;
            }
            out <<"\n";
        }
        out<<"</faces>\n";
    }
    //std::cout<<"Faces"<<std::endl;
    out<<"</mesh>\n";
    return out;
}
//====================================================================
} /* namespace mmx */
//====================================================================
# undef TMPL
# undef SELF

