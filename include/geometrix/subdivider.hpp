/**********************************************************************
 * PACKAGE  : geometrix
 * COPYRIGHT: (C) 2015, Bernard Mourrain, Inria
 **********************************************************************/
#pragma once

#include <list>
#include <geometrix/regularity.hpp>
#include <geometrix/node.hpp>


#define TMPL  template<class CONTROLER>
#define SELF  subdivider<CONTROLER>

//====================================================================
namespace mmx {
//--------------------------------------------------------------------

template<class CONTROLER >
class subdivider {
public:

    typedef CONTROLER                 Controler;
    typedef typename CONTROLER::Cell  Cell;

    typedef node<Cell>                Node;

    typedef Cell                      Input;

    subdivider(double e1=0.5, double e2=0.05);
    subdivider(Controler* ctrl, double e1=0.5, double e2=0.05);

    ~subdivider(void) ;

    void   set_input (Cell* bx);

    void   set_max_size(double epsilon) { m_max_size = epsilon; }
    void   set_min_size(double epsilon) { m_min_size = epsilon; }

    double max_size() { return m_max_size; }
    double min_size() { return m_min_size; }

    //Node*   input  (void)     { return m_root; }
    Controler* output (void)     { return this; }

    unsigned subdivide (Node * node);

    void run(void);

    Controler* controler() { return d; }

    /* Remove all vertices, edges and faces from the subdivision. */
    virtual void clear();

    template<class MESH> void get_tmsh(MESH* m);

private:
    Controler*          d;
    std::list<Node *>   m_nodes;
    Node*               m_root;
    double              m_max_size, m_min_size;
    //Output*             m_output;

};
//--------------------------------------------------------------------
TMPL SELF::subdivider(double e1, double e2):
    d(new Controler), m_max_size(e1), m_min_size(e2)
{
    //    m_output = new Output;
}
//--------------------------------------------------------------------
TMPL SELF::subdivider(Controler* ctrl, double e1, double e2):
    d(ctrl), m_max_size(e1), m_min_size(e2)
{ 
    //    m_output = new Output;
}
//--------------------------------------------------------------------
TMPL SELF::~subdivider(void) {
    //    delete m_output;
}
//--------------------------------------------------------------------
TMPL void SELF::set_input(Cell* c0)
{ 
    m_root = new Node;
    m_root->set_cell(c0);
    m_nodes.push_back(m_root) ;
}

//--------------------------------------------------------------------
TMPL void SELF::clear() {

}
//--------------------------------------------------------------------
TMPL unsigned SELF::subdivide(Node * node)
{
    Cell * cl = node->get_cell(), * left=0, * right=0;
    int v = d->split_direction(*cl);

    d->subdivide(v, cl, left, right);

    node->m_left = new Node(node, left,  Node::LEFT , v) ;
    node->m_right= new Node(node, right, Node::RIGHT, v) ;

    m_nodes.push_back(node->m_left);
    m_nodes.push_back(node->m_right);

    return v ;
}

//--------------------------------------------------------------------
TMPL void SELF::run() {

    double maxsz = m_max_size;//*this->output()->root()->get_cell()->size();
    double minsz = m_min_size;//*this->output()->root()->get_cell()->size();

    //mdebug()<<"Max size: "<< maxsz<< " Min size: "<<minsz;
    int r;
    while(!m_nodes.empty()) {
        Node* node = m_nodes.front() ;
        m_nodes.pop_front();
        mdebug()<<">> nodes"<<m_nodes.size();
        Cell* cl = node->get_cell() ;
        if((r = d->regularity(cl)) != OUTSIDE) {
            //mdebug()<<">> regularity"<<r<<d->size(cl)<<maxsz<<minsz;
            if( !ISINSIDE(r) && d->size(cl) > maxsz )
            {
                //mdebug()<<" subdiv > max, r:"<<r;
                this->subdivide(node);
            }
            else if( ISBOUNDARY(r) || ISINSIDE(r) )
            {
                mdebug()<<"process_regular";
                d->process_regular(cl);
            }
            else if(d->size(cl) > minsz)
            {
                //std::cout<<"Subdivide not regular cell "<<cl->size()<<std::endl;
                //mdebug()<<" subdiv min"<<r;
                this->subdivide(node);
            }
            else {
                mdebug()<<"process_singular"<<r;
                d->process_singular(cl);
            }
        }
    }
}

//--------------------------------------------------------------------
TMPL
template <class MESH>
void SELF::get_tmsh(MESH* m) {

    for(unsigned i=0;i< this->controler()->nbv();i++) {
        m->add_vertex(this->controler()->vertex(i)[0], this->controler()->vertex(i)[1], this->controler()->vertex(i)[2]);
    }

    for(unsigned i=0;i<this->controler()->nbv();i++) {
        for(unsigned j=0; j< SELF::Cell::dim;j++)
            if(this->controler()->vertex(i).m_neighbor[2*j]>=0)
                m->add_edge(i,this->controler()->vertex(i).m_neighbor[2*j]);
    }

}

//--------------------------------------------------------------------
} // namespace mmx
//====================================================================
#undef TMPL1
#undef TMPL
#undef SELF

