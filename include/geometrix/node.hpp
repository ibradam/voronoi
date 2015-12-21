/**********************************************************************
 * PACKAGE  : geometrix 
 * COPYRIGHT: (C) 2015, Bernard Mourrain, Inria
 **********************************************************************/
#pragma once

#include <iostream>
#include <vector>

#define NODE node<CELL>
#define TMPL template<class CELL>

 //====================================================================
namespace mmx {
//--------------------------------------------------------------------

template <class CELL>
struct node {
public:
    enum NODE_TYPE { LEFT, RIGHT } ;

    typedef CELL Cell;

public:
    node(void) ;
    node(CELL* cl) ;
    node(NODE* parent, CELL* cl, NODE_TYPE nodeType, int v=0) ;
    node(NODE* left, NODE * right, CELL cl, int v=0) ;

protected:
    node(NODE& node) ;

public:
    inline void set_cell(CELL* c) { m_cell = c ; }
    inline void set_parent(NODE * n)     { m_parent = n ; }

    inline void set_leftchild(NODE * n)  { m_left = n ; }
    inline void set_rightchild(NODE * n) { m_right = n ; }

    inline const CELL* get_cell(void) const { return m_cell; }
    inline       CELL* get_cell(void)       { return m_cell ; }

    //inline Object object(void) { return m_objects.front() ; }

    inline NODE_TYPE type(void) const { return m_type ; }
    inline int split_dir(void) const { return m_var; }

    inline NODE * left  (void) { return m_left ; }
    inline NODE * right (void) { return m_right ; }
    inline NODE * parent(void) { return m_parent ; }

    inline const NODE * left  (void) const { return m_left ; }
    inline const NODE * right (void) const { return m_right ; }
    inline const NODE * parent(void) const { return m_parent ; }

    int    var(void) const {return this->m_var;}
    bool   is_leaf(void) const;
    size_t leaf_distance() const;

public:
    CELL*               m_cell ;
    NODE_TYPE           m_type ;
    int                 m_var ;

    NODE * m_parent  ;
    NODE * m_left ;
    NODE * m_right ;

    int depth ;
    int index ;

} ;

//--------------------------------------------------------------------
TMPL NODE::node(void)
{
    m_type   = LEFT ;
    m_var    = 0 ;

    m_parent = NULL ;
    m_left   = NULL;
    m_right  = NULL;

    m_cell   = NULL ;

    depth = 0;
    index = 0;
}

TMPL NODE::node(CELL* cl)
{
    this->m_type   = LEFT ;
    this->m_var    = 0 ;

    this->m_parent = NULL ;
    this->m_left   = NULL;
    this->m_right  = NULL;

    this->m_cell   = cl ;

    depth = 0;
    index = 0;

}

TMPL NODE::node(NODE * left, NODE * right, CELL cl, int v) {

    this->m_parent = NULL ;
    this->m_cell   = cl ;
    this->m_var    = v ;

    left->type  = LEFT ;  left->parent = this ; m_left = left;
    right->type = RIGHT; right->parent = this ; m_right= right;

    depth = left->depth-1 ;

}

TMPL NODE::node(NODE * parent, CELL* cl, NODE_TYPE type, int v)
{
    this->m_cell   = cl ;
    this->m_type   = type ;
    this->m_var    = v ;
    this->m_parent = parent ;

    this->m_left   = NULL;
    this->m_right  = NULL;

    depth = parent->depth+1 ;

    switch(type) {
    case LEFT : parent->set_leftchild(this) ; break ;
    case RIGHT: parent->set_rightchild(this) ; break ;
    default: std::cerr << "Error : the node's type isn't appropriate \n" ; break ;
    }
}

TMPL bool NODE::is_leaf(void) const
{
    if((m_left == NULL) && (m_right == NULL) )
        return true;
    return false;
}

TMPL size_t NODE::leaf_distance(void) const
{
    if ( this->is_leaf() )
        return 0;

    struct inner {
        size_t operator()( const NODE* node ) {
            if ( node == 0 )
                return 0;
            else
                return node->leaf_distance();
        }
    } I;

    size_t d = 0;
    d = std::min( d, I(m_left)  );
    d = std::min( d, I(m_right) );

    return d+1;
}


//--------------------------------------------------------------------
} /* namespace mmx */
//====================================================================
#undef NODE
#undef TMPL

