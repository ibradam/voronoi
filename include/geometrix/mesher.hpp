/**********************************************************************
 * PACKAGE  : geometrix 
 * COPYRIGHT: (C) 2015, Bernard Mourrain, Inria
 **********************************************************************/
#pragma once

#include <fstream>
#include <vector>
#include <geometrix/mesh.hpp>
#include <geometrix/subdivision.hpp>

#define TMPL template<class Variant>
#define SELF mesher<Variant>

//====================================================================
namespace mmx {
//--------------------------------------------------------------------

/** @brief mesher of shapes.
 *  
 **/
template<class Variant>
struct mesher {

    typedef typename Variant::Controler          Controler;
    typedef typename Controler::Cell             Input;
    typedef subdivision<Controler>               Subdivisor;
    typedef typename Variant::Polygonizer        Polygonizer;
 
    typedef typename Variant::Mesh               Output;
    typedef typename Output::Point               Point;

    mesher(double e1= 0.5, double e2=0.1);
    mesher(Controler* dtr, double e1= 0.5, double e2=0.1);
    ~mesher(void) ;

    void   set_input(Input* c);
    void   set_max_size(double epsilon) { m_max_size = epsilon; }
    void   set_min_size(double epsilon) { m_min_size = epsilon; }

    Input*     input()               { return m_input; }
    Output*    output()              { return m_output; }

    Controler* controler()           { return d; }

    void run (void);

private:

    Controler*  d;
    double      m_max_size, m_min_size;
    Input*      m_input;
    Output*     m_output;

};

//--------------------------------------------------------------------
TMPL SELF::mesher(double e1, double e2): m_max_size(e1), m_min_size(e2)
{
    d = new Controler;
    m_input = NULL; 
    m_output = new Output;
}
//--------------------------------------------------------------------
TMPL SELF::mesher(Controler* f, double e1, double e2): m_max_size(e1), m_min_size(e2)
{
    d = f;
    m_input = NULL; //f->m_root;
    m_output = new Output;
}
//--------------------------------------------------------------------
TMPL SELF::~mesher(void) {

}

TMPL void SELF::set_input(Input* cl){
    //mdebug()<<"set_input";
     m_input = cl;
}

//--------------------------------------------------------------------
TMPL
void SELF::run(void) {

    if (m_input==NULL) return;

    Subdivisor s(d,m_max_size, m_min_size);
    s.set_input(this->input());
    s.run();
    mdebug()<<"Subdivision"<<this->controler()->m_regular.size() <<this->controler()->m_singular.size();


    Polygonizer p(this->controler());
    p.set_regular(&this->controler()->m_regular);
    p.set_singular(&this->controler()->m_singular);
    p.run();
    //mdebug()<<"Polygonizer done";

    //m_output=p.output();
    p.get_output(m_output);

    mdebug()<<"TMesh:"<<m_output->nbv()<<m_output->nbe()<<m_output->nbf();
}

//--------------------------------------------------------------------
} /* namespace mmx */
//====================================================================
#undef TMPL
#undef SELF
