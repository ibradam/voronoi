/**********************************************************************
 * PACKAGE  : geometrix 
 * COPYRIGHT: (C) 2015, Bernard Mourrain, Inria
 **********************************************************************/
#pragma once

# include <geometrix/point.hpp>

# define TMPL  template<class POINT>
# define SELF  edge<POINT>
//====================================================================
namespace mmx { 
//--------------------------------------------------------------------

struct edge 
{
public:
  

  edge(void) {}
  edge(int source, int destination) {
    this->m_source = source ;
    this->m_destination = destination ;
   }

  inline int  source(void)            { return m_source ; }
  inline int  source(void)      const { return m_source ; }
  inline int  destination(void)       { return m_destination ; }
  inline int  destination(void) const { return m_destination ; }
  inline void reverse() 
    { int t=m_source; m_source=m_destination; m_destination=t; }
  

private:
  int m_source, m_destination ;
};

} /* namespace mmx */
//====================================================================
# undef TMPL
# undef SELF  

