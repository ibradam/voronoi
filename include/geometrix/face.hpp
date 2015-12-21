/**********************************************************************
 * PACKAGE  : geometrix 
 * COPYRIGHT: (C) 2015, Bernard Mourrain, Inria
 **********************************************************************/
#pragma once

# include <algorithm>
# include <vector>
# include <geometrix/point.hpp> 
# include <geometrix/edge.hpp> 

//====================================================================
namespace mmx {
//--------------------------------------------------------------------

/** @brief face of a mesh.
 *
 * A face is a list of indices of points of a mesh.
 **/
class face 
{
public:
  
  typedef std::vector<int>::const_iterator const_iterator;

  face(void) : m_id(-1) { };
  face(int p1, int p2, int p3): m_id(-1) {
    m_points.push_back(p1);
    m_points.push_back(p2);
    m_points.push_back(p3);
  }; 

  const_iterator begin() const {return m_points.begin();}
  const_iterator end()   const {return m_points.end();}

  unsigned       size()  const {return m_points.size();}
  std::vector<int>& points()      {return m_points;}

  int            vertex(int i)       {return m_points[i];}
  int            vertex(int i) const {return m_points[i];}

  int            operator[](int i)       {return m_points[i];}
  const int      operator[](int i) const {return m_points[i];}

  void           insert     (int p) { m_points.push_back(p); }
  void           push_vertex(int p) { m_points.push_back(p); }
  void           push_back_vertex(int p) { m_points.push_back(p); }
 
  void           set_id(const int & i) { m_id=i; }
  int            id(void)        const { return m_id; }

  void           reverse() { std::reverse(m_points.begin(), m_points.end()); } 

private:
  std::vector<int> m_points;
  int m_id;
};


//--------------------------------------------------------------------
} /* namespace mmx */
//====================================================================
# undef TMPL
# undef TMPL1
# undef SELF
# undef VIEWER
