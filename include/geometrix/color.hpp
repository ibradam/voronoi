/**********************************************************************
 * PACKAGE  : geometrix 
 * COPYRIGHT: (C) 2015, Bernard Mourrain, Inria
 **********************************************************************/
#pragma once

#include <iostream>
#define SELF color
//====================================================================
namespace mmx {
//--------------------------------------------------------------------

struct color
{
  float r,g,b;
  color(): r(0.314),g(0.979),b(1){}
  color(double x, double y, double z): r(x),g(y),b(z){}
  color(int x, int y, int z): r(x/255.),g(y/255.),b(z/255.){}
  color(const SELF& c): r(c.r),g(c.g),b(c.b){}
  color& operator= (const SELF & c)
  {
    r=c.r; g=c.g; b=c.b; return *this;
  }
  bool operator ==(const color& c) const {return (this->r==c.r && this->g==c.g && this->b==c.b);}
  bool operator !=(const color& c) const {return !(*this==c); }

  unsigned red()   const {return r*255;}
  unsigned green() const {return g*255;}
  unsigned blue()  const {return b*255;}

};

inline std::ostream& operator<<(std::ostream& os, const SELF & c)
{
  os <<"Color("<<c.r<<","<<c.g<<","<<c.b<<")"; return os;
}

//--------------------------------------------------------------------
} /* namespace mmx */
//====================================================================
# undef SELF

