/**********************************************************************
 * PACKAGE  : geometrix 
 * COPYRIGHT: (C) 2015, Bernard Mourrain, Inria
 **********************************************************************/
#pragma once

#ifndef POINT_HPP
#define POINT_HPP

#include <cmath>

#define TMPL  template<class C, int N>
#define SELF  point<C,N>
#define Viewer viewer<axel,W>
#define Ostream std::ostream

#undef Scalar

//====================================================================
namespace mmx { 
//--------------------------------------------------------------------

template<class C, int N= 3>
class point 
{
    public:
    typedef C   Scalar;

    point(void) ;
    point(Scalar x, Scalar y=0, Scalar z= 0) ;
    point(const SELF & p) ;

    inline Scalar  x(void) const { return m_coord[0] ; }
    inline Scalar& x(void)       { return m_coord[0] ; }
    inline Scalar  y(void) const { if (N>1) return m_coord[1] ; else return 0; }
    inline Scalar& y(void)       { if (N>1) return m_coord[1] ; else return *new Scalar(0); }
    inline Scalar  z(void) const { if (N>2) return m_coord[2] ; else return 0; }
    inline Scalar& z(void)       { if (N>2) return m_coord[2] ; else return *new Scalar(0); }

    Scalar   operator [] (const int & i) const { return m_coord[i] ; }
    Scalar & operator [] (const int & i)       { return m_coord[i] ; }

    Scalar   at (const int & i) const { return m_coord[i] ; }
    Scalar & at (const int & i)       { return m_coord[i] ; }

    inline void setx(Scalar x) { this->m_coord[0] = x ; }
    inline void sety(Scalar y) { if (N>1) this->m_coord[1] = y ; }
    inline void setz(Scalar z) { if (N>2) this->m_coord[2] = z ; }

    bool     operator == (const point & other) const ;
    bool     operator != (const point & other) const ;
    bool     operator <  (const point & other) const ;
    SELF &   operator  = (const point & other) ;

    // Linear Algebra operations
    point operator+(const point &v)const; //Addition
    point operator-(const point &v)const; //Substraction
    point operator*(const Scalar&k)const;	//Multiplication by Scalar
    //Division by Scalar
    point operator/(const Scalar&k)const { return *this*( Scalar(1)/k); };	//Division by Scalar

    point operator*(const point &v)const; //Coordinate-wise Product

    point& operator+=(const point &v) //Addition
    { for(unsigned i=0;i<N;i++) this->at(i)+=v.at(i); return *this; }
    point& operator-=(const point &v) //Substraction
    { for(unsigned i=0;i<N;i++) this->at(i)-=v.at(i); return *this; }
    point& operator*=(const Scalar&k) //Multiplication by Scalar
    { for(unsigned i=0;i<N;i++) this->at(i)*=k; return *this; }
    point& operator/=(const Scalar&k) //Division by Scalar
    { for(unsigned i=0;i<N;i++) this->at(i)/=k; return *this; }

    point operator*=(const point &v)const //Coordinate-wise Product
    { for(unsigned i=0;i<N;i++) this->at(i)*=v.at(i); return *this; }

    point  operator-()const;		//Negative of a point
    point  cross(const point &v) const;   //Cross Product
    Scalar dot(point v2)const;		//Dot Product
    Scalar dist2(point v2)const;		//Distance Product
    Scalar norm() const;                  //L2 Norm

    private:
    Scalar m_coord[N] ;
    //Scalar m_coord[1] ;
    //Scalar m_coord[2] ;
};

TMPL
SELF::point(void)
{
    for(unsigned i=0;i<N;i++) this->m_coord[i] = (Scalar)0 ;
    //this->m_coord[1] = (Scalar)0 ;
    //this->m_coord[2] = (Scalar)0 ;
}

TMPL
SELF::point(Scalar x, Scalar y, Scalar z)
{
    this->m_coord[0] = x ;
    this->m_coord[1] = y ;
    this->m_coord[2] = z ;
}
TMPL
SELF::point(const SELF & other)
{
    for(unsigned i=0;i<N;i++)
        this->m_coord[i] = other[i] ;
}

TMPL
SELF & SELF::operator = (const SELF & other) 
{
    if(this == &other)
        return *this ;

    for(unsigned i=0;i<N;i++)
        this->m_coord[i] = other[i] ;

    return * this ;
}

TMPL bool 
SELF::operator == (const SELF & other) const {
    bool res = true;
    for(unsigned i=0;i<N && res; i++)
        res = (m_coord[i]==other[i]);
    return res;
}

TMPL bool 
SELF::operator != (const SELF & other) const {
    return !((*this)==other);
}


TMPL bool
SELF::operator < (const SELF & other) const {
    return ((this->x() < other.x())
            || ((this->x() == other.x()) && (this->y() < other.y()))
            || ((this->x() == other.x()) && (this->y() == other.y()) && (this->z() < other.z()))) ;
}

TMPL
SELF SELF::operator +(const SELF &v)const  //Addition
{
    SELF res;
    for(unsigned i=0;i<N;i++ )
        res[i] = m_coord[i] + v[i];
    return res;
}

TMPL
SELF SELF::operator -(const SELF &v)const  //Subtraction
{
    SELF res;
    for(unsigned i=0;i<N;i++ )
        res[i] = m_coord[i] - v[i];

    return res;
}

TMPL
SELF SELF::operator *(const Scalar &k)const  //Multiplication By Scalar
{
    SELF res;
    for(unsigned i=0;i<N;i++ )
        res[i] = m_coord[i] *k;
    return res;
}

TMPL
SELF operator *(typename SELF::Scalar k,const SELF &v) //Left multiplication By Scalar
{
    return v * k;
}

TMPL
SELF SELF::operator -()const   //Negative Of a point
{
    SELF res;
    for(unsigned i=0;i<N;i++ )
        res[i] = -m_coord[i];
    return res;
}

TMPL
SELF SELF::operator *(const SELF &v)const // Coordinate-wise Product
{
    SELF res;
    for(unsigned i=0;i<N;i++ )
        res[i] = m_coord[i] * v[i];
    return res;
}


TMPL
SELF SELF::cross(const SELF &v)const //Cross Product
{
    SELF res;
    res.x() = (this->y() * v.z()) - (this->z() * v.y());
    res.y() = (this->z() * v.x()) - (this->x() * v.z());
    res.z() = (this->x() * v.y()) - (this->y() * v.x());
    return res;
}

TMPL
SELF cross(const SELF v1, const SELF v2) //Cross Product of 2 points
{
    return v1.cross(v2);
}

TMPL
typename SELF::Scalar SELF::dot(const SELF v2) const //Dot product
{
    Scalar res=0;
    for(unsigned i=0;i<N;i++) res += m_coord[i] * v2[i];
    return res;
}

TMPL
typename SELF::Scalar SELF::dist2(const SELF v2) const //Dot product
{
    Scalar res=0;
    for(unsigned i=0;i<N;i++) res+=(m_coord[i] - v2[i])*(m_coord[i] - v2[i]);
            return res;
}
TMPL
typename SELF::Scalar dot(const SELF v1, const SELF v2) //Dot product of 2 points
{
    return v1.dot(v2);
}

TMPL
typename SELF::Scalar SELF::norm() const   //Norm
{
    return std::sqrt( this->dot(*this) );
}

TMPL inline typename SELF::Scalar read (const SELF& v, unsigned i) {
    return v[i];
}

TMPL inline typename SELF::Scalar distance (const SELF& p1, const SELF& p2) {
    return std::sqrt(p1.dist2(p2));
}

//--------------------------------------------------------------------
TMPL Ostream&
operator<<(Ostream& os, const SELF& p) {
    os <<(char *)"[";
    for(unsigned i=0;i<N;i++){
        if(i>0) os <<(char *)" ";
        os <<p[i];
    }
    os <<(char *)"]";
    return os;
}
//--------------------------------------------------------------------
template<class VIEWER, class C, class V, int N>
VIEWER&
operator<<(VIEWER& os, const SELF& p) {
    os<<"<point color=\""<<(int)(255*os.color.r)<<" "<<(int)(255*os.color.g)<<" "<<(int)(255*os.color.b)<<"\">"
     <<p.x()<<" "<<p.y()<<" "<<p.z()
    <<"</point>\n";
    return os;
}

//--------------------------------------------------------------------
} /* namespace mmx */
//====================================================================
#undef TMPL
#undef SELF
#undef Scalar
#undef Viewer
#undef Ostream

#endif //POINT_HPP
