#include <math.h>
#include <voronoi/mdebug.h>
#include <voronoi/hline.h>

/*!
 * \brief hline::distance2
 * \param p
 * \return the euclidean distance from the half-line to a point p.
 */
double hline::distance2(const mmx::point<double>& p) const {
    mmx::point<double> u = p - m_pt;
    double s = m_dir.dot(u);
    if (s<0)
        return m_pt.dist2(p);
    else
        return m_pt.dist2(p)-s*s;
}

double hline::distance2(double X, double Y, double Z) const {
    mdebug()<<"constructor"<<X<<Y<<Z;
    return this->distance2( mmx::point<double>(X,Y,Z) );
}

double hline::distance2(const mmx::point<double>& A, const mmx::point<double>& B) const {

    mmx::point<double> v = B - A;
    mmx::point<double> w0= A - m_pt;
    mmx::point<double> u = m_dir;
    mmx::point<double> w1;
    double a= u.dot(u);
    double b= u.dot(v);
    double c= v.dot(v);
    double d= u.dot(w0);
    double e= v.dot(w0);
    double s0,t0,s1,t1,s,t;
    s0 =(b*d-a*e)/(a*c-b*b);
    t0 =(c*d-b*e)/(a*c-b*b);
    mmx::point<double> w = w0 + s*v-t*u;
    if (a*c==b*b)
    {
        s1=0; t1=e/b;

        w1=w0+s1*v-t1*u;
        return sqrt(w1.dot(w1));
    }

  else if (s0>=0 && s0<=1 && t0>=0)
  {
  s1=s0;t1=t0;
  w1=w0+s1*v-t1*u;
  return sqrt(w1.dot(w1));
  }

}



