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
    mmx::point<double> w01= B - m_pt;
    mmx::point<double> u = m_dir;
    mmx::point<double> w1;
    double a= u.dot(u);
    double b= u.dot(v);
    double c= v.dot(v);
    double d= u.dot(w0);
    double e= v.dot(w0);
    double d1=u.dot(w01);
    double s0,t0,s1,t1,s,t,p0;
    s0 =(b*d-a*e)/(a*c-b*b);
    t0 =(c*d-b*e)/(a*c-b*b);
    p0=-e/c;
    mmx::point<double> w = w0 + s*v-t*u;
    mmx::point<double> H(A[0]+v[0]*p0, A[1]+v[1]*p0,A[2]+v[2]*p0);
    mmx::point<double> P = H-m_pt;

    using std::min;
    using std::max;
    if (a*c==b*b)
    {

        if (d>=0 || d1>=0)
        {
            s1=0; t1=e/b;
            w1=w0+s1*v-t1*u;
            return sqrt(w1.dot(w1));
        }
        else if(d<0 && d1<0)
        {   s1=0;t1=0;
            s=1;
            w=w0+s*v-t1*u;
            w1=w0+s1*v-t1*u;
            return min(sqrt(w1.dot(w1)),sqrt(w.dot(w)));
        }

    }
    if (a*c!=b*b)
    {
        if (d>=0 || d1>=0)
        {
            if (s0>=0 && s0<=1 && t0>=0)
            {
                s1=s0;t1=t0;
                w1=w0+s1*v-t1*u;
                return sqrt(w1.dot(w1));
            }

            if (s0<0 && t0>=0)
            {
                s1=0;t1=d/a;
                if(t1>=0)
                {
                    w1=w0+s1*v-t1*u;
                    return sqrt(w1.dot(w1));

                }
                else if(t1<0)
                {
                    s1=0;t1=0;
                    w1=w0+s1*v-t1*u;
                    return sqrt(w1.dot(w1));
                }

            }

            if (s0>1 && t0>=0)
            {
                s1=1;t1=(b+d)/a;
                if(t1>=0)
                {
                    w1=w0+s1*v-t1*u;
                    return sqrt(w1.dot(w1));

                }
                else if(t1<0)
                {
                    s1=1;t1=0;
                    w1=w0+s1*v-t1*u;
                    return sqrt(w1.dot(w1));
                }

            }



            if (s0>=0 && s0<=1 && t0<0)
            {
                s1=-e/c;t1=0;
                if(s1>=0 && s1<=1)
                {
                    w1=w0+s1*v-t1*u;
                    return sqrt(w1.dot(w1));
                }
                else if(s1<0)
                {
                    s1=0; t1=0;
                    w1=w0+s1*v-t1*u;
                    return sqrt(w1.dot(w1));
                }
                else if (s1>1)
                {
                    s1=1; t1=0;
                    w1=w0+s1*v-t1*u;
                    return sqrt(w1.dot(w1));

                }
            }
        }

        else if (d<0 && d1<0)
        {
            if (H[0]>=min(A[0],B[0]) && H[0]<=max(A[0],B[0]) && H[1]>=min(A[1],B[1]) && H[1]<=max(A[1],B[1]) && H[2]>=min(A[2],B[2]) && H[3]<=max(A[3],B[3]))
            {
                return  sqrt(P.dot(P));

            }
            else
                return min(sqrt(w0.dot(w0)), sqrt(w01.dot(w01)));



        }


    }
}








