#include <math.h>
#include <voronoi/hline.h>
#include <voronoi/point.hpp>
# include <voronoi/mdebug.h>
/*!
 * \brief equidist
 * \param H1
 * \param H2
 * \param A
 * \param B
 * \return the point equidistant to H1 and H2 on the segment [A,B]
 */
mmx::point<double> equidist(const hline& H1 , const hline& H2,
                            const mmx::point<double>& A,
                            const mmx::point<double>& B)
{
    mmx::point<double> p;
    double a,b,t3,t2,x1,x2,y1,y2,z1,z2,xa,xb,ya,yb,za,zb;
    x1=H1.m_pt[0]; y1=H1.m_pt[1]; z1=H1.m_pt[2];
    x2=H2.m_pt[0]; y2=H2.m_pt[1]; z2=H2.m_pt[2];
    xa=A[0]; ya=A[1]; za=A[2];
    xb=B[0]; yb=B[1]; zb=B[2];

    double t1 = 0.5*((x1*x1-2*x1*xb - x2*x2 + 2*x2*xb + y1*y1 - 2*y1*yb - y2*y2 + 2*y2*yb + z1*z1 - 2*z1*zb - z2*z2 + 2*z2*zb)/(x1*xa - x1*xb - x2*xa + x2*xb + y1*ya - y1*yb - y2*ya + y2*yb + z1*za - z1*zb - z2*za + z2*zb));
    double t21 = (x1 * xa - x1 * xb - x2 * xa + x2 * xb + y1 * ya - y1 * yb - y2 * ya + y2 * yb + z2 * za -  z2 * zb - za * zb + zb * zb + sqrt(x1 * x1 * xa * xa - 2 * x1 * x1 * xa * xb + x1 * x1 * xb * xb - x1 * x1 * za * za + 2 * x1 * x1 * za * zb - x1 * x1 * zb * zb - 2 * x1 * x2 * xa * xa + 4 * x1 * x2 * xa * xb - 2 * x1 * x2 * xb * xb + 2 * x1 * xa * y1 * ya - 2 * x1 * xa * y1 * yb - 2 * x1 * xa * y2 * ya + 2 * x1 * xa * y2 * yb + 2 * x1 * xa * z2 * za - 2 * x1 * xa * z2 * zb - 2 * x1 * xa * za * zb + 2 * x1 * xa * zb * zb - 2 * x1 * xb * y1 * ya + 2 * x1 * xb * y1 * yb + 2 * x1 * xb * y2 * ya - 2 * x1 * xb * y2 * yb - 2 * x1 * xb * z2 * za + 2 * x1 * xb * z2 * zb + 2 * x1 * xb * za * za - 2 * x1 * xb * za * zb + x2 * x2 * xa * xa - 2 * x2 * x2 * xa * xb + x2 * x2 * xb * xb + x2 * x2 * za * za - 2 * x2 * x2 * za * zb + x2 * x2 * zb * zb - 2 * x2 * xa * y1 * ya + 2 * x2 * xa * y1 * yb + 2 * x2 * xa * y2 * ya - 2 * x2 * xa * y2 * yb - 2 * x2 * xa * z2 * za + 2 * x2 * xa * z2 * zb + 2 * x2 * xa * za * zb - 2 * x2 * xa * zb * zb + 2 * x2 * xb * y1 * ya - 2 * x2 * xb * y1 * yb - 2 * x2 * xb * y2 * ya + 2 * x2 * xb * y2 * yb + 2 * x2 * xb * z2 * za - 2 * x2 * xb * z2 * zb - 2 * x2 * xb * za * za + 2 * x2 * xb * za * zb + y1 * y1 * ya * ya - 2 * y1 * y1 * ya * yb + y1 * y1 * yb * yb - y1 * y1 * za * za + 2 * y1 * y1 * za * zb - y1 * y1 * zb * zb - 2 * y1 * y2 * ya * ya + 4 * y1 * y2 * ya * yb - 2 * y1 * y2 * yb * yb + 2 * y1 * ya * z2 * za - 2 * y1 * ya * z2 * zb - 2 * y1 * ya * za * zb + 2 * y1 * ya * zb * zb - 2 * y1 * yb * z2 * za + 2 * y1 * yb * z2 * zb + 2 * y1 * yb * za * za - 2 * y1 * yb * za * zb + y2 * y2 * ya * ya - 2 * y2 * y2 * ya * yb + y2 * y2 * yb * yb + y2 * y2 * za * za - 2 * y2 * y2 * za * zb + y2 * y2 * zb * zb - 2 * y2 * ya * z2 * za + 2 * y2 * ya * z2 * zb + 2 * y2 * ya * za * zb - 2 * y2 * ya * zb * zb + 2 * y2 * yb * z2 * za - 2 * y2 * yb * z2 * zb - 2 * y2 * yb * za * za + 2 * y2 * yb * za * zb))/(za*za -  2*za*zb + zb*zb);
    double t22 = -(- x1 * xa + x1 * xb + x2 * xa - x2 * xb - y1 * ya + y1 * yb+ y2 * ya - y2 * yb - z2 * za + z2 * zb + za * zb - zb * zb + sqrt(x1 * x1 * xa * xa - 2 * x1 * x1 * xa * xb + x1 * x1 * xb * xb - x1 * x1 * za * za + 2 * x1 * x1 * za * zb - x1 * x1 * zb * zb - 2 * x1 * x2 * xa * xa + 4 * x1 * x2 * xa * xb - 2 * x1 * x2 * xb * xb + 2 * x1 * xa * y1 * ya - 2 * x1 * xa * y1 * yb - 2 * x1 * xa * y2 * ya + 2 * x1 * xa * y2 * yb + 2 * x1 * xa * z2 * za - 2 * x1 * xa * z2 * zb - 2 * x1 * xa * za * zb + 2 * x1 * xa * zb * zb - 2 * x1 * xb * y1 * ya + 2 * x1 * xb * y1 * yb + 2 * x1 * xb * y2 * ya - 2 * x1 * xb * y2 * yb - 2 * x1 * xb * z2 * za + 2 * x1 * xb * z2 * zb + 2 * x1 * xb * za * za - 2 * x1 * xb * za * zb + x2 * x2 * xa * xa - 2 * x2 * x2 * xa * xb + x2 * x2 * xb * xb + x2 * x2 * za * za - 2 * x2 * x2 * za * zb + x2 * x2 * zb * zb - 2 * x2 * xa * y1 * ya + 2 * x2 * xa * y1 * yb + 2 * x2 * xa * y2 * ya - 2 * x2 * xa * y2 * yb - 2 * x2 * xa * z2 * za + 2 * x2 * xa * z2 * zb + 2 * x2 * xa * za * zb - 2 * x2 * xa * zb * zb + 2 * x2 * xb * y1 * ya - 2 * x2 * xb * y1 * yb - 2 * x2 * xb * y2 * ya + 2 * x2 * xb * y2 * yb + 2 * x2 * xb * z2 * za - 2 * x2 * xb * z2 * zb - 2 * x2 * xb * za * za + 2 * x2 * xb * za * zb + y1 * y1 * ya * ya - 2 * y1 * y1 * ya * yb + y1 * y1 * yb * yb - y1 * y1 * za * za + 2 * y1 * y1 * za * zb - y1 * y1 * zb * zb - 2 * y1 * y2 * ya * ya + 4 * y1 * y2 * ya * yb - 2 * y1 * y2 * yb * yb + 2 * y1 * ya * z2 * za - 2 * y1 * ya * z2 * zb - 2 * y1 * ya * za * zb + 2 * y1 * ya * zb * zb - 2 * y1 * yb * z2 * za + 2 * y1 * yb * z2 * zb + 2 * y1 * yb * za * za - 2 * y1 * yb * za * zb + y2 * y2 * ya * ya - 2 * y2 * y2 * ya * yb + y2 * y2 * yb * yb + y2 * y2 * za * za - 2 * y2 * y2 * za * zb + y2 * y2 * zb * zb - 2 * y2 * ya * z2 * za + 2 * y2 * ya * z2 * zb + 2 * y2 * ya * za * zb - 2 * y2 * ya * zb * zb + 2 * y2 * yb * z2 * za - 2 * y2 * yb * z2 * zb - 2 * y2 * yb * za * za + 2 * y2 * yb * za * zb))/(za*za - 2*za*zb + zb*zb);
    double tt21 = (-x1 * xa + x1 * xb + x2 * xa - x2 * xb - y1 * ya + y1 * yb + y2 * ya - y2 * yb + z1 * za - z1 * zb - za * zb + zb * zb + sqrt(x1 * x1 * xa * xa - 2 * x1 * x1 * xa * xb + x1 * x1 * xb * xb + x1 * x1 * za * za - 2 * x1 * x1 * za * zb + x1 * x1 * zb * zb - 2 * x1 * x2 * xa * xa + 4 * x1 * x2 * xa * xb - 2 * x1 * x2 * xb * xb + 2 * x1 * xa * y1 * ya - 2 * x1 * xa * y1 * yb - 2 * x1 * xa * y2 * ya + 2 * x1 * xa * y2 * yb - 2 * x1 * xa * z1 * za + 2 * x1 * xa * z1 * zb + 2 * x1 * xa * za * zb - 2 * x1 * xa * zb * zb - 2 * x1 * xb * y1 * ya + 2 * x1 * xb * y1 * yb + 2 * x1 * xb * y2 * ya - 2 * x1 * xb * y2 * yb + 2 * x1 * xb * z1 * za - 2 * x1 * xb * z1 * zb - 2 * x1 * xb * za * za + 2 * x1 * xb * za * zb + x2 * x2 * xa * xa - 2 * x2 * x2 * xa * xb + x2 * x2 * xb * xb - x2 * x2 * za * za + 2 * x2 * x2 * za * zb - x2 * x2 * zb * zb - 2 * x2 * xa * y1 * ya + 2 * x2 * xa * y1 * yb + 2 * x2 * xa * y2 * ya - 2 * x2 * xa * y2 * yb + 2 * x2 * xa * z1 * za - 2 * x2 * xa * z1 * zb - 2 * x2 * xa * za * zb + 2 * x2 * xa * zb * zb + 2 * x2 * xb * y1 * ya - 2 * x2 * xb * y1 * yb - 2 * x2 * xb * y2 * ya + 2 * x2 * xb * y2 * yb - 2 * x2 * xb * z1 * za + 2 * x2 * xb * z1 * zb + 2 * x2 * xb * za * za - 2 * x2 * xb * za * zb + y1 * y1 * ya * ya - 2 * y1 * y1 * ya * yb + y1 * y1 * yb * yb + y1 * y1 * za * za - 2 * y1 * y1 * za * zb + y1 * y1 * zb * zb - 2 * y1 * y2 * ya * ya + 4 * y1 * y2 * ya * yb - 2 * y1 * y2 * yb * yb - 2 * y1 * ya * z1 * za + 2 * y1 * ya * z1 * zb + 2 * y1 * ya * za * zb - 2 * y1 * ya * zb * zb + 2 * y1 * yb * z1 * za - 2 * y1 * yb * z1 * zb - 2 * y1 * yb * za * za + 2 * y1 * yb * za * zb + y2 * y2 * ya * ya - 2 * y2 * y2 * ya * yb + y2 * y2 * yb * yb - y2 * y2 * za * za + 2 * y2 * y2 * za * zb - y2 * y2 * zb * zb + 2 * y2 * ya * z1 * za - 2 * y2 * ya * z1 * zb - 2 * y2 * ya * za * zb + 2 * y2 * ya * zb * zb - 2 * y2 * yb * z1 * za + 2 * y2 * yb * z1 * zb + 2 * y2 * yb * za * za - 2 * y2 * yb * za * zb))/(za*za - 2*za*zb + zb*zb);
    double tt22 = -(x1 * xa - x1 * xb - x2 * xa + x2 * xb + y1 * ya - y1 * yb - y2 * ya + y2 * yb - z1 * za + z1 * zb + za * zb - zb * zb + sqrt(x1 * x1 * xa * xa - 2 * x1 * x1 * xa * xb + x1 * x1 * xb * xb + x1 * x1 * za * za - 2 * x1 * x1 * za * zb + x1 * x1 * zb * zb - 2 * x1 * x2 * xa * xa + 4 * x1 * x2 * xa * xb - 2 * x1 * x2 * xb * xb + 2 * x1 * xa * y1 * ya - 2 * x1 * xa * y1 * yb - 2 * x1 * xa * y2 * ya + 2 * x1 * xa * y2 * yb - 2 * x1 * xa * z1 * za + 2 * x1 * xa * z1 * zb + 2 * x1 * xa * za * zb - 2 * x1 * xa * zb * zb - 2 * x1 * xb * y1 * ya + 2 * x1 * xb * y1 * yb + 2 * x1 * xb * y2 * ya - 2 * x1 * xb * y2 * yb + 2 * x1 * xb * z1 * za - 2 * x1 * xb * z1 * zb - 2 * x1 * xb * za * za + 2 * x1 * xb * za * zb + x2 * x2 * xa * xa - 2 * x2 * x2 * xa * xb + x2 * x2 * xb * xb - x2 * x2 * za * za + 2 * x2 * x2 * za * zb - x2 * x2 * zb * zb - 2 * x2 * xa * y1 * ya + 2 * x2 * xa * y1 * yb + 2 * x2 * xa * y2 * ya - 2 * x2 * xa * y2 * yb + 2 * x2 * xa * z1 * za - 2 * x2 * xa * z1 * zb - 2 * x2 * xa * za * zb + 2 * x2 * xa * zb * zb + 2 * x2 * xb * y1 * ya - 2 * x2 * xb * y1 * yb - 2 * x2 * xb * y2 * ya + 2 * x2 * xb * y2 * yb - 2 * x2 * xb * z1 * za + 2 * x2 * xb * z1 * zb + 2 * x2 * xb * za * za - 2 * x2 * xb * za * zb + y1 * y1 * ya * ya - 2 * y1 * y1 * ya * yb + y1 * y1 * yb * yb + y1 * y1 * za * za - 2 * y1 * y1 * za * zb + y1 * y1 * zb * zb - 2 * y1 * y2 * ya * ya + 4 * y1 * y2 * ya * yb - 2 * y1 * y2 * yb * yb - 2 * y1 * ya * z1 * za + 2 * y1 * ya * z1 * zb + 2 * y1 * ya * za * zb - 2 * y1 * ya * zb * zb + 2 * y1 * yb * z1 * za - 2 * y1 * yb * z1 * zb - 2 * y1 * yb * za * za + 2 * y1 * yb * za * zb + y2 * y2 * ya * ya - 2 * y2 * y2 * ya * yb + y2 * y2 * yb * yb - y2 * y2 * za * za + 2 * y2 * y2 * za * zb - y2 * y2 * zb * zb + 2 * y2 * ya * z1 * za - 2 * y2 * ya * z1 * zb - 2 * y2 * ya * za * zb + 2 * y2 * ya * zb * zb - 2 * y2 * yb * z1 * za + 2 * y2 * yb * z1 * zb + 2 * y2 * yb * za * za - 2 * y2 * yb * za * zb))/(za*za - 2*za*zb + zb*zb);
    double t4 = 0.5*((x1*x1-2*x1*xb - x2*x2 + 2*x2*xb + y1*y1 - 2 *y1*yb - y2*y2 + 2*y2*yb)/(x1*xa - x1*xb - x2*xa + x2*xb + y1*ya - y1*yb - y2*ya + y2*yb));


    if ((H1.distance2(A)-H2.distance2(A))*(H1.distance2(B)-H2.distance2(B))<=0)

    {
        if (H1.m_pt[2] < H2.m_pt[2]) {
            a = H1.m_pt[2]; b= H2.m_pt[2];
            t2=t21; t3=t22;
        } else {
            b = H1.m_pt[2]; a = H2.m_pt[2];
            t2=tt21; t3=tt22;
        }
        std::cout<<" "<<t1<<" "<<t2<<" "<<t3<< " "<<t4<<std::endl;

        p = t1*A+(1-t1)*B;
        if (t1>=0 && t1 <= 1 && p[2]<a)
            return p;
        else {
            p = t2*A+(1-t2)*B;
            if (t2>=0 && t2 <= 1 && p[2]>= a && p[2]<=b)
                return p;
            else {
                p = t3*A+(1-t3)*B;
                if (t3>=0 && t3 <= 1 && p[2]>= a && p[2]<=b)
                    return p;
                else {
                    p = t4*A+(1-t4)*B;
                    if (t4>=0 && t4 <= 1 && p[2]>= b)
                        return p;
                }
            }
        }


    }
    else
        std::cout<<"La mediatrice de " << H1 << "et de " << H2 << " ne coupe pas le segment" "[ " << A << "," << B <<"]" << std::endl;

}

/*!
 * \brief equidist
 * \param H1
 * \param H2
 * \param H3
 * \param A
 * \param B
 * \param info: output variable  (0 OK, 1 outside face, 2 segment in face)
 * \return the point on the trisectrice of H1, H2 and H3 and on the face [A,B]
 */
mmx::point<double> equidist(const hline& H1 , const hline& H2, const hline& H3,
                            const mmx::point<double>& A, const mmx::point<double>& B,
                            int& info)


{

    hline L1,L2,L3;
    mmx::point<double> q;
    info = 0;
    double a,b, c,t3,t2,x1,x2,x3,y1,y2,y3,z1,z2,z3,xa,xb,ya,yb,za,zb;

    if ( H1.m_pt[2]   <  H2.m_pt[2]  &&  H2.m_pt[2] <  H3.m_pt[2])
    {
        L1 = H1;
        L2 = H2;
        L3 = H3;
    }
    else if ( H1.m_pt[2]   <  H3.m_pt[2]  &&  H3.m_pt[2] <  H2.m_pt[2])
    {
        L1 = H1;
        L2 = H3;
        L3 = H2;
    }
    else if ( H2.m_pt[2]   <  H1.m_pt[2]  &&  H1.m_pt[2] <  H3.m_pt[2])
    {
        L1 = H2;
        L2 = H1;
        L3 = H3;
    }
    else if ( H2.m_pt[2]   <  H3.m_pt[2]  &&  H3.m_pt[2] <  H1.m_pt[2])
    {
        L1 = H2;
        L2 = H3;
        L3 = H1;
    }
    else if ( H3.m_pt[2]   <  H1.m_pt[2]  &&  H1.m_pt[2] <  H2.m_pt[2])
    {
        L1 = H3;
        L2 = H1;
        L3 = H2;
    }
    else if ( H3.m_pt[2]   <  H2.m_pt[2]  &&  H2.m_pt[2] <  H1.m_pt[2])
    {
        L1 = H3;
        L2 = H2;
        L3 = H1;
    }
    print (L1); print(L2);print (L3);
    x1=L1.m_pt[0]; y1=L1.m_pt[1]; z1=L1.m_pt[2];
    x2=L2.m_pt[0]; y2=L2.m_pt[1]; z2=L2.m_pt[2];
    x3=L3.m_pt[0]; y3=L3.m_pt[1]; z3=L3.m_pt[2];
    xa=A[0]; ya=A[1]; za=A[2];
    xb=B[0]; yb=B[1]; zb=B[2];
    double a1= (y1 * z2 - y1 * z3 - y2 * z1 + y2 * z3 + y3 * z1 - y3 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  b1  =(0.5) * (x1 * x1 * y2 - x1 * x1 * y3 - x2 * x2 * y1 + x2 * x2 * y3 + x3 * x3 * y1 - x3 * x3 * y2 + y1 * y1 * y2 - y1 * y1 * y3 - y1 * y2 * y2 + y1 * y3 * y3 - y1 * z2 * z2 + y1 * z3 * z3 + y2 * y2 * y3 - y2 * y3 * y3 + y2 * z1 * z1 - y2 * z3 * z3 - y3 * z1 * z1 + y3 * z2 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  c1=-(x1 * z2 - x1 * z3 - x2 * z1 + x2 * z3 + x3 * z1 - x3 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  d1=(-0.5) * (x1 * x1 * x2 - x1 * x1 * x3 - x1 * x2 * x2 + x1 * x3 * x3 - x1 * y2 * y2 + x1 * y3 * y3 - x1 * z2 * z2 + x1 * z3 * z3 + x2 * x2 * x3 - x2 * x3 * x3 + x2 * y1 * y1 - x2 * y3 * y3 + x2 * z1 * z1 - x2 * z3 * z3 - x3 * y1 * y1 + x3 * y2 * y2 - x3 * z1 * z1 + x3 * z2 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double   a2=(-0.5) * (y2 - y3) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  b2=(y1 * z2 - y1 * z3 + y2 * z3 - y3 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  c2=(0.5) * (x1 * x1 * y2 - x1 * x1 * y3 - x2 * x2 * y1 + x2 * x2 * y3 + x3 * x3 * y1 - x3 * x3 * y2 + y1 * y1 * y2 - y1 * y1 * y3 - y1 * y2 * y2 + y1 * y3 * y3 - y1 * z2 * z2 + y1 * z3 * z3 + y2 * y2 * y3 - y2 * y3 * y3 - y2 * z3 * z3 + y3 * z2 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double d2=(0.5) * (x2 - x3) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double e2=-(x1 * z2 - x1 * z3 + x2 * z3 - x3 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double f2=(-0.5) * (x1 * x1 * x2 - x1 * x1 * x3 - x1 * x2 * x2 + x1 * x3 * x3 - x1 * y2 * y2 + x1 * y3 * y3 - x1 * z2 * z2 + x1 * z3 * z3 + x2 * x2 * x3 - x2 * x3 * x3 + x2 * y1 * y1 - x2 * y3 * y3 - x2 * z3 * z3 - x3 * y1 * y1 + x3 * y2 * y2 + x3 * z2 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  a3=(0.5) * (y1 - y2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  b3=-z3 * (y1 - y2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double c3=(0.5) * (x1 * x1 * y2 - x1 * x1 * y3 - x2 * x2 * y1 + x2 * x2 * y3 + x3 * x3 * y1 - x3 * x3 * y2 + y1 * y1 * y2 - y1 * y1 * y3 - y1 * y2 * y2 + y1 * y3 * y3 + y1 * z3 * z3 + y2 * y2 * y3 - y2 * y3 * y3 - y2 * z3 * z3) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double d3=(-0.5) * (x1 - x2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  e3=z3 * (x1 - x2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  f3=(-0.5) * (x1 * x1 * x2 - x1 * x1 * x3 - x1 * x2 * x2 + x1 * x3 * x3 - x1 * y2 * y2 + x1 * y3 * y3 + x1 * z3 * z3 + x2 * x2 * x3 - x2 * x3 * x3 + x2 * y1 * y1 - x2 * y3 * y3 - x2 * z3 * z3 - x3 * y1 * y1 + x3 * y2 * y2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double a4=(0.5)* (x1 * x1 * y2 - x1 * x1 * y3 - x2 * x2 * y1 + x2 * x2 * y3 + x3 * x3 * y1 - x3 * x3 * y2 + y1 * y1 * y2 - y1 * y1 * y3 - y1 * y2 * y2 + y1 * y3 * y3 + y2 * y2 * y3 - y2 * y3 * y3) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double b4=(-0.5)* (x1 * x1 * x2 - x1 * x1 * x3 - x1 * x2 * x2 + x1 * x3 * x3 - x1 * y2 * y2 + x1 * y3 * y3 + x2 * x2 * x3 - x2 * x3 * x3 + x2 * y1 * y1 - x2 * y3 * y3 - x3 * y1 * y1 + x3 * y2 * y2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double z0;
    using std::min;
    using std::max;
    /*std::cout<<" "<<t1<<" "<<t2<<" "<<t3<< " "<<t4<<std::endl;*/

    if (xa==xb)

    {
        mdebug()<< "x constant";
        z0=(xa - b1)/a1;
        q=mmx::point<double>(xa,c1*z0+d1,z0);
        if (z0<z1 && min(ya,yb)<=c1*z0+d1 && c1*z0+d1<=max(ya,yb))
            return q;

        else if (0<=b2*b2 - (4.0)*a2*(c2 - xa))
        {
            z0=(0.5)*(-b2 - sqrt(b2*b2 - (4.0)*a2*(c2 - xa)))/a2;
            q=mmx::point<double>(xa,d2*z0*z0+e2*z0+f2,z0);
            if (min(ya,yb)<=d2*z0*z0+e2*z0+f2 && d2*z0*z0+e2*z0+f2<=max(ya,yb) && z1<=z0 and z0<=z2)
                return q;

            z0=(0.5)*(-b2+sqrt(b2*b2 - (4.0)*a2*(c2 - xa)))/a2;
            q=mmx::point<double>(xa,d2*z0*z0+e2*z0+f2,z0);
            if (min(ya,yb)<=d2*z0*z0+e2*z0+f2 && d2*z0*z0+e2*z0+f2<=max(ya,yb) and z1<=z0 && z0<=z2)
                return q;
        }

        else  if (0<=b3*b3 - (4.0)*a3*(c3 - xa))
        {
            z0=(0.5)*(-b3 - sqrt(b3*b3 - (4.0)*a3*(c3 - xa)))/a3;
            q==mmx::point<double>(xa,d3*z0*z0+e3*z0+f3,z0);
            if (min(ya,yb)<=d3*z0*z0+e3*z0+f3 && d3*z0*z0+e3*z0+f3<=max(ya,yb) && z2<=z0 && z0<=z3)
                return q;

            z0=(0.5)*(-b3+sqrt(b3*b3 - (4.0)*a3*(c3 - xa)))/a3;
            q=mmx::point<double>(xa,d3*z0*z0+e3*z0+f3,z0);
            if (min(ya,yb)<=d3*z0*z0+e3*z0+f3 && d3*z0*z0+e3*z0+f3<=max(ya,yb) && z2<=z0 && z0<=z3)
                return q;
        }
        else if (a4=xa && z3<=max(za,zb) && min(ya,yb)<=b4 && b4<=max(ya,yb))
        {
            info = 2;
            std::cout<<"L'intersection de la trissectrice de " << H1 << ", " << H2 << "et"  << H3 << "avec la face" "[ " << A << " " << B <<"]" "est le segment"  "[ " << z3 << " " << max(za,zb) <<"]"  << std::endl;
            return mmx::point<double>(0,0,0);
        }

        else
        {

            info = 1;
            std::cout<<"La trissectrice de " << H1 << ", " << H2 << "et"  << H3 << "ne coupe pas la face" "[ " << A << " " << B <<"]" << std::endl;
            return mmx::point<double>(0,0,0);

        }
    }


    else if(ya==yb)
    {
        mdebug()<< "y constant";
        z0=(ya - d1)/c1;
        if (z0<z1)
        {
            q=mmx::point<double>(a1*z0+b1,ya,z0);
            if (min(xa,xb)<=a1*z0+b1 && a1*z0+b1<=max(xa,xb))
                return q;
        }
        else
        {

            if (0<=e2*e2 - (4.0)*d2*(f2 - ya))
            {
                z0=0.5*(-e2 - sqrt(e2*e2 - (4.0)*d2*(f2 - ya)))/d2;
                q=mmx::point<double>(a2*z0*z0+b2*z0+c2,ya,z0);
                if (min(xa,xb)<=a2*z0*z0+b2*z0+c2 && a2*z0*z0+b2*z0+c2<=max(xa,xb) && z1<=z0 && z0<=z2)
                    return q;

                z0=0.5*(-e2+sqrt(e2*e2 - (4.0)*d2*(f2 - ya)))/d2;
                q=mmx::point<double>(a2*z0*z0+b2*z0+c2,ya,z0);
                if (min(xa,xb)<=a2*z0*z0+b2*z0+c2 && a2*z0*z0+b2*z0+c2<=max(xa,xb) && z1<=z0 && z0<=z2)
                    return  q;
            }

            else   if (0<=e3*e3 - (4.0)*d3*(f3 - ya))
            {
                z0=0.5*(-e3 - sqrt(e3*e3 - (4.0)*d3*(f3 - ya)))/d3;
                q=mmx::point<double>(a3*z0*z0+b3*z0+c3,ya,z0);
                if (min(xa,xb)<=a3*z0*z0+b3*z0+c3 && a3*z0*z0+b3*z0+c3<=max(xa,xb) && z2<=z0 && z0<=z3)
                    return q;

                z0=0.5*(-e3+sqrt(e3*e3 - (4.0)*d3*(f3 - ya)))/d3;
                q=mmx::point<double>(a3*z0*z0+b3*z0+c3,ya,z0);
                if (min(xa,xb)<=a3*z0*z0+b3*z0+c3 && a3*z0*z0+b3*z0+c3<=max(xa,xb) && z2<=z0 && z0<=z3)
                    return q;
            }
            else
            {
                if (b4=ya && z3<=max(za,zb) && min(xa,xb)<=a4 && a4<=max(xa,xb))
                {
                    info = 2;
                    std::cout<<"L'intersection de la trissectrice de " << H1 << ", " << H2 << "et"  << H3 << "avec la face" "[ " << A << " " << B <<"]" "est le segment"  "[ " << z3 << " " << max(za,zb) <<"]"  << std::endl;
                    return mmx::point<double>(0,0,0);
                }
                else
                {

                    info = 1;
                    std::cout<<"La trissectrice de " << H1 << ", " << H2 << "et"  << H3 << "ne coupe pas la face" "[ " << A << " " << B <<"]" << std::endl;
                    return mmx::point<double>(0,0,0);

                }
            }

        }

    }







    else  if (za==zb)

    {
        z0=za;
        q=mmx::point<double>(a1*z0+b1,c1*z0+d1,z0);
        if (z0<=z1 && min(xa,xb)<=a1*z0+b1 && a1*z0+b1<=max(xa,xb) && min(ya,yb)<=c1*z0+d1 && c1*z0+d1<=max(ya,yb))
            return q;

        q=mmx::point<double>(a2*z0*z0+b2*z0+c2,d2*z0*z0+e2*z0+f2,z0);
        if (z1<=z0 && z0<=z2 && min(xa,xb)<=a2*z0*z0+b2*z0+c2 && a2*z0*z0+b2*z0+c2<=max(xa,xb) && min(ya,yb)<=d2*z0*z0+e2*z0+f2 && d2*z0*z0+e2*z0+f2<=max(ya,yb))
            return q;


        q=mmx::point<double>(a3*z0*z0+b3*z0+c3,d3*z0*z0+e3*z0+f3,z0);
        if (z2<=z0 && z0<=z3 && min(xa,xb)<=a3*z0*z0+b3*z0+c3 && a3*z0*z0+b3*z0+c3<=max(xa,xb) && min(ya,yb)<=d3*z0*z0+e3*z0+f3 && d3*z0*z0+e3*z0+f3<=max(ya,yb))
            return q;

        q=mmx::point<double>(a4,b4,z0);
        if (z3<=z0 && min(xa,xb)<=a4 && a4<=max(xa,xb) && min(ya,yb)<=b4 && b4<=max(ya,yb))
            return q;

        info = 1;
        std::cout<<"La trissectrice de " << H1 << ", " << H2 << "et"  << H3 << "ne coupe pas la face" "[ " << A << " " << B <<"]" << std::endl;
        return mmx::point<double>(0,0,0);
    }







}








/*!
 * \brief equidist
 * \param H1
 * \param H2
 * \param H3
 * \param H4
 * \param info: output variable
 * \return the point equidistant point from the four hlines H1, H2,H3 and H4.
 */
mmx::point<double> equidist(const hline& H1 , const hline& H2, const hline& H3,const hline& H4, int& info)


{

    hline L1,L2,L3,L4;
    mmx::point<double> q;
    info = 0;
    double A1,B1,C1,A2,B2,C2,A3,B3,C3, D1,D2,D3,a, b, c,t3,t2,x1,x2,x3,y1,y2,y3,z1,z2,z3,x4,y4,z4;


    if (H1.m_pt[2]<H2.m_pt[2] && H2.m_pt[2]<H3.m_pt[2] && H3.m_pt[2]<H4.m_pt[2])
    {
        L1 = H1;
        L2 = H2;
        L3 = H3;
        L4=  H4;
    }
    else if (H1.m_pt[2]< H2.m_pt[2] && H2.m_pt[2]<H4.m_pt[2] && H4.m_pt[2]<H3.m_pt[2])
    {
        L1 = H1;
        L2 = H2;
        L3 = H4;
        L4 = H3;
    }
    else if (H1.m_pt[2]< H3.m_pt[2] && H3.m_pt[2]<H2.m_pt[2] && H2.m_pt[2]<H4.m_pt[2])
    {
        L1 = H1;
        L2 = H3;
        L3 = H2;
        L4= H4;
    }

    else if (H1.m_pt[2]< H3.m_pt[2] && H3.m_pt[2]<H4.m_pt[2] && H4.m_pt[2]<H2.m_pt[2])
    {
        L1 = H1;
        L2 = H3;
        L3 = H4;
        L4= H2;
    }

    else if (H1.m_pt[2]< H4.m_pt[2] && H4.m_pt[2]<H2.m_pt[2] && H2.m_pt[2]<H3.m_pt[2])
    {
        L1 = H1;
        L2 = H4;
        L3 = H2;
        L4 = H3;
    }


    else if (H1.m_pt[2]< H4.m_pt[2] && H4.m_pt[2]<H3.m_pt[2] && H3.m_pt[2]<H2.m_pt[2])
    {
        L1 = H1;
        L2 = H4;
        L3 = H3;
        L4 = H2;
    }


    else if (H2.m_pt[2]< H1.m_pt[2] && H1.m_pt[2]<H3.m_pt[2] && H3.m_pt[2]<H4.m_pt[2])
    {
        L1 = H2;
        L2 = H1;
        L3 = H3;
        L4 = H4;
    }


    else if (H2.m_pt[2]< H1.m_pt[2] && H1.m_pt[2]<H4.m_pt[2] && H4.m_pt[2]<H3.m_pt[2])
    {
        L1 = H2;
        L2 = H1;
        L3 = H4;
        L4 = H3;
    }


    else if (H2.m_pt[2]< H3.m_pt[2] && H3.m_pt[2]<H1.m_pt[2] && H1.m_pt[2]<H4.m_pt[2])
    {
        L1 = H2;
        L2 = H3;
        L3 = H1;
        L4 = H4;
    }



    else if (H2.m_pt[2]< H3.m_pt[2] && H3.m_pt[2]<H4.m_pt[2] && H4.m_pt[2]<H1.m_pt[2])
    {
        L1 = H2;
        L2 = H3;
        L3 = H4;
        L4 = H1;
    }



    else if (H2.m_pt[2]< H4.m_pt[2] && H4.m_pt[2]<H1.m_pt[2] && H1.m_pt[2]<H3.m_pt[2])
    {
        L1 = H2;
        L2 = H4;
        L3 = H1;
        L4 = H3;
    }


    else if (H2.m_pt[2]< H4.m_pt[2] && H4.m_pt[2]<H3.m_pt[2] && H3.m_pt[2]<H1.m_pt[2])
    {
        L1 = H2;
        L2 = H4;
        L3 = H3;
        L4 = H1;
    }


    else if (H3.m_pt[2]< H1.m_pt[2] && H1.m_pt[2]<H2.m_pt[2] && H2.m_pt[2]<H4.m_pt[2])
    {
        L1 = H3;
        L2 = H1;
        L3 = H2;
        L4 = H4;
    }


    else if (H3.m_pt[2]< H1.m_pt[2] && H1.m_pt[2]<H4.m_pt[2] && H4.m_pt[2]<H2.m_pt[2])
    {
        L1 = H3;
        L2 = H1;
        L3 = H4;
        L4 = H2;
    }


    else if (H3.m_pt[2]< H2.m_pt[2] && H2.m_pt[2]<H1.m_pt[2] && H1.m_pt[2]<H4.m_pt[2])
    {
        L1 = H3;
        L2 = H2;
        L3 = H1;
        L4 = H4;
    }


    else if (H3.m_pt[2]< H2.m_pt[2] && H2.m_pt[2]<H4.m_pt[2] && H4.m_pt[2]<H1.m_pt[2])
    {
        L1 = H3;
        L2 = H2;
        L3 = H4;
        L4 = H1;
    }


    else if (H3.m_pt[2]< H4.m_pt[2] && H4.m_pt[2]<H2.m_pt[2] && H2.m_pt[2]<H1.m_pt[2])
    {
        L1 = H3;
        L2 = H4;
        L3 = H2;
        L4 = H1;
    }


    else if (H3.m_pt[2]< H4.m_pt[2] && H4.m_pt[2]<H1.m_pt[2] && H1.m_pt[2]<H2.m_pt[2])
    {
        L1 = H3;
        L2 = H4;
        L3 = H1;
        L4 = H2;
    }


    else if (H4.m_pt[2]< H1.m_pt[2] && H1.m_pt[2]<H2.m_pt[2] && H2.m_pt[2]<H3.m_pt[2])
    {
        L1 = H4;
        L2 = H1;
        L3 = H2;
        L4 = H3;
    }


    else if (H4.m_pt[2]< H1.m_pt[2] && H1.m_pt[2]<H3.m_pt[2] && H3.m_pt[2]<H2.m_pt[2])
    {
        L1 = H4;
        L2 = H1;
        L3 = H3;
        L4 = H2;
    }


    else if (H4.m_pt[2]< H2.m_pt[2] && H2.m_pt[2]<H1.m_pt[2] && H1.m_pt[2]<H3.m_pt[2])
    {
        L1 = H4;
        L2 = H2;
        L3 = H1;
        L4 = H3;
    }


    else if (H4.m_pt[2]< H2.m_pt[2] && H2.m_pt[2]<H3.m_pt[2] && H3.m_pt[2]<H1.m_pt[2])
    {
        L1 = H4;
        L2 = H2;
        L3 = H3;
        L4 = H1;
    }

    else if (H4.m_pt[2]< H3.m_pt[2] && H3.m_pt[2]<H1.m_pt[2] && H1.m_pt[2]<H2.m_pt[2])
    {
        L1 = H4;
        L2 = H3;
        L3 = H1;
        L4 = H2;
    }


    else if (H4.m_pt[2]< H3.m_pt[2] && H3.m_pt[2]<H2.m_pt[2] && H2.m_pt[2]<H1.m_pt[2])
    {
        L1 = H4;
        L2 = H3;
        L3 = H2;
        L4 = H1;
    }

    print (L1); print(L2);print (L3);print (L4);
    x1=L1.m_pt[0]; y1=L1.m_pt[1]; z1=L1.m_pt[2];
    x2=L2.m_pt[0]; y2=L2.m_pt[1]; z2=L2.m_pt[2];
    x3=L3.m_pt[0]; y3=L3.m_pt[1]; z3=L3.m_pt[2];
    x4=L4.m_pt[0]; y4=L4.m_pt[1]; z4=L4.m_pt[2];
    double  a1= (y1 * z2 - y1 * z3 - y2 * z1 + y2 * z3 + y3 * z1 - y3 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  b1  =(0.5) * (x1 * x1 * y2 - x1 * x1 * y3 - x2 * x2 * y1 + x2 * x2 * y3 + x3 * x3 * y1 - x3 * x3 * y2 + y1 * y1 * y2 - y1 * y1 * y3 - y1 * y2 * y2 + y1 * y3 * y3 - y1 * z2 * z2 + y1 * z3 * z3 + y2 * y2 * y3 - y2 * y3 * y3 + y2 * z1 * z1 - y2 * z3 * z3 - y3 * z1 * z1 + y3 * z2 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  c1=-(x1 * z2 - x1 * z3 - x2 * z1 + x2 * z3 + x3 * z1 - x3 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  d1=(-0.5) * (x1 * x1 * x2 - x1 * x1 * x3 - x1 * x2 * x2 + x1 * x3 * x3 - x1 * y2 * y2 + x1 * y3 * y3 - x1 * z2 * z2 + x1 * z3 * z3 + x2 * x2 * x3 - x2 * x3 * x3 + x2 * y1 * y1 - x2 * y3 * y3 + x2 * z1 * z1 - x2 * z3 * z3 - x3 * y1 * y1 + x3 * y2 * y2 - x3 * z1 * z1 + x3 * z2 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  a2=(-0.5) * (y2 - y3) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  b2=(y1 * z2 - y1 * z3 + y2 * z3 - y3 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  c2=(0.5) * (x1 * x1 * y2 - x1 * x1 * y3 - x2 * x2 * y1 + x2 * x2 * y3 + x3 * x3 * y1 - x3 * x3 * y2 + y1 * y1 * y2 - y1 * y1 * y3 - y1 * y2 * y2 + y1 * y3 * y3 - y1 * z2 * z2 + y1 * z3 * z3 + y2 * y2 * y3 - y2 * y3 * y3 - y2 * z3 * z3 + y3 * z2 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  d2=(0.5) * (x2 - x3) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  e2=-(x1 * z2 - x1 * z3 + x2 * z3 - x3 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  f2=(-0.5) * (x1 * x1 * x2 - x1 * x1 * x3 - x1 * x2 * x2 + x1 * x3 * x3 - x1 * y2 * y2 + x1 * y3 * y3 - x1 * z2 * z2 + x1 * z3 * z3 + x2 * x2 * x3 - x2 * x3 * x3 + x2 * y1 * y1 - x2 * y3 * y3 - x2 * z3 * z3 - x3 * y1 * y1 + x3 * y2 * y2 + x3 * z2 * z2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  a3=(0.5) * (y1 - y2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  b3=-z3 * (y1 - y2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  c3=(0.5) * (x1 * x1 * y2 - x1 * x1 * y3 - x2 * x2 * y1 + x2 * x2 * y3 + x3 * x3 * y1 - x3 * x3 * y2 + y1 * y1 * y2 - y1 * y1 * y3 - y1 * y2 * y2 + y1 * y3 * y3 + y1 * z3 * z3 + y2 * y2 * y3 - y2 * y3 * y3 - y2 * z3 * z3) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  d3=(-0.5) * (x1 - x2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  e3=z3 * (x1 - x2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  f3=(-0.5) * (x1 * x1 * x2 - x1 * x1 * x3 - x1 * x2 * x2 + x1 * x3 * x3 - x1 * y2 * y2 + x1 * y3 * y3 + x1 * z3 * z3 + x2 * x2 * x3 - x2 * x3 * x3 + x2 * y1 * y1 - x2 * y3 * y3 - x2 * z3 * z3 - x3 * y1 * y1 + x3 * y2 * y2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  a4=(0.5)* (x1 * x1 * y2 - x1 * x1 * y3 - x2 * x2 * y1 + x2 * x2 * y3 + x3 * x3 * y1 - x3 * x3 * y2 + y1 * y1 * y2 - y1 * y1 * y3 - y1 * y2 * y2 + y1 * y3 * y3 + y2 * y2 * y3 - y2 * y3 * y3) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  b4=(-0.5)* (x1 * x1 * x2 - x1 * x1 * x3 - x1 * x2 * x2 + x1 * x3 * x3 - x1 * y2 * y2 + x1 * y3 * y3 + x2 * x2 * x3 - x2 * x3 * x3 + x2 * y1 * y1 - x2 * y3 * y3 - x3 * y1 * y1 + x3 * y2 * y2) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    double  z0;
    double t0= -(0.5)*((2 * b1 * x2 - 2 * b1 * x4 + 2 * d1 * y2 - 2 * d1 * y4 - x2 * x2 + x4 * x4 - y2 * y2 + y4 * y4 - z2 * z2 + z4 * z4) / (a1 * x2 - a1 * x4 + c1 * y2 - c1 * y4 + z2 - z4));
    A1=-2 * a2 * x2 + 2 * a2 * x4 - 2 * d2 * y2 + 2 * d2 * y4;
    B1=-2 * b2 * x2 + 2 * b2 * x4 - 2 * e2 * y2 + 2 * e2 * y4 - 2 * z2 + 2 * z4;
    C1=-2 * c2 * x2 + 2 * c2 * x4 - 2 * f2 * y2 + 2 * f2 * y4 + x2 * x2 - x4 * x4 + y2 * y2 - y4 * y4 + z2 * z2 - z4 * z4;
    D1= B1*B1-4*A1*C1;
    double  t11=(0.5)*(-B1 - sqrt(B1 * B1 - 4 * B1 * C1))/A1;
    double  t12= (0.5)*(-B1 + sqrt (B1 * B1 - 4 * B1 * C1))/A1;
    A2=-2 * a3 * x2 + 2 * a3 * x4 - 2 * d3 * y2 + 2 * d3 * y4 - 1;
    B2=-2 * b3 * x2 + 2 * b3 * x4 - 2 * e3 * y2 + 2 * e3 * y4 + 2 * z4;
    C2=-2 * c3 * x2 + 2 * c3 * x4 - 2 * f3 * y2 + 2 * f3 * y4 + x2 * x2 - x4 * x4 + y2 * y2 - y4 * y4 - z4 * z4;
    D2= B2*B2-4*A2*C2;
    A3=-1; B3=2*z4;
    C3=-2 * a4 * x2 + 2 * a4 * x4 - 2 * b4 * y2 + 2 * b4 * y4 + x2 * x2 - x4 * x4 + y2 * y2 - y4 * y4 - z4 * z4;
    D3= B3*B3-4*A3*C3;
    double  t21=(0.5)*(-B2 - sqrt(B2* B2 - 4 * B2 * C2))/A2;
    double  t22=(0.5)*(-B2 + sqrt(B2 * B2 - 4 * B2 * C2))/A2;
    double  t31= (0.5)*(-B3 - sqrt(B3* B3 - 4 * B3 * C3))/A3 ;
    double  t32= (0.5)*(- B3 + sqrt( B3 * B3 - 4 * B3 * C3))/A3;

    //#if (t0<z1) and t0 in T1 := [a1*t+b1, c1*t+d1, t];
    //#if (B1^2-4*B1*C1>=0) z1<t11, t12<z2 in T2 := [a2*t^2+b2*t+c2, d2*t^2+e2*t+f2, t];
    //#if (B2^2-4*A1*C1>=0) z2<t21, t22<z3 in T3 := [a3*t^2+b3*t+c3, d3*t^2+e3*t+f3, t];
    //#if (B3^2-4*A3*C3>=0) z3<t31 ,t32<z4  in T4 := [a4, b4, t];
    using std::min;
    using std::max;


    q=mmx::point<double>(a1*t0+b1, c1*t0+d1, t0);
    if (t0<= z1)
    {
        return q;
    }
    else
    {
        if(D1>=0)
        {
            q=mmx::point<double>(a2*t11*t11+b2*t11+c2, d2*t11*t11+e2*t11+f2, t11);
            if (z1<=t11 && t11<=z2)
                return q;
            q=mmx::point<double>(a2*t12*t12+b2*t12+c2, d2*t12*t12+e2*t12+f2, t12);
            if (z1<=t12 && t12<=z2)
                return q;
        }
        else       if(D2>=0)
        {
            q=mmx::point<double>(a3*t21*t21+b3*t21+c3, d3*t21*t21+e3*t21+f3, t21);
            if (z2<=t21 && t21<=z3)
                return q;
            q=mmx::point<double>(a3*t22*t22+b3*t22+c3, d3*t22*t22+e3*t22+f3, t22);
            if (z2<=t22 && t22<=z3)
                return q;
        }
        else if(D3>=0)
        {
            q=mmx::point<double>(a4, b4, t31);
            if (z3<=t31 && t31<=z4)
                return q;
            q=mmx::point<double>(a4, b4, t32);
            if (z3<=t32 && t32<=z4)
                return q;
        }
    }

}
































