#include <math.h>
#include <voronoi/hline.h>
#include <voronoi/point.hpp>

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
        //std::cout<<" "<<t1<<" "<<t2<<" "<<t3<< " "<<t4<<std::endl;
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
 * \return the point on the trisectrice of H1, H2 and H3 and on the face [A,B]
 */
mmx::point<double> equidist(const hline& H1 , const hline& H2, const hline& H3,
                            const mmx::point<double>& A, const mmx::point<double>& B)


{
    //mmx::point<double> C,D,q;


 /*   if (A[0]=B[0]) {

    } else if (A[1]=B[1]) {

    } else if (A[2]=B[2]) {

    }


*/

}







