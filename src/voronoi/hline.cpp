#include <math.h>
#include <voronoi/mdebug.h>
#include <voronoi/hline.h>

/*!
 * \brief hline::distance2
 * \param p
 * \return the squared euclidean distance from the half-line to a point p.
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
    //mdebug()<<"constructor"<<X<<Y<<Z;
    return this->distance2( mmx::point<double>(X,Y,Z) );
}

/*!
 * \brief hline::distance2
 * \param A
 * \param B
 * \return the euclidean distance from the half-line to a segment [AB].
 */
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



/*!
 * \brief hline::distance2
 * \param H
 * \param A
 * \param B
 * \return the euclidean distance from the half-line H to a face of [AB].
 */

double hline::distance2_to_face(const mmx::point<double>& A, const mmx::point<double>& B) const {



    mmx::point<double> O(m_pt[0], m_pt[1],m_pt[2]);
    mmx::point<double> u(m_dir[0],m_dir[1],m_dir[2]);
    mmx::point<double> F1,F2,F3,F4,F,v0,w0,p0,G,P,Q,R,S,h,D;
    double x1,y1,z1,xa,ya,za,xb,yb,zb,a1,b1,c1,a,b,xo,yo,zo,c,x,y,z,d,d1,d2,s,t,e,e1,e2,w1,v1,t0;
    xo=O[0]; yo=O[1]; zo=O[2];
    xa=A[0]; ya=A[1]; za=A[2];
    xb=B[0]; yb=B[1]; zb=B[2];
    a=u[0];b=u[1];c=u[2];
    using std::min;
    using std::max;



    if (A[0]==B[0])

    {
        x=A[0];
        F1 = mmx::point<double>(x,min(A[1],B[1]),min(A[2],B[2]));
        F2 = mmx::point<double>(x,min(A[1],B[1]),max(A[2],B[2]));
        F3 = mmx::point<double>(x,max(A[1],B[1]),max(A[2],B[2]));
        F4 = mmx::point<double>(x,max(A[1],B[1]),min(A[2],B[2]));
        std::cout<< " "<<F1<<" "<<F2<< " "<<F3<< " "<<F4<<std::endl;
    }

    if (A[1]==B[1])

    {

        y=A[1];
        F1 = mmx::point<double>(min(A[0],B[0]),y,min(A[2],B[2]));
        F2 = mmx::point<double>(min(A[0],B[0]),y,max(A[2],B[2]));
        F3 = mmx::point<double>(max(A[0],B[0]),y,max(A[2],B[2]));
        F4 = mmx::point<double>(max(A[0],B[0]),y,min(A[2],B[2]));
        std::cout<< " "<<F1<<" "<<F2<< " "<<F3<< " "<<F4<<std::endl;
    }

    if (A[2]==B[2])

    {

        z=A[2];
        F1 = mmx::point<double>(min(A[0],B[0]),min(A[1],B[1]),z);
        F2 = mmx::point<double>(min(A[0],B[0]),max(A[1],B[1]),z);
        F3 = mmx::point<double>(max(A[0],B[0]),max(A[1],B[1]),z);
        F4 = mmx::point<double>(max(A[0],B[0]),min(A[1],B[1]),z);
        std::cout<< " "<<F1<<" "<<F2<< " "<<F3<< " "<<F4<<std::endl;
    }

    v0=F2-F1;
    w0=F4-F1;
    p0=v0.cross(w0);
    G=F1-O;
    e1=G.dot(v0);
    e2=G.dot(w0);
    e=v0.dot(w0);
    v1=v0.dot(v0);
    w1=w0.dot(w0);
    F=B-A;
    R=O-A;
    S=F1-O;
    s=(e1*w1-e2*e)/(v1*w1-e*e);
    t=(e2*v1-e1*e)/(v1*w1-e*e);



    if((u.cross(F)).norm()==0)
    {

        if (s>=0 && s<=1 && t>=0 && t<=1)
        {
            if(p0.dot(R)!=0)
            {
                P=F1+s*v0+t*w0;
                Q=P-O;
                d=sqrt(Q.dot(Q));
                return d;
            }

            else
            {
                d=0;
                return d;
            }
        }
        else
        {
            d1= min(this->distance2(F1, F2),this->distance2(F2, F3));
            d2= min(this->distance2(F3, F4),this->distance2(F4, F1));
            d=min(d1,d2);
            std::cout<< " "<<d1<<" "<<d2<<std::endl;
            return d;

        }


    }

    else

    {
        t0 = (p0.dot(S))/(u.dot(p0));
        h=O+t0*u;
        D=h-F1;


        if (t0<0 && s>=0 && s<=1 && t>=0 && t<=1 )
        {
            P=F1+s*v0+t*w0;
            Q=P-O;
            d=sqrt(Q.dot(Q));
            return d;

        }
        else
        {
            d1= min(this->distance2(F1, F2),this->distance2(F2, F3));
            d2= min(this->distance2(F3, F4),this->distance2(F4, F1));
            d=min(d1,d2);
            std::cout<< " "<<d1<<" "<<d2<<std::endl;
            return d;
        }

        if(t0>=0 && p0.dot(D)!=0)
        {
            d1= min(this->distance2(F1, F2),this->distance2(F2, F3));
            d2= min(this->distance2(F3, F4),this->distance2(F4, F1));
            d=min(d1,d2);
            std::cout<< " "<<d1<<" "<<d2<<std::endl;
            return d;

        }

        if(t0>=0 && p0.dot(D)==0)
        {
            d=0;
            return d;
        }


    }


}



/*

      if( a!=0 || b != 0)
    {
        a1 = sqrt(a * a + b * b + c * c);
        b1= sqrt(a * a + b * b);
        c1= sqrt(a * a * c * c + (a * a + b * b)*(a * a + b * b) + b * b * c * c);
        U= mmx::point<double>(a / a1, b / a1,  c / a1);
        V= mmx::point<double>(-b / b1,a / b1,  0);
        W= mmx::point<double>(-a * c / c1,-b * c / c1,  (a * a + b * b) / c1);
        A1= mmx::point<double> ( c*za/a1 +((a*a*a) * a1 + a * a1 * b * b)*xa/(c1*c1) +(a * a * a1 * b + a1 * (b*b*b)) * ya/(c1*c1), (-a*a*b*b1 - (b*b*b)* b1 - b*b1 * c * c) * xa/(c1*c1) + (a*a*a * b1 + a * b * b * b1 + a * b1 * c * c)* ya /(c1*c1),  c1 *za/(a1*a1)- a * c * xa/c1 -b*c *ya/c1);
        B1= mmx::point<double> ( c*zb/a1 +((a*a*a) * a1 + a * a1 * b * b)*xb/(c1*c1) +(a * a * a1 * b + a1 * (b*b*b)) * yb/(c1*c1), (-a*a*b*b1 - (b*b*b)* b1 - b*b1 * c * c) * xb/(c1*c1) + (a*a*a * b1 + a * b * b * b1 + a * b1 * c * c)* yb /(c1*c1),  c1 *zb/(a1*a1)- a * c * xb/c1 -b*c *yb/c1);
        O1= mmx::point<double> ( c*zo/a1 +((a*a*a) * a1 + a * a1 * b * b)*xo/(c1*c1) +(a * a * a1 * b + a1 * (b*b*b)) * yo/(c1*c1), (-a*a*b*b1 - (b*b*b)* b1 - b*b1 * c * c) * xo/(c1*c1) + (a*a*a * b1 + a * b * b * b1 + a * b1 * c * c)* yo /(c1*c1),  c1 *zo/(a1*a1)- a * c * xo/c1 -b*c *yo/c1);


                d= fabs(O1[1]-y);
                return d;

            }
            else
            {
                d1= min(L1.distance2(F01, F02),L1.distance2(F02, F03));
                d2= min(L1.distance2(F03, F04),L1.distance2(F04, F01));
                d=min(d1,d2);
                return d;
            }


        }
        if (A1[2]==B1[2])

        {

            z=A1[2];
            F01 = mmx::point<double>(min(A1[0],B1[0]),min(A1[1],B1[1]),z);
            F02 = mmx::point<double>(min(A1[0],B1[0]),max(A1[1],B1[1]),z);
            F03 = mmx::point<double>(max(A1[0],B1[0]),max(A1[1],B1[1]),z);
            F04 = mmx::point<double>(max(A1[0],B1[0]),min(A1[1],B1[1]),z);
              F1  = F01[0]*U+F01[1]*V+F01[2]*W;
              F2  = F02[0]*U+F02[1]*V+F02[2]*W;
              F3  = F03[0]*U+F03[1]*V+F03[2]*W;
              F4  = F04[0]*U+F04[1]*V+F04[2]*W;
            hline L1 (O1[0],O1[1],O1[2]);
            std::cout<<" "<<O1<<" "<<F01<<" "<<F02<< " "<<F03<< " "<<F1<<" "<<F2<< " "<<F3<< " "<<F4<<std::endl;
            if(O1[0]<=max(A1[0],B1[0]) && O1[1]>=min(A1[1],B1[1]) && O1[1]<=max(A1[1],B1[1]))
            {
                d= fabs(O1[2]-z);
                return d;

            }
            else
            {
                d1= min(L1.distance2(F01, F02),L1.distance2(F02, F03));
                d2= min(L1.distance2(F03, F04),L1.distance2(F04, F01));
                d=min(d1,d2);
                return d;
            }

        }
    }






    if( a!=0 || c != 0)
    {
        a1 = sqrt(a * a + b * b + c * c);
        b1 = sqrt(a * a + c * c);
        c1= sqrt(a * a * b * b + (a * a + c * c)*(a * a + c * c) + b * b * c * c);
        U= mmx::point<double>(a / a1, b / a1,  c / a1);
        V= mmx::point<double>(c / b1,0, -a / b1);
        W= mmx::point<double>(-a * b / c1,(a * a + c * c) / c1,  -b * c / c1);
        A1= mmx::point<double> (b * ya / a1 + ((a*a*a) * a1 + a * a1 * c * c) * xa /(c1*c1) - (-a * a * a1 * c - a1 * (c*c*c)) * za /(c1*c1), (a * a * b1 * c + b * b * b1 * c + b1 * (c*c*c)) * xa / (c1*c1) - ((a*a*a) * b1 + a * b * b * b1 + a * b1 * c * c) * za / (c1*c1), c1 * ya / (a1*a1) - a * b * xa / c1 - b * c * za / c1);
        B1= mmx::point<double> ( b * yb / a1 + ((a*a*a) * a1 + a * a1 * c * c) * xb /(c1*c1) - (-a * a * a1 * c - a1 * (c*c*c)) * zb /(c1*c1), (a * a * b1 * c + b * b * b1 * c + b1 * (c*c*c)) * xb / (c1*c1) - ((a*a*a) * b1 + a * b * b * b1 + a * b1 * c * c) * zb / (c1*c1), c1 * yb / (a1*a1) - a * b * xb / c1 - b * c * zb / c1);
        O1= mmx::point<double> ( b * yo / a1 + ((a*a*a) * a1 + a * a1 * c * c) * xo /(c1*c1) - (-a * a * a1 * c - a1 * (c*c*c)) * zo /(c1*c1), (a * a * b1 * c + b * b * b1 * c + b1 * (c*c*c)) * xo / (c1*c1) - ((a*a*a) * b1 + a * b * b * b1 + a * b1 * c * c) * zo / (c1*c1), c1 * yo / (a1*a1) - a * b * xo / c1 - b * c * zo / c1);

        if (A1[0]==B1[0])

        {
            x=A1[0];
            F01 = mmx::point<double>(x,min(A1[1],B1[1]),min(A1[2],B1[2]));
            F02 = mmx::point<double>(x,min(A1[1],B1[1]),max(A1[2],B1[2]));
            F03 = mmx::point<double>(x,max(A1[1],B1[1]),max(A1[2],B1[2]));
            F04 = mmx::point<double>(x,max(A1[1],B1[1]),min(A1[2],B1[2]));
            F1  = F01[0]*U+F01[1]*V+F01[2]*W;
            F2  = F02[0]*U+F02[1]*V+F02[2]*W;
            F3  = F03[0]*U+F03[1]*V+F03[2]*W;
            F4  = F04[0]*U+F04[1]*V+F04[2]*W;
            std::cout<< " "<<F1<<" "<<F2<< " "<<F3<< " "<<F4<<std::endl;
            v0=F2-F1;
            w0=F4-F1;
            p0=v0.cross(w0);
            G=F1-O;
            e1=G.dot(v0);
            e2=G.dot(w0);
            e=v0.dot(w0);
            v1=v0.dot(v0);
            w1=w0.dot(w0);
            F=B-A;
            R=O-A;
            s=(e1*w1-e2*e)/(v1*w1-e*e);
            t=(e2*v1-e1*e)/(v1*w1-e*e);
            hline L1 (O[0],O[1],O[2]);
            if(u.dot(F)==0)
            {


                if (s>=0 && s<=1 && t>=0 && t<=1)
                {
                    P=F1+s*v0+t*w0;
                    if((P[0]-O[0])/u[0]<0)
                    {
                        Q=P-O;
                        d=sqrt(Q.dot(Q));
                        return d;
                    }
                    else
                    {
                        d=0;
                        return d;
                    }
                }
                else
                {
                    d1= min(L1.distance2(F1, F2),L1.distance2(F2, F3));
                    d2= min(L1.distance2(F3, F4),L1.distance2(F4, F1));
                    d=min(d1,d2);
                    std::cout<< " "<<d1<<" "<<d2<<std::endl;
                    std::cout<< " "<<H.distance2(F1, F2)<<" "<<H.distance2(F2, F3)<< " "<<H.distance2(F3, F4)<< " "<<H.distance2(F4, F1)<<std::endl;
                    return d;
                }


            }

            if(u.cross(F)==0)
            {

                if (s>=0 && s<=1 && t>=0 && t<=1)
                {
                    if(p0.dot(R)!=0)
                    {
                        P=F1+s*v0+t*w0;
                        Q=P-O;
                        d=sqrt(Q.dot(Q));
                        return d;
                    }

                    else
                    {
                        d=0;
                        return d;
                    }
                }
                else
                {
                    d1= min(L1.distance2(F1, F2),L1.distance2(F2, F3));
                    d2= min(L1.distance2(F3, F4),L1.distance2(F4, F1));
                    d=min(d1,d2);
                    std::cout<< " "<<d1<<" "<<d2<<std::endl;
                    std::cout<< " "<<L1.distance2(F1, F2)<<" "<<L1.distance2(F2, F3)<< " "<<L1.distance2(F3, F4)<< " "<<L1.distance2(F4, F1)<<std::endl;
                    return d;

                }




            }
        }


         hline L1 (O1[0],O1[1],O1[2]);


             std::cout<<" "<<O1<<" "<<F01<<" "<<F02<< " "<<F03<< " "<<F04<< " "<<F1<<" "<<F2<< " "<<F3<< " "<<F4<<std::endl;
            if(O1[1]>=min(A1[1],B1[1]) && O1[1]<=max(A1[1],B1[1]) && O1[2]>=min(A1[2],B1[2]) && O1[2]<=max(A1[2],B1[2]))
            {
                if (O1[0]> x)
                { d= fabs(O1[0]-x);
                    return d;
                }
                else if (O1[0]<= x)
                {

                    d=0;
                     std::cout<<" "<<d<<std::endl;
                    return d;
                }
            }
            else
            {
                d1= min(L1.distance2(F01, F02),L1.distance2(F02, F03));
                d2= min(L1.distance2(F03, F04),L1.distance2(F04, F01));
                d=min(d1,d2);
                 std::cout<<" "<<d1<<" "<<d2<<std::endl;
                return d;
            }
        }

        if (A1[1]==B1[1])

        {

            y=A1[1];
            F01 = mmx::point<double>(min(A1[0],B1[0]),y,min(A1[2],B1[2]));
            F02 = mmx::point<double>(min(A1[0],B1[0]),y,max(A1[2],B1[2]));
            F03 = mmx::point<double>(max(A1[0],B1[0]),y,max(A1[2],B1[2]));
            F04 = mmx::point<double>(max(A1[0],B1[0]),y,min(A1[2],B1[2]));
            F1  = F01[0]*U+F01[1]*V+F01[2]*W;
            F2  = F02[0]*U+F02[1]*V+F02[2]*W;
            F3  = F03[0]*U+F03[1]*V+F03[2]*W;
            F4  = F04[0]*U+F04[1]*V+F04[2]*W;
            std::cout<< " "<<F1<<" "<<F2<< " "<<F3<< " "<<F4<<std::endl;
            v0=F2-F1;
            w0=F4-F1;
            p0=v0.cross(w0);
            G=F1-O;
            e1=G.dot(v0);
            e2=G.dot(w0);
            e=v0.dot(w0);
            v1=v0.dot(v0);
            w1=w0.dot(w0);
            F=B-A;
            R=O-A;
            s=(e1*w1-e2*e)/(v1*w1-e*e);
            t=(e2*v1-e1*e)/(v1*w1-e*e);
            hline L1 (O[0],O[1],O[2]);
            if(u.dot(F)==0)
            {


                if (s>=0 && s<=1 && t>=0 && t<=1)
                {
                    P=F1+s*v0+t*w0;
                    if((P[0]-O[0])/u[0]<0)
                    {
                        Q=P-O;
                        d=sqrt(Q.dot(Q));
                        return d;
                    }
                    else
                    {
                        d=0;
                        return d;
                    }
                }
                else
                {
                    d1= min(L1.distance2(F1, F2),L1.distance2(F2, F3));
                    d2= min(L1.distance2(F3, F4),L1.distance2(F4, F1));
                    d=min(d1,d2);
                    std::cout<< " "<<d1<<" "<<d2<<std::endl;
                    std::cout<< " "<<L1.distance2(F1, F2)<<" "<<L1.distance2(F2, F3)<< " "<<L1.distance2(F3, F4)<< " "<<L1.distance2(F4, F1)<<std::endl;
                    return d;
                }


            }

            if(u.cross(F)==0)
            {

                if (s>=0 && s<=1 && t>=0 && t<=1)
                {
                    if(p0.dot(R)!=0)
                    {
                        P=F1+s*v0+t*w0;
                        Q=P-O;
                        d=sqrt(Q.dot(Q));
                        return d;
                    }

                    else
                    {
                        d=0;
                        return d;
                    }
                }
                else
                {
                    d1= min(L1.distance2(F1, F2),L1.distance2(F2, F3));
                    d2= min(L1.distance2(F3, F4),L1.distance2(F4, F1));
                    d=min(d1,d2);
                    std::cout<< " "<<d1<<" "<<d2<<std::endl;
                    std::cout<< " "<<L1.distance2(F1, F2)<<" "<<L1.distance2(F2, F3)<< " "<<L1.distance2(F3, F4)<< " "<<L1.distance2(F4, F1)<<std::endl;
                    return d;

                }




            }


        }


        hline L1 (O1[0],O1[1],O1[2]);
             std::cout<<" "<<O1<<" "<<F01<<" "<<F02<< " "<<F03<< " "<<F04<< " "<<F1<<" "<<F2<< " "<<F3<< " "<<F4<<std::endl;
            if(O1[0]<=max(A1[0],B1[0]) && O1[2]>=min(A1[2],B1[2]) && O1[2]<=max(A1[2],B1[2]))
            {
                d=fabs(O1[1]-y);
                 std::cout<<" "<<d<<std::endl;
                return d;

            }
            else
            {
                d1= min(L1.distance2(F01, F02),L1.distance2(F02, F03));
                d2=min(L1.distance2(F03, F04),L1.distance2(F04, F01));
                d=min(d1,d2);
                 std::cout<<" "<<L1.distance2(F01, F02)<<" "<<L1.distance2(F02, F03)<< " "<<L1.distance2(F03, F04)<< " "<<L1.distance2(F04, F01)<<std::endl;
                return d;
            }


        }

       if (A1[2]==B1[2])

        {

            z=A1[2];
            F01 = mmx::point<double>(min(A1[0],B1[0]),min(A1[1],B1[1]),z);
            F02 = mmx::point<double>(min(A1[0],B1[0]),max(A1[1],B1[1]),z);
            F03 = mmx::point<double>(max(A1[0],B1[0]),max(A1[1],B1[1]),z);
            F04 = mmx::point<double>(max(A1[0],B1[0]),min(A1[1],B1[1]),z);
            F1  = F01[0]*U+F01[1]*V+F01[2]*W;
            F2  = F02[0]*U+F02[1]*V+F02[2]*W;
            F3  = F03[0]*U+F03[1]*V+F03[2]*W;
            F4  = F04[0]*U+F04[1]*V+F04[2]*W;
            std::cout<< " "<<F1<<" "<<F2<< " "<<F3<< " "<<F4<<std::endl;

            v0=F2-F1;
            w0=F4-F1;
            p0=v0.cross(w0);
            G=F1-O;
            e1=G.dot(v0);
            e2=G.dot(w0);
            e=v0.dot(w0);
            v1=v0.dot(v0);
            w1=w0.dot(w0);
            F=B-A;
            R=O-A;
            s=(e1*w1-e2*e)/(v1*w1-e*e);
            t=(e2*v1-e1*e)/(v1*w1-e*e);
            hline L1 (O[0],O[1],O[2]);
            if(u.dot(F)==0)
            {


                if (s>=0 && s<=1 && t>=0 && t<=1)
                {
                    P=F1+s*v0+t*w0;
                    if((P[0]-O[0])/u[0]<0)
                    {
                        Q=P-O;
                        d=sqrt(Q.dot(Q));
                        return d;
                    }
                    else
                    {
                        d=0;
                        return d;
                    }
                }
                else
                {
                    d1= min(L1.distance2(F1, F2),L1.distance2(F2, F3));
                    d2= min(L1.distance2(F3, F4),L1.distance2(F4, F1));
                    d=min(d1,d2);
                    std::cout<< " "<<d1<<" "<<d2<<std::endl;
                    std::cout<< " "<<L1.distance2(F1, F2)<<" "<<L1.distance2(F2, F3)<< " "<<L1.distance2(F3, F4)<< " "<<L1.distance2(F4, F1)<<std::endl;
                    return d;
                }


            }

            if(u.cross(F)==0)
            {

                if (s>=0 && s<=1 && t>=0 && t<=1)
                {
                    if(p0.dot(R)!=0)
                    {
                        P=F1+s*v0+t*w0;
                        Q=P-O;
                        d=sqrt(Q.dot(Q));
                        return d;
                    }

                    else
                    {
                        d=0;
                        return d;
                    }
                }
                else
                {
                    d1= min(L1.distance2(F1, F2),L1.distance2(F2, F3));
                    d2= min(L1.distance2(F3, F4),L1.distance2(F4, F1));
                    d=min(d1,d2);
                    std::cout<< " "<<d1<<" "<<d2<<std::endl;
                    std::cout<< " "<<L1.distance2(F1, F2)<<" "<<L1.distance2(F2, F3)<< " "<<L1.distance2(F3, F4)<< " "<<L1.distance2(F4, F1)<<std::endl;
                    return d;

                }




            }

        }


         hline L1 (O1[0],O1[1],O1[2]);
             std::cout<<" "<<O1<<" "<<F01<<" "<<F02<< " "<<F03<< " "<<F04<< " "<<F1<<" "<<F2<< " "<<F3<< " "<<F4<<std::endl;
            if(O1[0]<=max(A1[0],B1[0]) && O1[1]>=min(A1[1],B1[1]) && O1[1]<=max(A1[1],B1[1]))
            {
                d = fabs(O1[2]-z);
                return d;

            }
            else
            {
                d1= min(L1.distance2(F01, F02),L1.distance2(F02, F03));
                d2= min(L1.distance2(F03, F04),L1.distance2(F04, F01));
                d=min(d1,d2);
                std::cout<<" "<<L1.distance2(F01, F02)<<" "<<L1.distance2(F02, F03)<< " "<<L1.distance2(F03, F04)<< " "<<L1.distance2(F04, F01)<<std::endl;
                return d;
            }


      }


        // The projection of O on the face [AB]


    //}



//}



    if( b!=0 || c != 0)
    {
        a1 = sqrt(a * a + b * b + c * c);
        b1= sqrt(b * b + c * c);
        c1 = sqrt(a * a * b * b + (b * b + c * c)*(b * b + c * c) + a * a * c * c);
        U= mmx::point<double>(a / a1, b / a1,  c / a1);
        V= mmx::point<double>(0,-c/b1, b / b1);
        W= mmx::point<double>((b * b + c * c) / c1,-a * b / c1,  -a * c / c1);
        A1= mmx::point<double> ( a * xa / a1 - (-a1 * (b*b*b) - a1 * b * c * c) * ya / (c1*c1) + (a1 * b * b * c + a1 * (c*c*c)) * za / (c1*c1), -(a * a * b1 * c + b * b * b1 * c + b1 * (c*c*c)) * ya / (c1*c1) + (a * a * b * b1 + (b*b*b) * b1 + b * b1 * c * c) * za /(c1*c1), c1 * xa / (a1*a1) - a * b * ya / c1 - a * c * za / c1);
        B1= mmx::point<double> ( a * xb / a1 - (-a1 * (b*b*b) - a1 * b * c * c) * yb / (c1*c1) + (a1 * b * b * c + a1 * (c*c*c)) * zb / (c1*c1), -(a * a * b1 * c + b * b * b1 * c + b1 * (c*c*c)) * yb / (c1*c1) + (a * a * b * b1 + (b*b*b) * b1 + b * b1 * c * c) * zb /(c1*c1), c1 * xb / (a1*a1) - a * b * yb / c1 - a * c * zb / c1);
        O1= mmx::point<double> ( a * xo / a1 - (-a1 * (b*b*b) - a1 * b * c * c) * yo / (c1*c1) + (a1 * b * b * c + a1 * (c*c*c)) * zo / (c1*c1), -(a * a * b1 * c + b * b * b1 * c + b1 * (c*c*c)) * yo / (c1*c1) + (a * a * b * b1 + (b*b*b) * b1 + b * b1 * c * c) * zo /(c1*c1), c1 * xo / (a1*a1) - a * b * yo / c1 - a * c * zo / c1);

        if (A1[0]==B1[0])

        {
            x=A1[0];
            F01 = mmx::point<double>(x,min(A1[1],B1[1]),min(A1[2],B1[2]));
            F02 = mmx::point<double>(x,min(A1[1],B1[1]),max(A1[2],B1[2]));
            F03 = mmx::point<double>(x,max(A1[1],B1[1]),max(A1[2],B1[2]));
            F04 = mmx::point<double>(x,max(A1[1],B1[1]),min(A1[2],B1[2]));
            F1  = F01[0]*U+F01[1]*V+F01[2]*W;
            F2  = F02[0]*U+F02[1]*V+F02[2]*W;
            F3  = F03[0]*U+F03[1]*V+F03[2]*W;
            F4  = F04[0]*U+F04[1]*V+F04[2]*W;
            std::cout<< " "<<F1<<" "<<F2<< " "<<F3<< " "<<F4<<std::endl;

            v0=F2-F1;
            w0=F4-F1;
            p0=v0.cross(w0);
            G=F1-O;
            e1=G.dot(v0);
            e2=G.dot(w0);
            e=v0.dot(w0);
            v1=v0.dot(v0);
            w1=w0.dot(w0);
            F=B-A;
            R=O-A;
            s=(e1*w1-e2*e)/(v1*w1-e*e);
            t=(e2*v1-e1*e)/(v1*w1-e*e);
            hline L1 (O[0],O[1],O[2]);
            if(u.dot(F)==0)
            {


                if (s>=0 && s<=1 && t>=0 && t<=1)
                {
                    P=F1+s*v0+t*w0;
                    if((P[0]-O[0])/u[0]<0)
                    {
                        Q=P-O;
                        d=sqrt(Q.dot(Q));
                        return d;
                    }
                    else
                    {
                        d=0;
                        return d;
                    }
                }
                else
                {
                    d1= min(L1.distance2(F1, F2),L1.distance2(F2, F3));
                    d2= min(L1.distance2(F3, F4),L1.distance2(F4, F1));
                    d=min(d1,d2);
                    std::cout<< " "<<d1<<" "<<d2<<std::endl;
                    std::cout<< " "<<L1.distance2(F1, F2)<<" "<<L1.distance2(F2, F3)<< " "<<L1.distance2(F3, F4)<< " "<<L1.distance2(F4, F1)<<std::endl;
                    return d;
                }


            }

            if(u.cross(F)==0)
            {

                if (s>=0 && s<=1 && t>=0 && t<=1)
                {
                    if(p0.dot(R)!=0)
                    {
                        P=F1+s*v0+t*w0;
                        Q=P-O;
                        d=sqrt(Q.dot(Q));
                        return d;
                    }

                    else
                    {
                        d=0;
                        return d;
                    }
                }
                else
                {
                    d1= min(L1.distance2(F1, F2),L1.distance2(F2, F3));
                    d2= min(L1.distance2(F3, F4),L1.distance2(F4, F1));
                    d=min(d1,d2);
                    std::cout<< " "<<d1<<" "<<d2<<std::endl;
                    std::cout<< " "<<L1.distance2(F1, F2)<<" "<<L1.distance2(F2, F3)<< " "<<L1.distance2(F3, F4)<< " "<<L1.distance2(F4, F1)<<std::endl;
                    return d;

                }




            }
        }
        hline L1 (O1[0],O1[1],O1[2]);
             std::cout<<" "<<O1<<" "<<F01<<" "<<F02<< " "<<F03<< " "<<F04<< " "<<F1<<" "<<F2<< " "<<F3<< " "<<F4<<std::endl;
            if(O1[1]>=min(A1[1],B1[1]) && O1[1]<=max(A1[1],B1[1]) && O1[2]>=min(A1[2],B1[2]) && O1[2]<=max(A1[2],B1[2]))
            {
                if (O1[0]> x)
                {d= fabs(O1[0]-x);
                    return d;
                }
                else if (O1[0]<= x)
                {
                    d=0;
                    return d;
                }
            }
            else
            {
                d1= min(L1.distance2(F01, F02),L1.distance2(F02, F03));
                d2= min(L1.distance2(F03, F04),L1.distance2(F04, F01));
                  std::cout<<" "<<d1<<" "<<d2<<std::endl;
                d=min(d1,d2);
                return d;
            }
        }

      if (A1[1]==B1[1])

        {
            y=A1[1];
            F01 = mmx::point<double>(min(A1[0],B1[0]),y,min(A1[2],B1[2]));
            F02 = mmx::point<double>(min(A1[0],B1[0]),y,max(A1[2],B1[2]));
            F03 = mmx::point<double>(max(A1[0],B1[0]),y,max(A1[2],B1[2]));
            F04 = mmx::point<double>(max(A1[0],B1[0]),y,min(A1[2],B1[2]));
            F1  = F01[0]*U+F01[1]*V+F01[2]*W;
            F2  = F02[0]*U+F02[1]*V+F02[2]*W;
            F3  = F03[0]*U+F03[1]*V+F03[2]*W;
            F4  = F04[0]*U+F04[1]*V+F04[2]*W;
            std::cout<< " "<<F1<<" "<<F2<< " "<<F3<< " "<<F4<<std::endl;



            v0=F2-F1;
            w0=F4-F1;
            p0=v0.cross(w0);
            G=F1-O;
            e1=G.dot(v0);
            e2=G.dot(w0);
            e=v0.dot(w0);
            v1=v0.dot(v0);
            w1=w0.dot(w0);
            F=B-A;
            R=O-A;
            s=(e1*w1-e2*e)/(v1*w1-e*e);
            t=(e2*v1-e1*e)/(v1*w1-e*e);
            hline L1 (O[0],O[1],O[2]);
            if(u.dot(F)==0)
            {


                if (s>=0 && s<=1 && t>=0 && t<=1)
                {
                    P=F1+s*v0+t*w0;
                    if((P[0]-O[0])/u[0]<0)
                    {
                        Q=P-O;
                        d=sqrt(Q.dot(Q));
                        return d;
                    }
                    else
                    {
                        d=0;
                        return d;
                    }
                }
                else
                {
                    d1= min(L1.distance2(F1, F2),L1.distance2(F2, F3));
                    d2= min(L1.distance2(F3, F4),L1.distance2(F4, F1));
                    d=min(d1,d2);
                    std::cout<< " "<<d1<<" "<<d2<<std::endl;
                    std::cout<< " "<<L1.distance2(F1, F2)<<" "<<L1.distance2(F2, F3)<< " "<<L1.distance2(F3, F4)<< " "<<L1.distance2(F4, F1)<<std::endl;
                    return d;
                }


            }

            if(u.cross(F)==0)
            {

                if (s>=0 && s<=1 && t>=0 && t<=1)
                {
                    if(p0.dot(R)!=0)
                    {
                        P=F1+s*v0+t*w0;
                        Q=P-O;
                        d=sqrt(Q.dot(Q));
                        return d;
                    }

                    else
                    {
                        d=0;
                        return d;
                    }
                }
                else
                {
                    d1= min(L1.distance2(F1, F2),L1.distance2(F2, F3));
                    d2= min(L1.distance2(F3, F4),L1.distance2(F4, F1));
                    d=min(d1,d2);
                    std::cout<< " "<<d1<<" "<<d2<<std::endl;
                    std::cout<< " "<<L1.distance2(F1, F2)<<" "<<L1.distance2(F2, F3)<< " "<<L1.distance2(F3, F4)<< " "<<L1.distance2(F4, F1)<<std::endl;
                    return d;

                }




            }

        }


 hline L1 (O1[0],O1[1],O1[2]);
             std::cout<<" "<<O1<<" "<<F01<<" "<<F02<< " "<<F03<< " "<<F04<< " "<<F1<<" "<<F2<< " "<<F3<< " "<<F4<<std::endl;
            if(O1[0]<=max(A1[0],B1[0]) && O1[2]>=min(A1[2],B1[2]) && O1[2]<=max(A1[2],B1[2]))
            {
                d= fabs(O1[1]-y);
                return d;

            }
            else
            {
                d1= min(L1.distance2(F01, F02),L1.distance2(F02, F03));
                d2=min(L1.distance2(F03, F04),L1.distance2(F04, F01));
                d=min(d1,d2);
                  std::cout<<" "<<d1<<" "<<d2<<std::endl;
                return d;
            }


        }
        if (A1[2]==B1[2])

        {

            z=A1[2];
            F01 = mmx::point<double>(min(A1[0],B1[0]),min(A1[1],B1[1]),z);
            F02 = mmx::point<double>(min(A1[0],B1[0]),max(A1[1],B1[1]),z);
            F03 = mmx::point<double>(max(A1[0],B1[0]),max(A1[1],B1[1]),z);
            F04 = mmx::point<double>(max(A1[0],B1[0]),min(A1[1],B1[1]),z);
            F1  = F01[0]*U+F01[1]*V+F01[2]*W;
            F2  = F02[0]*U+F02[1]*V+F02[2]*W;
            F3  = F03[0]*U+F03[1]*V+F03[2]*W;
            F4  = F04[0]*U+F04[1]*V+F04[2]*W;
            std::cout<< " "<<F1<<" "<<F2<< " "<<F3<< " "<<F4<<std::endl;

            v0=F2-F1;
            w0=F4-F1;
            p0=v0.cross(w0);
            G=F1-O;
            e1=G.dot(v0);
            e2=G.dot(w0);
            e=v0.dot(w0);
            v1=v0.dot(v0);
            w1=w0.dot(w0);
            F=B-A;
            R=O-A;
            s=(e1*w1-e2*e)/(v1*w1-e*e);
            t=(e2*v1-e1*e)/(v1*w1-e*e);
            hline L1 (O[0],O[1],O[2]);
            if(u.dot(F)==0)
            {


                if (s>=0 && s<=1 && t>=0 && t<=1)
                {
                    P=F1+s*v0+t*w0;
                    if((P[0]-O[0])/u[0]<0)
                    {
                        Q=P-O;
                        d=sqrt(Q.dot(Q));
                        return d;
                    }
                    else
                    {
                        d=0;
                        return d;
                    }
                }
                else
                {
                    d1= min(L1.distance2(F1, F2),L1.distance2(F2, F3));
                    d2= min(L1.distance2(F3, F4),L1.distance2(F4, F1));
                    d=min(d1,d2);
                    std::cout<< " "<<d1<<" "<<d2<<std::endl;
                    std::cout<< " "<<L1.distance2(F1, F2)<<" "<<L1.distance2(F2, F3)<< " "<<L1.distance2(F3, F4)<< " "<<L1.distance2(F4, F1)<<std::endl;
                    return d;
                }


            }

            if(u.cross(F)==0)
            {

                if (s>=0 && s<=1 && t>=0 && t<=1)
                {
                    if(p0.dot(R)!=0)
                    {
                        P=F1+s*v0+t*w0;
                        Q=P-O;
                        d=sqrt(Q.dot(Q));
                        return d;
                    }

                    else
                    {
                        d=0;
                        return d;
                    }
                }
                else
                {
                    d1= min(L1.distance2(F1, F2),L1.distance2(F2, F3));
                    d2= min(L1.distance2(F3, F4),L1.distance2(F4, F1));
                    d=min(d1,d2);
                    std::cout<< " "<<d1<<" "<<d2<<std::endl;
                    std::cout<< " "<<L1.distance2(F1, F2)<<" "<<L1.distance2(F2, F3)<< " "<<L1.distance2(F3, F4)<< " "<<L1.distance2(F4, F1)<<std::endl;
                    return d;

                }




            }



        }

         hline L1 (O1[0],O1[1],O1[2]);
            std::cout<<" "<<O1<<" "<<F01<<" "<<F02<< " "<<F03<< " "<<F04<< " "<<F1<<" "<<F2<< " "<<F3<< " "<<F4<<std::endl;
            if(O1[0]<=max(A1[0],B1[0]) && O1[1]>=min(A1[1],B1[1]) && O1[1]<=max(A1[1],B1[1]))
            {
                d= fabs(O1[2]-z);
                return d;

            }
            else
            {
                d1= min(L1.distance2(F01, F02),L1.distance2(F02, F03));
                d2= min(L1.distance2(F03, F04),L1.distance2(F04, F01));
                d=min(d1,d2);
                  std::cout<<" "<<d1<<" "<<d2<<std::endl;
                return d;
            }





    }



}*/
