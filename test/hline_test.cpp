
#include <voronoi/hline.h>
#include <voronoi/point.hpp>
#include <voronoi/mediatrice.h>

int main()
{
    hline L1(1,2,-3);
    print (L1);
    hline L2(2,1,2);
    print (L2);
    hline L3(-1,2,-1);
    print (L3);
    hline L4(2,-3,3);
    print (L4);
    std::cout<<"distance "<< L2.distance2(4,-4,4)<<std::endl;
    mmx::point<double> p(2,2,5);
    std::cout<<"distance "<< L1.distance2(p)<<std::endl;

    mmx::point< double> A(-3,-2,1), B(-3,2,3);
    std::cout<<A<<" "<<B<<std::endl;
    std::cout<<"point-equidistant "<< equidist(L1,L2,A,B) <<std::endl;

    mmx::point< double> C(0,-12.5,2), D(0,2.5,4);
    std::cout<<C<<" "<<D<<std::endl;
    int info;
    //std::cout<<"point-equidistant "<< equidist(L1,L2,L3,C,D,info) <<std::endl;
    // std::cout<<"point-equidistant "<< equidist(L1,L2,L3,L4,info) <<std::endl;
    mmx::point<double> q0=equidist(L1,L2,L3,C,D,info);
    std::cout<<"point-equidistant "<< q0 <<std::endl;
    std:: cout << L1.distance2(q0)<< std::endl;
    std:: cout << L2.distance2(q0)<< std::endl;
    std:: cout << L3.distance2(q0)<< std::endl;


    mmx::point<double> q1=equidist(L1,L2,A,B);
    std::cout<<"point-equidistant "<< q1 <<std::endl;
    std:: cout << L1.distance2(q1)<< std::endl;
    std:: cout << L2.distance2(q1)<< std::endl;


    mmx::point<double> q2=equidist(L1,L2,L3,C,D,info);
    std::cout<<"point-equidistant "<< q2 <<std::endl;
    std:: cout << L1.distance2(q2)<< std::endl;
    std:: cout << L2.distance2(q2)<< std::endl;
    std:: cout << L3.distance2(q2)<< std::endl;

    mmx::point<double> q=equidist(L1,L2,L3,L4,info);
    std::cout<<"point-equidistant "<< q <<std::endl;
    std:: cout << L1.distance2(q)<< std::endl;
    std:: cout << L2.distance2(q)<< std::endl;
    std:: cout << L3.distance2(q)<< std::endl;
    std:: cout << L4.distance2(q)<< std::endl;

     std::cout<<"distance "<< L2.distance2(A, B)<<std::endl;

     std:: cout << L2.distance2(A)<< std::endl;
     std:: cout << L2.distance2(B)<< std::endl;
 std::cout<<"distance "<< L2.distance2(L2, A, B)<<std::endl;
    return 0;

}
