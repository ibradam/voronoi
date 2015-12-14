
#include <voronoi/hline.h>
#include <voronoi/point.hpp>
#include <voronoi/mediatrice.h>

int main()
{
    hline L1(-1,1,-3);
    print (L1);
    hline L2(1,-1,-1);
    print (L2);
     hline L3(2,1,1);
    std::cout<<"distance "<< L2.distance2(4,-4,4)<<std::endl;
    mmx::point<double> p(2,2,5);
    std::cout<<"distance "<< L1.distance2(p)<<std::endl;

    mmx::point< double> A(-0.22,0.94,3.14), B(-1.34,-0.21,10);
    std::cout<<A<<" "<<B<<std::endl;
    std::cout<<"point-equidistant "<< equidist(L1,L2,A,B) <<std::endl;

    mmx::point< double> C(-1,-1,3), D(2,1,3);
    std::cout<<C<<" "<<D<<std::endl;
    int info;
    std::cout<<"point-equidistant "<< equidist(L1,L2,L3,C,D,info) <<std::endl;
}
