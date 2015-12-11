
#include <voronoi/hline.h>
#include <voronoi/point.hpp>
#include <voronoi/mediatrice.h>

int main()
{
    hline L1(2,5,-1);
    print (L1);
    hline L2(-2,4,1);
    print (L2);
    std::cout<<"distance "<< L2.distance2(4,-4,4)<<std::endl;
    mmx::point<double> p(2,2,5);
    std::cout<<"distance "<< L1.distance2(p)<<std::endl;

    mmx::point< double> A(6,6,-3), B(-6,3,3);
    std::cout<<A<<" "<<B<<std::endl;
    std::cout<<"point-equidistant "<< equidist(L1,L2,A,B) <<std::endl;


}
