#include <voronoi/hline.h>
#include <voronoi/point.hpp>
#include <voronoi/mediatrice.h>
#include <voronoi/voronoi_hline.hpp>
#include <voronoi/distfield.hpp>

int main()
{
    using std::cout;
    using std::endl;

    hline L1(1,2,-3);
    hline L2(1,-2,1);
    hline L3(-1,2,-1);
    hline L4(2,-3,3);

    distfield<hline>* f = new distfield<hline>;
    f->add(L1);
    f->add(L2);
    f->add(L3);
    f->add(L4);

    int i;
    cout<<"distance: "<<f->distance2(-2,-3,0, i)<<" #"<<i<<" "<< f->site(i)<<endl;

    mmx::mesher< voronoi_hline<double> > mshr;
    mshr.set_input(mshr.controler()->init_cell(f, -1, 1, -1, 1, -1, 1));
}

