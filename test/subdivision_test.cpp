#include <geometrix/axldebug.hpp>
#include <voronoi/hline.h>
#include <voronoi/point.hpp>
#include <voronoi/mediatrice.h>
#include <voronoi/voronoi_hline.hpp>
#include <voronoi/distfield.hpp>

int main()
{
    using std::cout;
    using std::endl;
    hline L[4] = {
        hline(1,2,-3),
        hline(1,-2,1),
        hline(-0.9,2.1,-1),
        hline(2,-3,3)
    };
    int N = 4;

    distfield<hline>* f = new distfield<hline>;
    for(unsigned i=0;i<N;i++) f->add(L[i]);

    int t;
    cout<<"distance^2: "<<f->distance2(-2,-3,0, t)<<" #"<<t<<" "<< f->site(t)<<endl;

    voronoi_hline<double>::Controler* ctrl = new voronoi_hline<double>::Controler;

    mmx::subdivision< voronoi_hline<double>::Controler > mshr(ctrl);
    mshr.set_max_size(1);
    mshr.set_min_size(0.1);
    mshr.set_input(mshr.controler()->init_cell(f, -5, 5, -5, 5, -5, 5));
    mshr.run();

    //mdebug()<<"Mesh:"<<mshr.output()->nbv()<<mshr.output()->nbe()<<mshr.output()->nbf();
    mdebug()<<"Subdivision"
           <<mshr.controler()->m_regular.size()
           <<mshr.controler()->m_singular.size();


    axldebug f0("voronoi_mesh.axl", axldebug::init);
    for(unsigned i=0;i<N;i++) axl_print(f0,L[i],5);
    //f0<<*mshr.controler();
    f0<< *ctrl;
    f0.close();

    //axldebug f1("voronoi_surface.axl", axldebug::init);
    //f1<<*mshr.output();
    //f1.close();
    //std::cout<<"axel voronoi_mesh.axl voronoi_surface.axl "<<std::endl;
    system("axel voronoi_mesh.axl&");

}

