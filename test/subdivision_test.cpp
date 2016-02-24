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
        hline(-1,2,-1),
        hline(2,-3,3)
    };
    int N = 4;

    distfield<hline>* f = new distfield<hline>;
    for(unsigned i=0;i<N;i++) f->add(L[i]);


    int t;
    cout<<"distance^2: "<<f->distance2(-2,-3,0, t)<<" #"<<t<<" "<< f->site(t)<<endl;

    voronoi_hline<double>::Controler* ctrl = new voronoi_hline<double>::Controler;

    mmx::subdivider< voronoi_hline<double>::Controler > mshr(ctrl);
    mshr.set_max_size(1);
    mshr.set_min_size(0.1);
    mshr.set_input(mshr.controler()->init_cell(f, -6, 6, -6, 5, -5, 5));
    mshr.run();

    mdebug()<<"Subdivision"
           <<mshr.controler()->m_regular.size()
           <<mshr.controler()->m_singular.size();

    axldebug f0("voronoi_mesh.axl", axldebug::init);
    for(unsigned i=0;i<N;i++) axl_print(f0,L[i],5);

    voronoi_hline<double>::Mesh* mt= new voronoi_hline<double>::Mesh;
    mshr.get_tmsh(mt);
    f0.set_color(0,0,255);
    f0<< *mt;
    f0.close();


    mmx::polygonizer3d_mc< voronoi_hline<double>::Controler > mc(mshr.controler());
    mc.set_regular(&mshr.controler()->m_regular);
    mc.run();

    voronoi_hline<double>::Mesh* m= new voronoi_hline<double>::Mesh;
    mc.get_mesh(m);
    mdebug()<<"Mesh:"<<m->nbv()<<m->nbf();
    axldebug a1("voronoi_surface.axl",axldebug::init);
    a1<<*m;
    a1.close();

    system("axel voronoi_mesh.axl voronoi_surface.axl &");

}

