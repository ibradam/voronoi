#ifndef HLINE_H
#define HLINE_H
#include <iostream>
#include <voronoi/point.hpp>

struct hline {

    hline(double X,double Y, double Z): m_pt(X,Y,Z), m_dir(0,0,1) {}
    hline(double X,double Y, double Z, double U, double V, double W): m_pt(X,Y,Z),m_dir(U,V,W){
        m_dir /= m_dir.norm();
    }

    double distance2(double X,double Y, double Z) const;
    double distance2(const mmx::point<double>& p) const;

    mmx::point<double> m_pt;
    mmx::point<double> m_dir;

};

inline void print(const hline& L) {
    std::cout<<"hline ("<<L.m_pt[0]<<" "<<L.m_pt[1] <<" "<<L.m_pt[2]<<" -> "<<L.m_dir[0]<<" "<<L.m_dir[1] <<" "<<L.m_dir[2]<<")"<<std::endl;
}

inline std::ostream& operator<<( std::ostream& os , const hline& L) {
    os<<"hline ("<<L.m_pt[0]<<" "<<L.m_pt[1] <<" "<<L.m_pt[2]<<" -> "<<L.m_dir[0]<<" "<<L.m_dir[1] <<" "<<L.m_dir[2]<<")";
    return os;
}


#endif // HLINE_H
