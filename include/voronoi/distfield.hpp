/**********************************************************************
 * PACKAGE  : voronoi
 * COPYRIGHT: (C) 2015, Ibrahim Adamou, Bernard Mourrain, Inria
 **********************************************************************/
#pragma once

#include <vector>

#define TMPL template<class SITE>
#define SELF distfield<SITE>

//====================================================================
template<class SITE>
struct distfield {

    void add(const SITE& L) { m_sites.push_back(L) ;}

    double distance2(double x, double y, double z, int& i);

    unsigned nbs() const { return m_sites.size(); }

    SITE  site(unsigned i) const { return m_sites[i]; }
    SITE& site(unsigned i)       { return m_sites[i]; }

private:
    std::vector<SITE> m_sites;
};

//--------------------------------------------------------------------
TMPL
double SELF::distance2(double x, double y, double z, int& i) {
    double r = m_sites[0].distance2(x,y,z), d;
    i = 0;
    for(unsigned k=1;k<this->nbs();k++) {
        d = m_sites[k].distance2(x,y,z);
        if(d<r) {
            r = d;
            i=k;
        }
    }
    //mdebug()<<"closest"<<i;
    return r;
}
//====================================================================
#undef SELF
#undef TMPL
