/**********************************************************************
 * PACKAGE  : geometrix 
 * COPYRIGHT: (C) 2015, Bernard Mourrain, Inria
 **********************************************************************/
#ifndef TUPLE_HPP
#define TUPLE_HPP

//====================================================================
namespace mmx {
//--------------------------------------------------------------------

template<int N, typename T>
struct tuple {

    tuple(int a, int b) {
        d_tab[0]=a;
        d_tab[1]=b;
    }

    tuple(int a, int b, int c) {
        d_tab[0]=a;
        d_tab[1]=b;
        d_tab[2]=c;
    }

    tuple(int a, int b, int c, int d) {
        d_tab[0]=a;
        d_tab[1]=b;
        d_tab[2]=c;
        d_tab[3]=d;
    }
    int operator[] (int i) const { return d_tab[i]; }

    int d_tab[N];

    static int size;
};

template<int N, typename T>
int tuple<N,T>::size =N;


//--------------------------------------------------------------------
} /* namespace mmx */
//====================================================================
#undef TMPL
#endif //REALGEO_TUPLE_HPP
