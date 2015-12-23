/**********************************************************************
 * PACKAGE  : geometrix
 * COPYRIGHT: (C) 2015, Bernard Mourrain, Inria
 **********************************************************************/
#pragma once

#include <vector>
#include <geometrix/mdebug.hpp>
#include <geometrix/tmsh_cell.hpp>
#include <geometrix/tmsh_vertex.hpp>

//====================================================================
namespace mmx {
//--------------------------------------------------------------------
template <class M, class P>
void set_middle(M& m, const P& a, const P& b) {
    for(unsigned i=0;i<M::dim;i++) m[i]=(a[i]+b[i])/2.;
}
//--------------------------------------------------------------------
} /* namespace mmx */
//====================================================================

#include <geometrix/tmsh_controler2.hpp>
#include <geometrix/tmsh_controler3.hpp>

//--------------------------------------------------------------------
#undef SELF
#undef TMPL

