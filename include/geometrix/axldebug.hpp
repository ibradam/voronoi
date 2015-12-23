/*********************************************************************
 * geometrix package
 *********************************************************************
 * Bernard Mourrain, Inria
 *********************************************************************/
#ifndef AXLDEBUG_HPP
#define AXLDEBUG_HPP
#include <fstream>
#include <string>
#include <geometrix/mesh.hpp>

#ifndef AXL_DEBUG_PREFIX
#define AXL_DEBUG_PREFIX ""
#endif

struct axldebug {

    enum mode {init, app};

    axldebug( mode md = app) {
        m_file = std::string(AXL_DEBUG_PREFIX)+"tmp.axl";

        if (md == init) {
            m_f.open (m_file.c_str(), std::fstream::out);
            m_f << "<axl>\n";
        } else {
            m_f.open (m_file.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
        }
    }

    axldebug(const char* file, mode md = app) {
        m_file = std::string(AXL_DEBUG_PREFIX)+file;

        if (md == init) {
            m_f.open (m_file.c_str(), std::fstream::out);
            m_f << "<axl>\n";
        } else {
            m_f.open (m_file.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
        }
    }

    ~axldebug() {
        m_f<<"\n";
        m_f.close();
    }

     void close() {
        m_f << "</axl>\n";
        m_f.close();
    }

    int view(void) {
        m_f<< "</axl>\n";
        m_f.close();
        std::string cmd = std::string("axel ")+m_file+" &";
        return system(cmd.c_str());
    }

    axldebug& operator << (const char* x)  { m_f<< x; return *this; }
    axldebug& operator << (int x)          { m_f<< x; return *this; }
    axldebug& operator << (unsigned x)     { m_f<< x; return *this; }
    axldebug& operator << (double x)       { m_f<< x; return *this; }

    template<class C, class V>
    axldebug& operator<< (const mmx::mesh<C,V>& msh) {
        mmx::axl_print(*this,msh,mmx::color(255,0,0));
        return *this;
    }
    std::fstream m_f;
    std::string  m_file;

};



#endif // AXLDEBUG_HPP
