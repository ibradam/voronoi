/**********************************************************************
 * PACKAGE  : mdebug for easy print of debug  message 
 * COPYRIGHT: (C) 2015, Bernard Mourrain, Inria
 **********************************************************************/
#ifndef MDEBUG_HPP
#define MDEBUG_HPP
#include <iostream>
#include <fstream>
#include <string>

#ifndef MDEBUG_PREFIX
#define MDEBUG_PREFIX ""
#endif

struct mdebug {

    enum mode {init, app};
    mdebug(int l) {
#ifdef MMDEBUG_COUT
        m_level=l;
#endif
    }

    mdebug( mode md = app) {
#ifdef MDEBUG_COUT
        m_level=0;
#endif
#ifdef MDEBUG_FILE
        m_file = std::string(MDEBUG_PREFIX)+"debug.txt";

        if (md == init) {
            m_f.open (m_file.c_str(), std::fstream::out);
        } else {
            m_f.open (m_file.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
        }
#endif
    }

    mdebug(const char* file, mode md = app) {
#ifdef MDEBUG_COUT
        m_level=0;
#endif
#ifdef DEBUG_FILE
        m_file = std::string(MDEBUG_PREFIX)+file;

        if (md == init) {
            m_f.open (m_file.c_str(), std::fstream::out);
        } else {
            m_f.open (m_file.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
        }
#endif
    }

    ~mdebug() {
#ifdef MDEBUG_FILE
        if(m_level<=verbose) m_f<<"\n";
        m_f.close();
#endif
#ifdef MDEBUG_COUT
        //if(m_level <= verbose)
            std::cout<<std::endl;
#endif
    }

    void clear() {
#ifdef MDEBUG_FILE
        m_f.close();
        m_f.open (m_file.c_str(), std::fstream::out);
#endif
    }

    template<class X>
    mdebug& operator << (const X& x)     {
#ifdef MDEBUG_FILE
        if(m_level<=verbose) m_f<<x<<" ";
#endif
#ifdef MDEBUG_COUT
        //if(m_level<=verbose)
            std::cout<<x<<" ";
#endif
        return *this;
    }
#ifdef MDEBUG_FILE
    std::fstream m_f;
    std::string  m_file;
#endif
#ifdef MDEBUG_COUT
    int m_level;
//    static int verbose;
#endif
};


#endif // MDEBUG_HPP
