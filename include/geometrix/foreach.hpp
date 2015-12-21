/**********************************************************************
 * PACKAGE  : geometrix 
 * COPYRIGHT: (C) 2015, Bernard Mourrain, Inria
 **********************************************************************/
#ifndef MMX_FOREACH
 
struct MmxForeachContainerBase {};

template <typename T> class MmxForeachContainer : public MmxForeachContainerBase {
public:
    inline MmxForeachContainer(const T& t): c(t), brk(0), i(c.begin()), e(c.end()) {} ;
    const T c ;
    mutable int brk ;
    mutable typename T::const_iterator i, e ;
    inline bool condition() const { 
        return (!brk++ && i != e) ; 
    }
} ;

template <typename T>
inline T * mmxForeachpointer(const T &) {
    return 0 ;
}

template <typename T>
inline MmxForeachContainer<T>
mmxForeachContainerNew(const T& t) {
    return MmxForeachContainer<T>(t) ;
}

template <typename T> inline const MmxForeachContainer<T> *mmxForeachContainer(const MmxForeachContainerBase *base, const T *) {
    return static_cast<const MmxForeachContainer<T> *>(base) ;
}

# define MMX_FOREACH(variable, container) \
     if(0) {} else for (const MmxForeachContainerBase & _container_ = mmxForeachContainerNew(container); \
          mmxForeachContainer(&_container_, true ? 0 : mmxForeachpointer(container))->condition(); \
          ++mmxForeachContainer(&_container_, true ? 0 : mmxForeachpointer(container))->i) \
         for (variable = *mmxForeachContainer(&_container_, true ? 0 : mmxForeachpointer(container))->i; \
              mmxForeachContainer(&_container_, true ? 0 : mmxForeachpointer(container))->brk; \
              --mmxForeachContainer(&_container_, true ? 0 : mmxForeachpointer(container))->brk)

# endif

#ifndef foreach
# define foreach MMX_FOREACH
#endif