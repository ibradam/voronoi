#ifndef MEDIATRICE_H
#define MEDIATRICE_H

mmx::point<double> equidist(const hline& H1 , const hline& H2,
                            const mmx::point<double>& A,
                            const mmx::point<double>& B);

mmx::point<double> equidist(const hline& H1 , const hline& H2, const hline& H3,
                            const mmx::point<double>& A, const mmx::point<double>& B,
                            int& info);


#endif // MEDIATRICE_H

