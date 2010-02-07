#ifndef FORTRAN_H
#define FORTRAN_H

#if defined(NOF77UNDERSCORE)
#define F77(s) s
#elif defined (F77UPPERCASE)
#define F77(s) s
#define dgemv DGEMV
#define dgetrf DGETRF
#define dgetrs DGETRS
#define dscal DSCAL
#define ilaenv ILAENV
#define dsyev DSYEV
#define dgemm DGEMM
#define daxpy DAXPY
#define dnrm2 DNRM2
#define dcopy DCOPY
#define ddot DDOT
#define dlaset DLASET
#define dlacpy DLACPY
#define dlarnv DLARNV
#else
#define F77(s) s##_
#endif

#endif

