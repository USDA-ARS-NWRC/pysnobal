#ifndef	_SORT_H
#define	_SORT_H

#ifndef SORT_ALG
/*
 * By default, use qsort.  On systems where qsort is slow, add _quicksort
 * to the Makefile in the $IPW/src/lib/libipw/util subdirectory, and 
 * issue the command "make install".
 * Then change SORT_ALG below to _quicksort
 */
#define SORT_ALG	qsort
#endif

/*
 *  GCC complains about the 4th argument to qsort of being incompatible
 *  pointer types.  Some how it doesn't consider a function prototype of
 *
 *		int compare_func(double *x, double *y);
 *
 *  to be compatible with qsort's 4th parameter, which is:
 *
 *		int (*) (const void *, const void *)
 *
 *  So we need to cast the 4th argument to silence the warning.
 */
#if SORT_ALG == qsort
#undef SORT_ALG
#define SORT_ALG(array,n,size,func)	qsort(array,n,size,(int (*)()) (func))

#elif SORT_ALG == _quicksort
#undef SORT_ALG
#define SORT_ALG(array,n,size,func) _quicksort(array,n,size,(int (*)()) (func))

#else
#error Unknown value for SORT_ALG macro
#endif

#endif  /* _SORT_H */
