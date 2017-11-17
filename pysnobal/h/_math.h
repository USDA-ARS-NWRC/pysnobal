#ifndef	IPW_MATH_H
#define	IPW_MATH_H

/* ------------------------------------------------------------------------ */

/*
 *  typedef for complex numbers (from UCSB's header file: complex.h)
 */
typedef struct {
        double  re;
        double  im;
} COMPLEX_T;
 
/*
 * degrees in circle
 */
#define DEGS_IN_CIRCLE  3.6e2

/* ------------------------------------------------------------------------ */

/*
 *  Macros
 */

#define	BIT(i)			( 1 << (i) )
/*
 *  Returns 2 raised to the i'th power (which is a bit mask for i'th bit)
 */

#define	_MASK(n)		( ~(~0 << (n)) )
#define	MASK(n)			( (_MASK((n) - 1) << 1) | 1 )
/*
 *  Returns a mask of the lower n bits.
 */

#define ABS(x)			( (x) < 0 ? -(x) : (x) )
#define MAX(a, b)               ( (a) > (b) ? (a) : (b) )
#define	MIN(a, b)		( (a) < (b) ? (a) : (b) )
#define	ODD(n)			( (n) & 1 )

/*
 *  Angles, and arcs.
 */

#define DEG_TO_RAD(d)		( (d) * (M_PI / 180.) )
#define RAD_TO_DEG(r)		( (r) * (180. / M_PI) )
#define MIN_TO_DEG(m)		( (m) / 60. )
#define SEC_TO_DEG(s)		( (s) / 3600. )

/* ------------------------------------------------------------------------ */

extern int	akcoef(double * x, double * y, int nx, double * c);
extern void	apfit(double *arx, double *ary, long n, long d, double *c,
		      const int savemem);
extern char   * dtoa(char * s, double d);
extern int	hbit(unsigned i);
extern int	ipow2(int exponent);
extern char   * itoa(char * s, int i);
extern int	ltoi(long i);
extern void	msolve(double **A, double *x, double *b, long n);
extern int	ndig(int i);
extern int	rotate(double mu, double azm, double mu_r, double lam_r,
		       double * muPrime, double * aPrime);
extern int	seval(double x[], double y[], int nx, double c[],
		      double  * u, double * s, int m);
extern double	splint(double * x, double * y, int nx, double * c,
		       double a, double b);
extern double	zerobr(double a, double b, double t, double (*f)() );


#endif  /* IPW_MATH_H */
