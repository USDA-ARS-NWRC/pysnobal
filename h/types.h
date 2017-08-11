#ifndef	IPW_TYPES_H
#define	IPW_TYPES_H

/*
 * basic IPW types
 */

/* ------------------------------------------------------------------------- */

/*
 *  For backward compatibility, we define this type.
 *  New code should *not* use this type, but rather directly refer to
 *  "void *".
 */

typedef void   *addr_t;                 /* generic memory address	*/

/* ------------------------------------------------------------------------- */

/*
 *  Boolean type
 */

typedef int     bool_t;

#define	TRUE	1
#define	FALSE	0

/* ------------------------------------------------------------------------- */

/*
 *  Vectors of strings
 */

typedef struct {			/* vector of strings:		 */
	int             n;		/* -- # element in v		 */
	int             curr;		/* -- v[curr] is current string	 */
	char          **v;		/* -- -> strings		 */
} STRVEC_T;

extern STRVEC_T *addsv(STRVEC_T *strvec, char *string);
extern STRVEC_T *delsv(STRVEC_T *strvec, int i);
extern STRVEC_T *dupsv(STRVEC_T *strvec);
extern int       freesv(STRVEC_T *strvec);
extern char     *walksv(STRVEC_T *strvec, bool_t reset);

/* ------------------------------------------------------------------------- */

/*
 *  Integer and floating-point pixel types
 */

#define PIXEL	unsigned int
#define FPIXEL	float

typedef PIXEL   pixel_t;
typedef FPIXEL  fpixel_t;

#define PIXEL_BIT	(sizeof(PIXEL) * CHAR_BIT)

#define FPIXEL_MAX	FLT_MAX
#define FPIXEL_MIN	(-(FLT_MAX))

#endif  /* IPW_TYPES_H */
