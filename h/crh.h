#ifndef	CRH_H
#define	CRH_H

/*
 *	crh.h -- header file for class range (CRH) header
 */

typedef struct {
	int		cls;		/* class identifier		 */
	fpixel_t	lo;		/* low value in class		 */
	fpixel_t	hi;		/* high value in class		 */
	fpixel_t	rep;		/* representative value for class*/
} CLASS;

typedef struct {
	int             nclass;		/* # classes       		 */
	char	       *units;		/* -> units of classified data   */
	pixel_t		floor;		/* floor class (< lowest class)  */
	pixel_t		ceil;		/* floor class (> highest class) */
	char           *annot;		/* annotation      		 */
	CLASS	       *class;		/* class = class[class #]	 */
} CRH_T;

/* ------------------------------------------------------------------------ */
 
#define	CRH_HNAME	"cr"		/* header name within IPW	 */
#define	CRH_VERSION	"$Revision: 1.1 $"	/* SCCS revsion #	 */

/* field keywords */
#define	CRH_NCLASS	"nclass"
#define	CRH_ANNOT	"annot"
#define CRH_UNITS	"units"
#define	CRH_CLASS	"class"
#define CRH_FLOOR	"floor"
#define CRH_CEIL	"ceil"

/* field access macros */
#define	crh_nclass(p)	( (p)->nclass )
#define	crh_class(p)	( (p)->class )
#define	crh_lo(p,c)	( (p)->class[c].lo )
#define	crh_hi(p,c)	( (p)->class[c].hi )
#define	crh_rep(p,c)	( (p)->class[c].rep )
#define	crh_cls(p,c)	( (p)->class[c].cls )
#define	crh_floor(p)	( (p)->floor )
#define	crh_ceil(p)	( (p)->ceil )
#define	crh_annot(p)	( (p)->annot )
#define crh_units(p)	( (p)->units )

/* ------------------------------------------------------------------------ */
 
/* function declarations */

/* sort keys for crhsort() */
#define SORT_BY_CLASS	0
#define SORT_BY_RANGE	1

extern CRH_T  ** crh(int fd);
extern CRH_T  ** crhdup(CRH_T **oldhpp, int nbands);
extern CRH_T   * crhmake(int nclass, fpixel_t *lo, fpixel_t *hi,
			 fpixel_t *rep, pixel_t pfloor, pixel_t pceil, 
			 char *annot, char *units);
extern CRH_T  ** crhread(int fd);
extern int       crhsort(CRH_T *crhp, int sortkey);
extern int       crhwrite(int fd, CRH_T **crhpp);

#endif
