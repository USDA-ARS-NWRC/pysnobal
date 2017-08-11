#ifndef	LQH_H
#define	LQH_H

/*
 * Linear-quantization (LQ) headers for IPW images.
 */

/* ----------------------------------------------------------------------- */

typedef struct {
	int            *bkpt;		/* breakpoint array		 */
	fpixel_t       *map;		/* LQ = map[pixel]		 */
	char           *units;		/* units of measure		 */
	char           *interp;		/* interpolation function name	 */
 /*
  * stuff below here is derived, not read from header
  */
	int             nbkpts;		/* # breakpoints		 */
	int             maplen;		/* # elements / map array	 */
	bool_t          mapped;		/* ? map has been interpolated	 */
	double         *lininv;		/* linear inverse map coefs.	 */
} LQH_T;

#define	LQH_HNAME	"lq"		/* header name within IPW	 */
#define	LQH_VERSION	"$Revision: 1.6 $"	/* RCS revsion #	 */

/* field keywords */
#define	LQH_MAP		"map"
#define	LQH_UNITS	"units"
#define	LQH_INTERP	"interp"

/* field access macros */
#define	lqh_bkpt(p)	( (p)->bkpt )
#define	lqh_map(p)	( (p)->map )
#define	lqh_units(p)	( (p)->units )
#define	lqh_interp(p)	( (p)->interp )
#define	lqh_nbkpts(p)	( (p)->nbkpts )
#define	lqh_maplen(p)	( (p)->maplen )
#define	lqh_mapped(p)	( (p)->mapped )
#define	lqh_lininv(p)	( (p)->lininv )

/* "iname" values */
#define	LQH_LIN_INTERP	"linear"

/* ----------------------------------------------------------------------- */

/* function declarations */

extern LQH_T  ** lqh(int fd);
extern LQH_T  ** lqhdup(LQH_T **oldhpp, int nbands);
extern LQH_T   * lqhmake(int nbits, int nbkpts, pixel_t *ival, 
                        fpixel_t *fval, char *units, char *interp);
extern LQH_T  ** lqhread(int fd);
extern int       lqhwrite(int fd, LQH_T **lqhpp);

#endif
