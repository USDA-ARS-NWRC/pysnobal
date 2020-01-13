#ifndef	GEO_H
#define	GEO_H

/*
 * geodetic header
 */

typedef struct {
	double          bline;		/* begin line #			 */
	double          bsamp;		/* begin sample #		 */
	double          dline;		/* line increment		 */
	double          dsamp;		/* sample increment		 */
	char           *units;		/* units of measurement		 */
	char           *csys;		/* coord system ID		 */
} GEOH_T;

#define	GEOH_HNAME	"geo"		/* header name within IPW	 */
#define	GEOH_VERSION	"$Revision: 1.7 $"	/* RCS revsion #	 */

/* field keywords */
#define	GEOH_BLINE	"bline"
#define	GEOH_BSAMP	"bsamp"
#define	GEOH_DLINE	"dline"
#define	GEOH_DSAMP	"dsamp"
#define	GEOH_UNITS	"units"
#define	GEOH_CSYS	"coord_sys_ID"

/* field access macros */
#define	geoh_bline(p)	( (p)->bline )
#define	geoh_bsamp(p)	( (p)->bsamp )
#define	geoh_dline(p)	( (p)->dline )
#define	geoh_dsamp(p)	( (p)->dsamp )
#define	geoh_units(p)	( (p)->units )
#define	geoh_csys(p)	( (p)->csys )

/* coordinate conversion macros */
#define	GEO_LINE(p, L)	( ((L) * (p)->dline) + (p)->bline )
#define	GEO_SAMP(p, s)	( ((s) * (p)->dsamp) + (p)->bsamp )

/* ------------------------------------------------------------------------ */

/* function declarations */

extern GEOH_T ** geohdup(GEOH_T **oldhpp, int nbands);
extern int       geohfree(GEOH_T **geohpp, int nbands);
extern GEOH_T  * geohmake(double bline, double bsamp, double dline,
                          double dsamp, char *units, char *csys);
extern GEOH_T ** geohread(int fd);
extern double    geohspace(GEOH_T *geoh);
extern int       geohwrite(int fd, GEOH_T **geohpp);

#endif
