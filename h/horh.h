#ifndef	HOR_H
#define	HOR_H

/*
 * horizon header
 */

typedef struct {
	double          azimuth;	/* direction toward horizons	 */
} HORH_T;

#define	HORH_HNAME	"hor"		/* header name within IPW	 */
#define	HORH_VERSION	"$Revision: 1.3 $"	/* RCS revsion #	 */

/* field keywords */
#define	HORH_AZM	"azimuth"

/* field access macros */
#define	horh_azm(p)	( (p)->azimuth )

/* ------------------------------------------------------------------------ */

/* function declarations */

extern HORH_T  * horhmake(double azm);
extern HORH_T ** horhread(int fd);
extern int       horhwrite(int fd, HORH_T **horhpp);

#endif
