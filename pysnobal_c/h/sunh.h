#ifndef	SUNH_H
#define	SUNH_H

/*
**	sunh -- sun geometry header
*/

/* ------------------------------------------------------------------------ */
 
typedef struct {
	double          cos_sun;	/* cosine solar zenith ang	 */
	double          zen_sun;	/* solar zenith angle (radians)	 */
	double          azm_sun;	/* solar azimuth (rad from S)	 */
}               SUNH_T;

#define	SUNH_HNAME	"sun"		/* header name within IPW	 */
#define	SUNH_VERSION	"$Revision: 1.1 $"	/* RCS revsion #	 */

/* field keywords */
#define	SUNH_COS	"cos_sun"
#define	SUNH_ZEN	"zen_sun"
#define	SUNH_AZM	"azm_sun"

/* field access macros */
#define sunh_cos(p)	((p)->cos_sun)
#define sunh_zen(p)	((p)->zen_sun)
#define sunh_azm(p)	((p)->azm_sun)

/* ------------------------------------------------------------------------ */
 
/* function declarations */

extern SUNH_T  * sunhmake(double cos_zen, double azm);
extern SUNH_T ** sunhread(int fd);
extern int       sunhwrite(int fd, SUNH_T **sunhpp);

#endif
