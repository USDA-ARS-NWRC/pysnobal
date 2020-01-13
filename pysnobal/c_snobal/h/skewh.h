#ifndef	SKEWH_H
#define	SKEWH_H

/*
**	The presence of as skew header indicates that an image has been skewed
**	horizontally by the skew program.  The skew angle in the header is
**	used by skew to undo the skewing.
*/

/* ------------------------------------------------------------------------ */
 
typedef struct {
	double          angle;		/* skew angle (degrees)		 */
} SKEWH_T;

#define	SKEWH_HNAME	"skew"		/* header name within IPW	 */
#define	SKEWH_VERSION	"$Revision: 1.3 $"	/* RCS revsion #	 */

/* field keywords */
#define	SKEWH_ANGLE		"angle"

/* field access macros */
#define	skewh_angle(p)	( (p)->angle )

/* misc. */
#define	SKEW_MAX	45.0		/* maximum skew angle		 */
#define	SKEW_MIN	(-45.0)		/* minimum skew angle		 */

/* ------------------------------------------------------------------------ */
 
/* function declarations */

extern SKEWH_T  * skewhmake(double angle);
extern SKEWH_T ** skewhread(int fd);
extern int        skewhwrite(int fd, SKEWH_T **skewhpp);

#endif
