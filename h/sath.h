#ifndef	SAT_H
#define	SAT_H

/* ------------------------------------------------------------------------ */
 
typedef struct {
	char           *platform;	/* e.g. "Landsat", "ER2", ...	*/
	char           *sensor;		/* e.g. "TM", "AVIRIS", ... 	*/
	char           *location;	/* e.g. Landsat path,row,quad	*/
	char           *gmdate;		/* YYYYMMDD		 	*/
	char           *gmtime;		/* hhmmss.sss...	 	*/
} SATH_T;

#define	SATH_HNAME	"sat"		/* header name within IPW	 */
#define	SATH_VERSION	"$Revision: 1.2 $"	/* RCS revsion #	 */

/* field keywords */
#define	SATH_PLATFORM		"platform"
#define	SATH_SENSOR		"sensor"
#define	SATH_LOCATION		"location"
#define	SATH_GMDATE		"gmdate"
#define	SATH_GMTIME		"gmtime"

/* field access macros */
#define	sath_platform(p)	( (p)->platform )
#define	sath_sensor(p)	( (p)->sensor )
#define	sath_location(p) ( (p)->location )
#define	sath_gmdate(p)	( (p)->gmdate )
#define	sath_gmtime(p)	( (p)->gmtime )

/* ------------------------------------------------------------------------ */
 
/* function declarations */

extern SATH_T  ** sathdup(SATH_T **oldhpp, int nbands);
extern SATH_T   * sathmake(char *platform, char *sensor, char *location,
                           char *gmdate, char *gmtime);
extern SATH_T  ** sathread(int fd);
extern int        sathwrite(int fd, SATH_T **sathpp);

#endif
