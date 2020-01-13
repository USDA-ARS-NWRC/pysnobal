#ifndef	WIN_H
#define	WIN_H

/*
 * window header
 */

/* ------------------------------------------------------------------------ */
 
typedef struct {
	double          bline;		/* begin line #			 */
	double          bsamp;		/* begin sample #		 */
	double          dline;		/* line increment		 */
	double          dsamp;		/* sample increment		 */
} WINH_T;

#define	WINH_HNAME	"win"		/* header name within IPW	 */
#define	WINH_VERSION	"$Revision: 1.7 $"	/* RCS revsion #	 */

/* field keywords */
#define	WINH_BLINE	"bline"
#define	WINH_BSAMP	"bsamp"
#define	WINH_DLINE	"dline"
#define	WINH_DSAMP	"dsamp"

/* field access macros */
#define	winh_bline(p)	( (p)->bline )
#define	winh_bsamp(p)	( (p)->bsamp )
#define	winh_dline(p)	( (p)->dline )
#define	winh_dsamp(p)	( (p)->dsamp )

/* coordinate conversion macros */
#define	WIN_LINE(p, L)	( ((L) * (p)->dline) + (p)->bline )
#define	WIN_SAMP(p, s)	( ((s) * (p)->dsamp) + (p)->bsamp )

/* ------------------------------------------------------------------------ */
 
/* function declarations */

extern WINH_T ** winhdup(WINH_T **oldhpp, int nbands);
extern int       winhfree(WINH_T **winhpp, int nbands);
extern WINH_T  * winhmake(double bline, double bsamp, double dline,
                          double dsamp);
extern WINH_T ** winhread(int fd);
extern int       winhwrite(int fd, WINH_T **winhpp);

#endif
