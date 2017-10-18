#ifndef	BIH_H
#define	BIH_H

/*
 * basic image header
 */

/* ------------------------------------------------------------------------ */

typedef struct {
	char           *byteorder;	/* byte order			 */
	int             nlines;		/* # lines			 */
	int             nsamps;		/* # samples / line		 */
	int             nbands;		/* # bands			 */
} BIHI_T;

#define	BIHI_HNAME	"basic_image_i"	/* per-image component name	 */

typedef struct {
	BIHI_T         *img;		/* -> per-image component	 */
	int             nbytes;		/* # bytes / pixel		 */
	int             nbits;		/* # bits / pixel		 */
	STRVEC_T       *history;	/* -> processing history strings */
	STRVEC_T       *annot;		/* -> annotation strings	 */
} BIH_T;

#define	BIH_HNAME	"basic_image"	/* header name within IPW	 */
#define	BIH_VERSION	"$Revision: 1.11 $"	/* RCS revsion #	 */

/* field keywords */
#define	BIH_BYTEORDER	"byteorder"
#define	BIH_NLINES	"nlines"
#define	BIH_NSAMPS	"nsamps"
#define	BIH_NBANDS	"nbands"
#define	BIH_NBYTES	"bytes"
#define	BIH_NBITS	"bits"
#define	BIH_HISTORY	"history"
#define	BIH_ANNOT	"annot"

/* field access macros */
#define bih_byteorder(p)	( (p)->img->byteorder )
#define	bih_nlines(p)		( (p)->img->nlines )
#define	bih_nsamps(p)		( (p)->img->nsamps )
#define	bih_nbands(p)		( (p)->img->nbands )
#define	bih_nbytes(p)		( (p)->nbytes )
#define	bih_nbits(p)		( (p)->nbits )
#define	bih_history(p)		( (p)->history )
#define	bih_annot(p)		( (p)->annot )

/* ------------------------------------------------------------------------ */

extern char * o_byteorder;

/* ------------------------------------------------------------------------ */

extern BIH_T  ** bih(int fd);
extern BIH_T  ** bihdup(BIH_T **oldhpp);
extern BIH_T   * bihmake(int nbytes, int nbits, STRVEC_T *history,
			STRVEC_T *annot, BIH_T *oldp, int nlines,
			int nsamps, int nbands);
extern BIH_T  ** bihread(int fd);
extern int       bihwrite(int fd, BIH_T **bihpp);
extern char    * hbyteorder(int fd);
extern int       hnbands(int fd);
extern int       hnbits(int fd, int band);
extern int       hnbytes(int fd, int band);
extern int       hnlines(int fd);
extern int       hnsamps(int fd);
extern char    * hostorder(void);
extern long      imgsize(int fd);
extern void	 no_history(int fd);
extern int       sampsize(int fd);

#endif
