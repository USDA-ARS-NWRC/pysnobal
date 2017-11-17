#ifndef	OR_H
#define	OR_H

/*
 * orientation header
 */

typedef struct {
	char           *orient;		/* orientation			 */
	char           *origin;		/* corner of origin		 */
} ORH_T;

#define	ORH_HNAME	"or"		/* header name within IPW	 */
#define	ORH_VERSION	"$Revision: 1.4 $"	/* RCS revsion #	 */

/* field keywords */
#define	ORH_OR		"orientation"
#define	ORH_ORIGIN	"origin"

/* field access macros */
#define orh_orient(p)	( (p)->orient )
#define orh_origin(p)	( (p)->origin )

/* definitions */
#define ROW		"row-major"
#define COLUMN		"column-major"
#define ORIG_1		"upper left"
#define ORIG_2		"upper right"
#define ORIG_3		"lower right"
#define ORIG_4		"lower left"
#define IPW_ORIENT	ROW
#define IPW_ORIGIN	ORIG_1
#define XT_ORIENT	COLUMN

/* ------------------------------------------------------------------------ */

/* function declarations */

extern ORH_T  ** orhdup(ORH_T **oldhpp, int nbands);
extern int       orhfree(ORH_T **orhpp, int nbands);
extern ORH_T   * orhmake(char *orient, char *origin);
extern ORH_T  ** orhread(int fd);
extern int       orhwrite(int fd, ORH_T **orhpp);

#endif
