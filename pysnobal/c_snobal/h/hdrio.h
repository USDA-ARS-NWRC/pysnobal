#ifndef	HDRIO_H
#define	HDRIO_H

/*
 * header file for IPW image header I/O subsystem
 */

#define	BOIMAGE		"image"		/* name for image data preamble	 */
#define	HGOT_DATA	1		/* got a header data record	 */
#define	HGOT_PRMB	0		/* got a header preamble record	 */
#define	HREC_MAX	LINE_MAX	/* max hdr rec size (incl '\n')	 */
#define	HREC_SIZ	(HREC_MAX + 1)	/* max hdr rec size plus EOS	 */
#define	NO_BAND		(-1)

/* ------------------------------------------------------------------------ */
 
/*
 * '\f' at end of version string will stop "more" before image data
 */
#define	boimage(fd)	hwprmb(fd, BOIMAGE, NO_BAND, "$Revision: 1.5 $\f")

/* ------------------------------------------------------------------------ */
 
extern int       hcopy(int fdi, int fdo);
extern addr_t	 hdralloc(int n, int size, int fd, const char *name);
extern int       hgetrec(int fd, char *comment, char *key, char *value);
extern int       hpass(int fdi, int fdo);
extern int       hputrec(int fd, const char *comment, const char *key,
                         const char *value);
extern int       hrband(int fd);
extern char    * hrname(int fd);
extern int       hrskip(int fdi);
extern char    * hrvers(int fd);
extern int       hwprmb(int fd, const char *name, int band, 
                        const char *version);

#endif
