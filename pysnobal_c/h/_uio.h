#ifndef	IPW_UIO_H
#define	IPW_UIO_H

/*
 * IPW uio (Unix IO)
 */

#define NO_FD		(-1)		/* no file descriptor (file
					   descriptors are >= 0) */

#define ASSERT_OK_FD(fd)	assert(((fd) >= 0) && ((fd) <= _ipwmaxfd))

extern int	_ipwmaxfd;
/*
 *	Maximum allowable IPW file descriptor.
 */

/* ------------------------------------------------------------------------ */

extern char    *mktemplate(const char *prefix);
extern bool_t   ubof(int fd);
extern int      uclose(int fd);
extern long     ucopy(int fdi, int fdo, long ncopy);
extern bool_t   ueof(int fd);
extern char    *ufilename(int fd);
extern char    *ugets(int fd, char *buf, int nbytes);
extern int      uputs(int fd, const char *buf);
extern int      uread(int fd, addr_t buf, int nbytes);
extern int      uremove(const char *filename);
extern int      uropen(const char *name);
extern long     urskip(int fd, long nbytes);
extern int      ustdin(void);
extern int      ustdout(void);
extern int      uwflush(int fd);
extern int      uwopen(const char *name);
extern int      uwopen_temp(char *name);
extern int      uwrite(int fd, const addr_t buf, int nbytes);

#endif  /* IPW_UIO_H */
