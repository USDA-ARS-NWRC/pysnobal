#ifndef	IPW_UTIL_H
#define	IPW_UTIL_H

/*
 * IPW utility functions
 */

#define SAFE_FREE(p)		if ((p) != NULL) free(p)

/* ------------------------------------------------------------------------ */

extern addr_t	allocnd(int elsize, int ndim, ...);
extern addr_t	ecalloc(int nelem, int elsize);
extern int	imgcopy(int fdi, int fdo);
extern void	no_tty(int fd);

#endif  /* IPW_UTIL_H */
