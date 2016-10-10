#ifndef	IPW_ERROR_H
#define	IPW_ERROR_H

/*
 * IPW error handling
 */

//#include "_types.h"

#define	assert(expr)		if (! (expr)) \
					bug("assertion \"" #expr "\" failed")
#define	bug(s)			_bug(s, __FILE__, __LINE__)

/* ------------------------------------------------------------------------ */

extern void	_bug(const char *msg, const char *file, int line);

/* ------------------------------------------------------------------------ */

extern void	error(const char *format, ...);
extern void	usrerr(const char *format, ...);
extern void	warn(const char *format, ...);
extern void	syserr(void);
extern void	uferr(int fd);

#endif  /* IPW_ERROR_H */
