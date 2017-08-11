#ifndef	IPW_STR_H
#define	IPW_STR_H

/*
 * Character and string macros and functions.
 */

#define	BLANK		' '
#define	EOS		'\0'
#define WHITE_SPACE     " \t\n"

/* ------------------------------------------------------------------------ */

#define	strdiff(s1, s2)		( strcmp(s1, s2) != 0 )
#define	strdiffn(s1, s2, n)	( strncmp(s1, s2, n) != 0 )
#define	streq(s1, s2)		( strcmp(s1, s2) == 0 )
#define	streqn(s1, s2, n)	( strncmp(s1, s2, n) == 0 )

#define	STRDIFF(s1, s2)		( (s1)[0] != (s2)[0] && strdiff(s1, s2) )
#define	STRDIFFN(s1, s2, n)	( (s1)[0] != (s2)[0] && strdiffn(s1, s2, n) )
#define	STREQ(s1, s2)		( (s1)[0] == (s2)[0] && streq(s1, s2) )
#define	STREQN(s1, s2, n)	( (s1)[0] == (s2)[0] && streqn(s1, s2, n) )

#if 0
#define	SKIP_DIGITS(p)		( (p) += strspn(p, DEC_DIGITS) )
#define	SKIP_SPACE(p)		( (p) += strspn(p, WHITE_SPACE) )
#endif
/* ------------------------------------------------------------------------ */

extern char   *	rmlead(char * str);
extern void     rmtrail(char * str);

#endif  /* IPW_STR_H */
