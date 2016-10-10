#ifndef GCTP_H
#define GCTP_H


#if defined(__STDC__) || defined(__cplusplus)
# define P_(s) s
#else
# define P_(s) ()
#endif

/* gtrnz0.c */
int gtrnz0_ P_((double *crdin, int *indef, double *tparin, double *crdio, int *iodef, double *tpario, FILE *ipfile, int *iflg));

#undef P_

#endif  /* GCTP_H */
