#ifndef IPW_ERRNO_H
#define IPW_ERRNO_H

/*
 *  Header file for global variable with error code: "ipw_errno"
 */

/* ------------------------------------------------------------------------ */
 
/*
 *  Possible error codes
 */

#define IPWERR_MESSAGE	1	/* (message saved via "usrerr" routine)	*/

#define IPWERR_MEMORY	2	/* Insufficient memory			*/

/*
 *  header I/O
 */
#define IPWERR_RD_BIH	3	/* Cannot read BI image header		*/
#define IPWERR_WR_BIH	4	/* Cannot write BI header		*/
#define IPWERR_WR_BOI	5	/* Cannot write image-data header	*/
#define IPWERR_WR_GEOH	6	/* Cannot write GEO header		*/
#define IPWERR_WR_LQH	7	/* Cannot write LQ header		*/

/* ------------------------------------------------------------------------ */
 
extern int ipw_errno;

/* ------------------------------------------------------------------------ */
 
extern char * ipw_strerr(int errno);

#endif  /* IPW_ERRNO_H */
