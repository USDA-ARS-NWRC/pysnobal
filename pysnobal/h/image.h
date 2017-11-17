#ifndef	IMAGE_H
#define	IMAGE_H

#include "bih.h"
#include "lqh.h"
#include "geoh.h"

#include "ipw_errno.h"

/*
 * header file for IPW image data structures
 */

/* ------------------------------------------------------------------------ */

typedef struct {
	char	       *name;		/* filename of image 		 */
	int             fd;		/* image file descriptor	 */
	unsigned int	flags;		/* flags (access, IO status)	 */

#define	IPW_READ	BIT(0)		/* read access to image		 */
#define IPW_WRITE	BIT(1)		/* write access to image	 */

#define	PIXEL_ACCESS	BIT(2)		/* access pixels as integers	 */
#define	FPIXEL_ACCESS	BIT(3)		/* access pixels as floating-pt	 */

#define	IPW_DATA_IO	BIT(4)		/* read/written any image data?  */

	int		nlines;		/* # of lines in image		 */
	int		nsamples;	/* # of samples per line	 */
	int		nbands;		/* # of bands per sample	 */
	int	       *nbits;		/* array of # of bits per band   */

	BIH_T	      **bih;		/* array of BI headers for bands */
	LQH_T	      **lqh;		/*   "   "  LQ headers  "    "   */
	GEOH_T	      **geoh;		/*   "   "  GEO headers "    "   */

	FPIXEL	       *lqh_mins;	/* array of minimums for LQ hdrs */
	FPIXEL	       *lqh_maxs;	/* array of minimums for LQ hdrs */
} IMAGE_T;

#define	IPW_STDIN	"-"
#define IPW_STDOUT	"-"

/* ------------------------------------------------------------------------ */

/*
 *  Basic image routines
 */

extern int        access(IMAGE_T *image);
extern STRVEC_T * annotation(IMAGE_T *image, int band);
extern bool_t     close_image(IMAGE_T *image);
extern GEOH_T *   get_geoh(IMAGE_T *image, int band);
extern bool_t     has_lqh(IMAGE_T *image, int band);
extern STRVEC_T * history(IMAGE_T *image, int band);
extern char *     image_name(IMAGE_T *image);
extern FPIXEL     lqh_max(IMAGE_T *image, int band);
extern FPIXEL     lqh_min(IMAGE_T *image, int band);
extern int        nbands(IMAGE_T *image);
extern FPIXEL *   new_fpbuf(int nsamples, int nbands);
extern IMAGE_T *  new_image(char *name, int nlines, int nsamples, int nbands,
			    int *nbits);
extern PIXEL *    new_pbuf(int nsamples, int nbands);
extern int        nlines(IMAGE_T *image);
extern void       not_a_tty(IMAGE_T *image);
extern int        nsamples(IMAGE_T *image);
extern IMAGE_T *  open_image(char *name, int access);
extern bool_t     read_fpbuf(IMAGE_T *image, FPIXEL *buffer, int nsamples);
extern bool_t     read_pbuf(IMAGE_T *image, PIXEL *buffer, int nsamples);
extern void       set_annotation(IMAGE_T *image, int band,
				 STRVEC_T *annotation);
extern GEOH_T *   set_geoh(IMAGE_T *image, int band, GEOH_T *geoh);
extern void       set_history(IMAGE_T *image, int band, STRVEC_T *history);
extern void       set_lqh(IMAGE_T *image, FPIXEL mins[], FPIXEL maxs[],
			  char **units);
extern bool_t     write_fpbuf(IMAGE_T *image, FPIXEL *buffer, int nsamples);
extern bool_t     write_pbuf(IMAGE_T *image, PIXEL *buffer, int nsamples);

/* ------------------------------------------------------------------------ */

/*
 *  Mask image routines
 */

extern IMAGE_T *  open_mask(char *file, int inLines, int inSamples);

#endif
