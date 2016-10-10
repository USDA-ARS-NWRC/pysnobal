#ifndef	FPIO_H
#define	FPIO_H

extern int         fpclose(int fd);
extern fpixel_t  * fpfmax(int fd);
extern fpixel_t  * fpfmin(int fd);
extern void        fphdrs(int fdi, int nbands, int fdo);
extern fpixel_t ** fpmap(int fd);
extern int       * fpmaplen(int fd);
extern int         fpvread(int fd, fpixel_t *buf, int npixv);
extern int         fpvwrite(int fd, fpixel_t *buf, int npixv);
extern int         mnxfp(fpixel_t *x, int npixv, int nbands, fpixel_t *mmval);

#endif
