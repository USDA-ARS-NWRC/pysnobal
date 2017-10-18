#ifndef	PIXIO_H
#define	PIXIO_H

extern int      pvread(int fd, pixel_t *buf, int npixv);
extern int      pvwrite(int fd, pixel_t *buf, int npixv);
extern int      pxclose(int fd);

#endif
