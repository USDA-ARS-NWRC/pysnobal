#ifndef	TOPO_H
#define	TOPO_H

/*
 *  Topographic calculations
 */

extern fpixel_t	basin(int seedX, int seedY, fpixel_t * slope,
		      fpixel_t * exposure, fpixel_t * elev, pixel_t * bmask,
		      int x_count, int y_count, fpixel_t flat, int clipFlag,
		      fpixel_t clipElev);

#endif
