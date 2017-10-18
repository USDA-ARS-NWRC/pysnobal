#ifndef	IPW_DATETIME_H
#define	IPW_DATETIME_H

/*
 * Date and time functions for IPW
 */

typedef struct {
	int	year;		/* year (note: 90 = 90 AD, not 1990) 	*/
	int	month;		/* month [1 - 12]			*/
	int	day;		/* day of month [1 - 31]		*/
	int	hour;		/* time of day: hour [0 - 23]		*/
	int	min;		/*  "   "   " : minutes [0 - 59]	*/
	int	sec;		/*  "   "   " : seconds [0 - 59]	*/
} datetime_t;

#define SECS_IN_MIN		60
#define SECS_IN_HR		3600
#define SECS_IN_DAY		86400
#define MINS_IN_HR		60
#define HRS_IN_DAY		24

/* ------------------------------------------------------------------------ */

#define HR_TO_MIN(h)		( (h) * 60 )
#define HR_TO_SEC(h)		( (h) * 3600 )
#define MIN_TO_SEC(m)		( (m) * 60 )
#define MIN_TO_HR(m)		( (m) / (float) MINS_IN_HR )
#define SEC_TO_MIN(s)		( (s) / (float) SECS_IN_MIN )
#define SEC_TO_HR(s)		( (s) / (float) SECS_IN_HR )
#define SEC_TO_DAY(s)		( (s) / (float) SECS_IN_DAY )
#define DAY_TO_SEC(d)		( (d) * SECS_IN_DAY )

#define HMS_TO_SEC(h,m,s)	( HR_TO_SEC(h) + MIN_TO_SEC(m) + (s) )
/*
 *  Convert a time in hour:minutes:seconds into total seconds.
 */

#define DAYS_IN_YR(y)		( leapyear(y) ? 366 : 365 )
/*
 *  Return number of days in a given year
 */

/* ------------------------------------------------------------------------ */

extern void	    add_to_dt(datetime_t * dt, int days, int seconds);
extern bool_t	    dt2fmtstr(datetime_t * dt, char * format, char * buffer);
extern char *	    dt2str(datetime_t * dt);
extern void	    dt_diff(datetime_t * dt1, datetime_t * dt2,
			    int * days, int * seconds);
extern bool_t	    dt_in_order(datetime_t * dt1, datetime_t * dt2);
extern void	    gmt2local(datetime_t * dt, int zone, bool_t isDST);
extern bool_t       leapyear(int year);
extern void	    local2gmt(datetime_t * dt, int zone, bool_t isDST);
extern datetime_t * make_dt(int year, int month, int day, int hour, int min,
			    int sec);
extern datetime_t * now_dt(void);
extern int	    numdays(int year, int month);
extern void	    sec2hms(int total_secs, int * hrs, int * mins, int * secs);
extern bool_t	    str2dt(char * str, char * format, datetime_t * dt);
extern int	    waterday(int year, int month, int day, int * wateryear);
extern void	    wday2mday(int wateryear, int waterday, int * year,
			      int * month, int * day);
extern int	    weekday(int year, int month, int day);
extern void	    yday2mday(int year, int yrday, int * month, int * day);
extern int	    yearday(int year, int month, int day);
extern char *	    zone2str(int zone, int isdst);

#endif  /* IPW_DATETIME_H */
