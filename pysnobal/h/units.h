#ifndef UNITS_H
#define UNITS_H

/*
 *  Units of measurement
 */

/* ------------------------------------------------------------------------ */

typedef int UNITS_T;


#define	U_NONE				0

#define	U_PERCENT			1

#define	U_CELSIUS			2
#define	U_KELVIN			3	
#define	U_FARENHEIT			4

#define U_METERS			5
#define U_KILOMETERS			6

#define	U_KILOGRAMS			7

#define U_PASCALS			8

#define U_METERS_PER_SECOND		9

#define U_KILOGRAMS_PER_SQUARE_METER	10
#define U_KILOGRAMS_PER_CUBIC_METER	11
#define U_WATTS_PER_SQUARE_METER	12
#define U_JOULES_PER_SQUARE_METER	13

#define _U_MAXIMUM_UNITS_ID		13

/* ------------------------------------------------------------------------ */

#define VALID_UNITS_ID(id)	((0 <= (id)) && ((id) <= _U_MAXIMUM_UNITS_ID))

/* ------------------------------------------------------------------------ */

extern int    units_match(char * string, UNITS_T units_id);
extern char * units_as_str(UNITS_T units_id);


#endif  /* UNITS_H */
