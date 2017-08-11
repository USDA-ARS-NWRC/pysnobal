#include <errno.h>
#include <math.h>

//#include "ipw.h"
#include "envphys.h"

double
sati(
		double  tk)		/* air temperature (K)	*/
{
	double  l10;
	double  x;

	if (tk <= 0.) {
		fprintf(stderr, "tk=%f\n, less than zero", tk);
//		assert(tk > 0.);
		exit(EXIT_FAILURE);
	}

	if (tk > FREEZE) {
		x = satw(tk);
		return(x);
	}

	errno = 0;
	l10 = log(1.e1);

	x = pow(1.e1,-9.09718*((FREEZE/tk)-1.) - 3.56654*log(FREEZE/tk)/l10 +
			8.76793e-1*(1.-(tk/FREEZE)) + log(6.1071)/l10);

	if (errno) {
		perror("sati: bad return from log or pow");
		exit(EXIT_FAILURE);
//		syserr();
//		error("sati: bad return from log or pow");
	}

	return(x*1.e2);
}
