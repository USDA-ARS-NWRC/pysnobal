#ifndef	_MAIN_H
#define	_MAIN_H

/*
 *  Macros and routines for IPW main functions and options
 */

/* ------------------------------------------------------------------------- */

/*
 *  IPW options and their arguments
 */

typedef int     aint_t;
typedef long    along_t;
typedef double  areal_t;
typedef char   *astr_t;

typedef union {
	aint_t         *aint_p;
	along_t        *along_p;
	areal_t        *areal_p;
	astr_t         *astr_p;
} ARGS_T;

typedef struct {

#ifdef STRING_OPTIONS
	const char     *option;
#else
	char            option;
#endif

	const char     *descrip;
	int             type;
	const char     *arg_descrip;
	bool_t          required;
	int             min_nargs;
	int             max_nargs;
	bool_t		commasInArg;	/* allow "," in a single STR_OPTARG? */
	int             nargs;
	ARGS_T          args;
} OPTION_T;

/* values of "option" */
#ifdef STRING_OPTIONS
#define OPERAND		NULL
#else
#define	OPERAND		0		/* DON'T CHANGE!		 */
#endif

/* values of "type" */
#define	NO_OPTARGS	0		/* DON'T CHANGE!		 */
#define	INT_OPTARGS	1
#define	REAL_OPTARGS	2
#define	STR_OPTARGS	3
#define	LONG_OPTARGS	4

#define	INT_OPERANDS	INT_OPTARGS
#define	REAL_OPERANDS	REAL_OPTARGS
#define	STR_OPERANDS	STR_OPTARGS
#define	LONG_OPERANDS	LONG_OPTARGS

/* values of "required" */
#define	OPTIONAL	FALSE		/* DON'T CHANGE			 */
#define	REQUIRED	TRUE		/* DON'T CHANGE			 */

/*
 * special constant for an option that takes a single string argument
 * which may contain commas
 */
#define ONE_ARG_WITH_COMMAS	1, 1, TRUE


#define	n_args(opt)		( (opt).nargs )
#define	got_opt(opt)		( n_args(opt) > 0 )

#define	int_argp(opt)		( (opt).args.aint_p )
#define	long_argp(opt)		( (opt).args.along_p )
#define	real_argp(opt)		( (opt).args.areal_p )
#define	str_argp(opt)		( (opt).args.astr_p )

#define	int_arg(opt, i)		( ((opt).args.aint_p)[i] )
#define	long_arg(opt, i)	( ((opt).args.along_p)[i] )
#define	real_arg(opt, i)	( ((opt).args.areal_p)[i] )
#define	str_arg(opt, i)		( ((opt).args.astr_p)[i] )

#define	n_opnds(opt)		n_args(opt)
#define	got_opnds(opt)		got_opt(opt)

#define	str_opnd(opt, i)	str_arg(opt, i)

/* ------------------------------------------------------------------------- */
 
/*
 *  Functions for IPW main functions and options
 */

extern void     ipwenter(int argc, char **argv, OPTION_T **optv,
                         const char *descrip);
extern void     ipwexit(int status);
extern void	opt_check(int n_min, int n_max, int n_opts,
			  OPTION_T *opt, ...);
extern void     usage(void);

#endif
