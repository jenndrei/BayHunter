/******************************************************************
 *                                                                *
 *   Copyright by Joachim Saul <saul@gfz-potsdam.de>              *
 *                                                                *
 ******************************************************************/

/* Uncomment this if you prefer to specify
 * the slowness (ray parameter) in sec/km
 * Otherwise the slowness is specified in
 * sec/degree.
 */
/* # define SECPERKM */

/* 
 * **********************************************
 * * Default values for command line parameters *
 * **********************************************
 */

/* Default value for horizontal slowness.
 * For an epicentral distance of 80 degrees:
 * 5.4 sec/degree  or 0.0486 sec/km
 */
# define PDEF_P      5.4

# define PDEF_FSAMP  20.0
# define PDEF_NSAMP  1024
# define PDEF_WATER  0.00001
# define PDEF_GAUSS  8.0
# define PDEF_TSHIFT 5.0
# define PDEF_SMOOTH 0.02
# define PDEF_ISPEED 0.5
# define PDEF_MFILE  "mod"
# define PDEF_RFFILE "rf"
# define PDEF_PREFIX PDEF_MFILE

# define DEF_PR 0.25
# define DEF_QP 450.0
# define DEF_QS 200.0
