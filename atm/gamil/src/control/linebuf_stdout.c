/* 
*** FORTRAN-callable routine to line buffer stdout.  Should only be used in
*** debug mode due to system overhead, and the fact that in theory stdout
*** should not have already been written to when this routine is called.
***
*** Original version: Jim Rosinski
***
 $Id: linebuf_stdout.c,v 1.4 1999/11/12 22:04:21 rosinski Exp $
 $Author: rosinski $
*/

#include <misc.h>
#include <cfort.h>

#include <stdio.h>

#if ( defined FORTRANCAPS )

#define linebuf_stdout LINEBUF_STDOUT

#elif ( defined FORTRANUNDERSCORE )

#define linebuf_stdout linebuf_stdout_

#elif ( defined FORTRANDOUBLEUNDERSCORE )

#define linebuf_stdout linebuf_stdout__

#endif

void linebuf_stdout ()
{
  setvbuf (stdout, NULL, _IOLBF, 0);
  return;
}
