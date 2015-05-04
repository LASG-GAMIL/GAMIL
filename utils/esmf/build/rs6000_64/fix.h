/*$Id: fix.h,v 1.3.2.2 2002/04/27 15:39:06 erik Exp $*/

/*
    This fixes various things in system files that are incomplete, for 
  instance many systems don't properly prototype all system functions.
  It is not intended to DUPLICATE anything in the system include files;
  if the compiler reports a conflict between a prototye in a system file
  and this file then the prototype in this file should be removed.

    This is included by files in src/sys/src
*/

#if !defined(ESMC_FIX_H)
#define ESMC_FIX_H

#include "MC.h"

/*
  This prototype lets us resolve the datastructure 'rusage' only in
  the source files using getrusage, and not in other source files.
*/
typedef struct rusage* s_rusage;

/* ----------------------IBM RS6000 ----------------------------------------*/
/* Some of the following prototypes are present in AIX 4.2 but not in AIX 3.X */
#if defined(__cplusplus)
extern "C" {
extern char   *mktemp(char *);
extern char   *getwd(char *);
extern int    getdomainname(char *,int);
extern int    strcasecmp(const char *, const char *);
extern int    getrusage(int,s_rusage);
}
#else
extern char   *mktemp(char *);
extern int    strcasecmp(const char *, const char *);
extern int    getrusage(int,s_rusage);
#endif
#endif
