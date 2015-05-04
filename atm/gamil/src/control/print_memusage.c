#include <misc.h>

#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>

#ifdef LINUX
#include <unistd.h>
#endif

#include <cfort.h>

#if ( defined FORTRANCAPS )
#define print_memusage PRINT_MEMUSAGE
#elif ( defined FORTRANUNDERSCORE )
#define print_memusage print_memusage_
#elif ( defined FORTRANDOUBLEUNDERSCORE )
#define print_memusage print_memusage__
#endif

void print_memusage ()
{
#ifndef CRAY
  struct rusage usage;

  if (getrusage (RUSAGE_SELF, &usage) < 0) {
    fprintf (stderr, "print_memusage: bad return from getrusage\n");

  } else {

    fprintf (stderr, "max rss=%ld shared mem=%ld unshared data=%ld unshared stack=%ld\n", 
	     usage.ru_maxrss, usage.ru_ixrss, usage.ru_idrss, usage.ru_isrss);

  }
#endif
}

                 
