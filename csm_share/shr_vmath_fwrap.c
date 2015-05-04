/* =============================================================================
** CVS: $Id: shr_vmath_fwrap.c,v 1.1 2001/10/15 22:28:52 kauff Exp $
** CVS: $Source: /fs/cgd/csm/models/CVS.REPOS/shared/csm_share/shr_vmath_fwrap.c,v $
** CVS: $Name: cam2_0_1_brnchT_release3 $
** =============================================================================
** Fortran wrappers for shr_vmath calls for systems that only
** provide a C interface to the system vector math routines.
** =============================================================================
*/

#if (defined IRIX64)

void shr_vmath_fwrap_vsqrt_(double *X, double *Y, int *n)
{
   vsqrt(X, Y, *n, 1, 1);
}

void shr_vmath_fwrap_vexp_(double *X, double *Y, int *n)
{
   vexp(X, Y, *n, 1, 1);
}

void shr_vmath_fwrap_vlog_(double *X, double *Y, int *n)
{
   vlog(X, Y, *n, 1, 1);
}

void shr_vmath_fwrap_vsin_(double *X, double *Y, int *n)
{
   vsin(X, Y, *n, 1, 1);
}

void shr_vmath_fwrap_vcos_(double *X, double *Y, int *n)
{
   vcos(X, Y, *n, 1, 1);
}

#endif
