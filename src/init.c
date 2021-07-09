/* > tools::package_native_routine_registration_skeleton(".") */
#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Fortran calls */
extern void F77_NAME(cubicbsplines_general)(void *, void *, void *, void *, void *);
extern void F77_NAME(d1_cubicbsplines_general)(void *, void *, void *, void *, void *);
extern void F77_NAME(d2_cubicbsplines_general)(void *, void *, void *, void *, void *);
extern void F77_NAME(integrated_cubicbsplines_general)(void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
  {"cubicbsplines_general",            (DL_FUNC) &F77_NAME(cubicbsplines_general),            5},
  {"d1_cubicbsplines_general",         (DL_FUNC) &F77_NAME(d1_cubicbsplines_general),         5},
  {"d2_cubicbsplines_general",         (DL_FUNC) &F77_NAME(d2_cubicbsplines_general),         5},
  {"integrated_cubicbsplines_general", (DL_FUNC) &F77_NAME(integrated_cubicbsplines_general), 6},
  {NULL, NULL, 0}
};

void R_init_cubicBsplines(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
