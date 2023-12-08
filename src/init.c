#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _psbcSpeedUp_psbcSpeedUp_internal(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"_psbcSpeedUp_psbcSpeedUp_internal", (DL_FUNC)&_psbcSpeedUp_psbcSpeedUp_internal, 13},
    {NULL, NULL, 0}};

void R_init_psbcSpeedUp(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
