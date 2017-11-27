#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
  extern void poly_mult(void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"poly_mult", (DL_FUNC) &poly_mult, 5},
  {NULL, NULL, 0}
};

void R_init_PolynomF(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
