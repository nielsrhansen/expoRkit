#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(dgpadm)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(r_dgexpv)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(r_dgphiv)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(r_dmexpv)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(r_zgexpv)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(r_zgphiv)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zgpadm)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"dgpadm",   (DL_FUNC) &F77_NAME(dgpadm),   11},
    {"r_dgexpv", (DL_FUNC) &F77_NAME(r_dgexpv), 15},
    {"r_dgphiv", (DL_FUNC) &F77_NAME(r_dgphiv), 16},
    {"r_dmexpv", (DL_FUNC) &F77_NAME(r_dmexpv), 15},
    {"r_zgexpv", (DL_FUNC) &F77_NAME(r_zgexpv), 15},
    {"r_zgphiv", (DL_FUNC) &F77_NAME(r_zgphiv), 16},
    {"zgpadm",   (DL_FUNC) &F77_NAME(zgpadm),   11},
    {NULL, NULL, 0}
};

void R_init_expoRkit(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
