#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP DGrowthR_exact_gaussian_fit_state(SEXP X_, SEXP y_, SEXP lengthscale_,
                                       SEXP sigma2_, SEXP nugget_);
SEXP DGrowthR_exact_gaussian_fit_state_dist(SEXP dist_sq_, SEXP y_,
                                            SEXP lengthscale_, SEXP sigma2_,
                                            SEXP nugget_);
SEXP DGrowthR_exact_gaussian_lbfgsb_dist(SEXP start_, SEXP lower_, SEXP upper_,
                                         SEXP dist_sq_, SEXP y_,
                                         SEXP nugget_weights_, SEXP maxit_);

static const R_CallMethodDef CallEntries[] = {
  {"DGrowthR_exact_gaussian_fit_state",
   (DL_FUNC) &DGrowthR_exact_gaussian_fit_state, 5},
  {"DGrowthR_exact_gaussian_fit_state_dist",
   (DL_FUNC) &DGrowthR_exact_gaussian_fit_state_dist, 5},
  {"DGrowthR_exact_gaussian_lbfgsb_dist",
   (DL_FUNC) &DGrowthR_exact_gaussian_lbfgsb_dist, 7},
  {NULL, NULL, 0}
};

void R_init_DGrowthR(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
