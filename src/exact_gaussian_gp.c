#include <R.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>
#include <float.h>
#include <math.h>
#include <string.h>

static void fill_gaussian_covariance(double *K, const double *X,
                                     const double *lengthscale,
                                     const double *nugget,
                                     int n, int p, double sigma2,
                                     int nugget_len, double jitter) {
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      double dist = 0.0;
      for (int k = 0; k < p; ++k) {
        double diff = X[i + n * k] - X[j + n * k];
        dist += (diff * diff) / lengthscale[k];
      }
      K[i + n * j] = sigma2 * exp(-dist);
    }
  }

  for (int i = 0; i < n; ++i) {
    double nugget_i = nugget[(nugget_len == 1) ? 0 : i];
    K[i + n * i] += sigma2 * nugget_i + jitter;
  }
}

static void fill_gaussian_covariance_dist(double *K, const double *dist_sq,
                                          const double *lengthscale,
                                          const double *nugget,
                                          int n, int p, double sigma2,
                                          int nugget_len, double jitter) {
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      double dist = 0.0;
      if (p == 1) {
        dist = dist_sq[i + n * j] / lengthscale[0];
      } else {
        for (int k = 0; k < p; ++k) {
          dist += dist_sq[i + n * j + n * n * k] / lengthscale[k];
        }
      }
      K[i + n * j] = sigma2 * exp(-dist);
    }
  }

  for (int i = 0; i < n; ++i) {
    double nugget_i = nugget[(nugget_len == 1) ? 0 : i];
    K[i + n * i] += sigma2 * nugget_i + jitter;
  }
}

typedef struct {
  const double *dist_sq;
  const double *y;
  const double *nugget_weights;
  int n;
  int p;
  int nugget_len;
} ExactProfileData;

static int fit_profile_correlation(const ExactProfileData *data,
                                   const double *par,
                                   double *chol,
                                   double *alpha,
                                   double *Rinv,
                                   double *lengthscale,
                                   double *nugget_scale,
                                   double *q,
                                   double *logdet_half) {
  int n = data->n;
  int p = data->p;

  for (int k = 0; k < p; ++k) {
    lengthscale[k] = exp(par[k]);
    if (!R_FINITE(lengthscale[k]) || lengthscale[k] <= 0.0) {
      return 1;
    }
  }

  *nugget_scale = exp(par[p]);
  if (!R_FINITE(*nugget_scale) || *nugget_scale <= 0.0) {
    return 1;
  }

  fill_gaussian_covariance_dist(chol, data->dist_sq, lengthscale,
                                data->nugget_weights, n, p, 1.0,
                                data->nugget_len, 0.0);
  for (int i = 0; i < n; ++i) {
    chol[i + n * i] += *nugget_scale * data->nugget_weights[(data->nugget_len == 1) ? 0 : i] -
      data->nugget_weights[(data->nugget_len == 1) ? 0 : i];
  }

  double jitter = 0.0;
  int info = 1;
  for (int attempt = 0; attempt < 8; ++attempt) {
    if (attempt > 0) {
      fill_gaussian_covariance_dist(chol, data->dist_sq, lengthscale,
                                    data->nugget_weights, n, p, 1.0,
                                    data->nugget_len, jitter);
      for (int i = 0; i < n; ++i) {
        chol[i + n * i] += *nugget_scale * data->nugget_weights[(data->nugget_len == 1) ? 0 : i] -
          data->nugget_weights[(data->nugget_len == 1) ? 0 : i];
      }
    }
    F77_CALL(dpotrf)("U", &n, chol, &n, &info FCONE);
    if (info == 0) {
      break;
    }
    jitter = (attempt == 0) ? sqrt(DBL_EPSILON) : jitter * 10.0;
  }

  if (info != 0) {
    return 1;
  }

  for (int i = 0; i < n; ++i) {
    alpha[i] = data->y[i];
  }

  int nrhs = 1;
  F77_CALL(dpotrs)("U", &n, &nrhs, chol, &n, alpha, &n, &info FCONE);
  if (info != 0) {
    return 1;
  }

  *q = 0.0;
  for (int i = 0; i < n; ++i) {
    *q += data->y[i] * alpha[i];
  }
  if (!R_FINITE(*q) || *q <= 0.0) {
    return 1;
  }

  *logdet_half = 0.0;
  for (int i = 0; i < n; ++i) {
    *logdet_half += log(chol[i + n * i]);
  }

  if (Rinv != NULL) {
    memset(Rinv, 0, ((size_t) n) * ((size_t) n) * sizeof(double));
    for (int i = 0; i < n; ++i) {
      Rinv[i + n * i] = 1.0;
    }
    F77_CALL(dpotrs)("U", &n, &n, chol, &n, Rinv, &n, &info FCONE);
    if (info != 0) {
      return 1;
    }
  }

  return 0;
}

static double exact_profile_lbfgsb_fn(int npar, double *par, void *ex) {
  ExactProfileData *data = (ExactProfileData *) ex;
  int n = data->n;
  int p = data->p;
  if (npar != p + 1) {
    return DBL_MAX / 4.0;
  }

  double *chol = Calloc(((size_t) n) * ((size_t) n), double);
  double *alpha = Calloc(n, double);
  double *lengthscale = Calloc(p, double);
  double nugget_scale = 0.0;
  double q = 0.0;
  double logdet_half = 0.0;

  int failed = fit_profile_correlation(data, par, chol, alpha, NULL,
                                       lengthscale, &nugget_scale, &q,
                                       &logdet_half);
  Free(chol);
  Free(alpha);
  Free(lengthscale);

  if (failed) {
    return DBL_MAX / 4.0;
  }

  return 0.5 * ((double) n) * (log(q / ((double) n)) + 1.0 + log(2.0 * M_PI)) +
    logdet_half;
}

static void exact_profile_lbfgsb_gr(int npar, double *par, double *gr, void *ex) {
  ExactProfileData *data = (ExactProfileData *) ex;
  int n = data->n;
  int p = data->p;

  for (int i = 0; i < npar; ++i) {
    gr[i] = 0.0;
  }
  if (npar != p + 1) {
    return;
  }

  double *chol = Calloc(((size_t) n) * ((size_t) n), double);
  double *alpha = Calloc(n, double);
  double *Rinv = Calloc(((size_t) n) * ((size_t) n), double);
  double *lengthscale = Calloc(p, double);
  double nugget_scale = 0.0;
  double q = 0.0;
  double logdet_half = 0.0;

  int failed = fit_profile_correlation(data, par, chol, alpha, Rinv,
                                       lengthscale, &nugget_scale, &q,
                                       &logdet_half);
  if (failed) {
    Free(chol);
    Free(alpha);
    Free(Rinv);
    Free(lengthscale);
    return;
  }

  for (int k = 0; k < p; ++k) {
    double trace_term = 0.0;
    double quad_term = 0.0;
    double inv_lengthscale = 1.0 / lengthscale[k];

    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < n; ++i) {
        double dist = 0.0;
        if (p == 1) {
          dist = data->dist_sq[i + n * j] * inv_lengthscale;
        } else {
          for (int kk = 0; kk < p; ++kk) {
            dist += data->dist_sq[i + n * j + n * n * kk] / lengthscale[kk];
          }
        }
        double dk = (p == 1) ? data->dist_sq[i + n * j] :
          data->dist_sq[i + n * j + n * n * k];
        double dR = exp(-dist) * dk * inv_lengthscale * inv_lengthscale;

        trace_term += Rinv[j + n * i] * dR;
        quad_term += alpha[i] * dR * alpha[j];
      }
    }

    gr[k] = lengthscale[k] *
      (0.5 * trace_term - 0.5 * ((double) n) * quad_term / q);
  }

  double trace_term = 0.0;
  double quad_term = 0.0;
  for (int i = 0; i < n; ++i) {
    double weight = data->nugget_weights[(data->nugget_len == 1) ? 0 : i];
    trace_term += Rinv[i + n * i] * weight;
    quad_term += alpha[i] * alpha[i] * weight;
  }
  gr[p] = nugget_scale *
    (0.5 * trace_term - 0.5 * ((double) n) * quad_term / q);

  Free(chol);
  Free(alpha);
  Free(Rinv);
  Free(lengthscale);
}

static SEXP fit_state_from_covariance_builder(SEXP y_, SEXP lengthscale_,
                                              SEXP sigma2_, SEXP nugget_,
                                              int n, int p,
                                              void (*fill_covariance)(double *, const double *,
                                                                      const double *, const double *,
                                                                      int, int, double, int, double),
                                              const double *covariate_data) {
  if (!isReal(y_)) {
    error("y_centered must be a numeric vector.");
  }
  if (LENGTH(y_) != n) {
    error("Length of y_centered must match covariance dimension.");
  }
  if (LENGTH(lengthscale_) != p) {
    error("Length of lengthscale must match number of covariates.");
  }

  const double *y = REAL(y_);
  const double *lengthscale = REAL(lengthscale_);
  double sigma2 = REAL(sigma2_)[0];
  const double *nugget = REAL(nugget_);
  int nugget_len = LENGTH(nugget_);

  if (!R_FINITE(sigma2) || sigma2 <= 0.0) {
    error("sigma2 must be positive and finite.");
  }
  if (nugget_len != 1 && nugget_len != n) {
    error("nugget must have length 1 or covariance dimension.");
  }
  for (int i = 0; i < nugget_len; ++i) {
    if (!R_FINITE(nugget[i]) || nugget[i] <= 0.0) {
      error("All nugget values must be positive and finite.");
    }
  }
  for (int k = 0; k < p; ++k) {
    if (!R_FINITE(lengthscale[k]) || lengthscale[k] <= 0.0) {
      error("All lengthscale values must be positive and finite.");
    }
  }

  SEXP chol_ = PROTECT(allocMatrix(REALSXP, n, n));
  SEXP alpha_ = PROTECT(allocVector(REALSXP, n));

  double *K = REAL(chol_);
  double *alpha = REAL(alpha_);

  double jitter = 0.0;
  int info = 1;

  for (int attempt = 0; attempt < 8; ++attempt) {
    fill_covariance(K, covariate_data, lengthscale, nugget, n, p, sigma2,
                    nugget_len, jitter);
    F77_CALL(dpotrf)("U", &n, K, &n, &info FCONE);
    if (info == 0) {
      break;
    }
    jitter = (attempt == 0) ? sqrt(DBL_EPSILON) : jitter * 10.0;
  }

  if (info != 0) {
    UNPROTECT(2);
    error("Exact Gaussian covariance matrix is not positive definite, even after adding jitter.");
  }

  for (int i = 0; i < n; ++i) {
    alpha[i] = y[i];
  }

  int nrhs = 1;
  F77_CALL(dpotrs)("U", &n, &nrhs, K, &n, alpha, &n, &info FCONE);
  if (info != 0) {
    UNPROTECT(2);
    error("Failed to solve exact Gaussian GP linear system.");
  }

  double quad = 0.0;
  for (int i = 0; i < n; ++i) {
    quad += y[i] * alpha[i];
  }

  double logdet_half = 0.0;
  for (int i = 0; i < n; ++i) {
    logdet_half += log(K[i + n * i]);
  }

  double neg_loglik = 0.5 * quad + logdet_half +
    0.5 * ((double) n) * log(2.0 * M_PI);

  SEXP jitter_ = PROTECT(ScalarReal(jitter));
  SEXP neg_loglik_ = PROTECT(ScalarReal(neg_loglik));
  SEXP out = PROTECT(allocVector(VECSXP, 4));
  SEXP names = PROTECT(allocVector(STRSXP, 4));

  SET_VECTOR_ELT(out, 0, chol_);
  SET_VECTOR_ELT(out, 1, alpha_);
  SET_VECTOR_ELT(out, 2, jitter_);
  SET_VECTOR_ELT(out, 3, neg_loglik_);

  SET_STRING_ELT(names, 0, mkChar("chol"));
  SET_STRING_ELT(names, 1, mkChar("alpha"));
  SET_STRING_ELT(names, 2, mkChar("jitter"));
  SET_STRING_ELT(names, 3, mkChar("neg_loglik"));
  setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(6);
  return out;
}

SEXP DGrowthR_exact_gaussian_fit_state(SEXP X_, SEXP y_, SEXP lengthscale_,
                                       SEXP sigma2_, SEXP nugget_) {
  if (!isReal(X_) || !isMatrix(X_)) {
    error("X must be a numeric matrix.");
  }
  SEXP dim = getAttrib(X_, R_DimSymbol);
  int n = INTEGER(dim)[0];
  int p = INTEGER(dim)[1];

  const double *X = REAL(X_);
  return fit_state_from_covariance_builder(y_, lengthscale_, sigma2_,
                                           nugget_, n, p,
                                           fill_gaussian_covariance, X);
}

SEXP DGrowthR_exact_gaussian_fit_state_dist(SEXP dist_sq_, SEXP y_,
                                            SEXP lengthscale_, SEXP sigma2_,
                                            SEXP nugget_) {
  if (!isReal(dist_sq_)) {
    error("dist_sq must be numeric.");
  }

  SEXP dim = getAttrib(dist_sq_, R_DimSymbol);
  if (dim == R_NilValue) {
    error("dist_sq must be a matrix or array.");
  }

  int ndim = LENGTH(dim);
  if (ndim != 2 && ndim != 3) {
    error("dist_sq must be a matrix or 3-dimensional array.");
  }

  int n = INTEGER(dim)[0];
  if (INTEGER(dim)[1] != n) {
    error("dist_sq must be square in its first two dimensions.");
  }
  int p = (ndim == 2) ? 1 : INTEGER(dim)[2];

  const double *dist_sq = REAL(dist_sq_);
  return fit_state_from_covariance_builder(y_, lengthscale_, sigma2_,
                                           nugget_, n, p,
                                           fill_gaussian_covariance_dist,
                                           dist_sq);
}

SEXP DGrowthR_exact_gaussian_lbfgsb_dist(SEXP start_, SEXP lower_, SEXP upper_,
                                         SEXP dist_sq_, SEXP y_,
                                         SEXP nugget_weights_, SEXP maxit_) {
  if (!isReal(start_) || !isReal(lower_) || !isReal(upper_)) {
    error("start, lower, and upper must be numeric.");
  }
  if (!isReal(dist_sq_)) {
    error("dist_sq must be numeric.");
  }
  if (!isReal(y_)) {
    error("y_centered must be numeric.");
  }
  if (!isReal(nugget_weights_)) {
    error("nugget_weights must be numeric.");
  }

  SEXP dim = getAttrib(dist_sq_, R_DimSymbol);
  if (dim == R_NilValue) {
    error("dist_sq must be a matrix or 3-dimensional array.");
  }

  int ndim = LENGTH(dim);
  if (ndim != 2 && ndim != 3) {
    error("dist_sq must be a matrix or 3-dimensional array.");
  }

  int n = INTEGER(dim)[0];
  if (INTEGER(dim)[1] != n) {
    error("dist_sq must be square in its first two dimensions.");
  }
  int p = (ndim == 2) ? 1 : INTEGER(dim)[2];
  int npar = p + 1;

  if (LENGTH(start_) != npar || LENGTH(lower_) != npar || LENGTH(upper_) != npar) {
    error("start, lower, and upper must have length ncol(X) + 1.");
  }
  if (LENGTH(y_) != n) {
    error("Length of y_centered must match covariance dimension.");
  }

  int nugget_len = LENGTH(nugget_weights_);
  if (nugget_len != 1 && nugget_len != n) {
    error("nugget_weights must have length 1 or covariance dimension.");
  }

  for (int i = 0; i < nugget_len; ++i) {
    if (!R_FINITE(REAL(nugget_weights_)[i]) || REAL(nugget_weights_)[i] <= 0.0) {
      error("All nugget weights must be positive and finite.");
    }
  }

  ExactProfileData data;
  data.dist_sq = REAL(dist_sq_);
  data.y = REAL(y_);
  data.nugget_weights = REAL(nugget_weights_);
  data.n = n;
  data.p = p;
  data.nugget_len = nugget_len;

  double *x = Calloc(npar, double);
  double *lower = Calloc(npar, double);
  double *upper = Calloc(npar, double);
  int *nbd = Calloc(npar, int);

  for (int i = 0; i < npar; ++i) {
    x[i] = REAL(start_)[i];
    lower[i] = REAL(lower_)[i];
    upper[i] = REAL(upper_)[i];
    nbd[i] = 2;
  }

  double fmin = 0.0;
  int fail = 0;
  int fncount = 0;
  int grcount = 0;
  int maxit = asInteger(maxit_);
  if (maxit <= 0) {
    maxit = 100;
  }
  char msg[60];
  msg[0] = '\0';

  lbfgsb(npar, 5, x, lower, upper, nbd, &fmin,
         exact_profile_lbfgsb_fn, exact_profile_lbfgsb_gr, &fail,
         &data, 1e7, 0.0, &fncount, &grcount, maxit, msg, 0, 10);

  SEXP par_ = PROTECT(allocVector(REALSXP, npar));
  for (int i = 0; i < npar; ++i) {
    REAL(par_)[i] = x[i];
  }

  SEXP value_ = PROTECT(ScalarReal(fmin));
  SEXP conv_ = PROTECT(ScalarInteger(fail));
  SEXP counts_ = PROTECT(allocVector(INTSXP, 2));
  INTEGER(counts_)[0] = fncount;
  INTEGER(counts_)[1] = grcount;
  SEXP msg_ = PROTECT(mkString(msg));

  SEXP out = PROTECT(allocVector(VECSXP, 5));
  SEXP names = PROTECT(allocVector(STRSXP, 5));

  SET_VECTOR_ELT(out, 0, par_);
  SET_VECTOR_ELT(out, 1, value_);
  SET_VECTOR_ELT(out, 2, conv_);
  SET_VECTOR_ELT(out, 3, counts_);
  SET_VECTOR_ELT(out, 4, msg_);

  SET_STRING_ELT(names, 0, mkChar("par"));
  SET_STRING_ELT(names, 1, mkChar("value"));
  SET_STRING_ELT(names, 2, mkChar("convergence"));
  SET_STRING_ELT(names, 3, mkChar("counts"));
  SET_STRING_ELT(names, 4, mkChar("message"));
  setAttrib(out, R_NamesSymbol, names);

  Free(x);
  Free(lower);
  Free(upper);
  Free(nbd);

  UNPROTECT(7);
  return out;
}
