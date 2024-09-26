#ifndef thresholds_H
#define thresholds_H

//' Intercepts reparameterisation
//'
//' @param X Intercepts or unconstrained parameter values for grades low-to-high.
//' @param CON2UN TRUE if going from the constrained space to the unconstrained one.
//'
//' @returns It allows to go back and forth from constrained intercepts
//' to unconstrained parameters.
//'
//' @export
// [[Rcpp::export]]
Eigen::VectorXd reparThr(const Eigen::VectorXd& X, bool CON2UN = false) {
  const int n_thr = X.size();
  Eigen::VectorXd out(n_thr);

  if (!CON2UN) {
    // From unconstrained to constrained
    out = X;
    out.segment(1, n_thr-1) = X.segment(1, n_thr-1).array().exp();
    for (unsigned int i = 1; i < n_thr; i++) {
      out[i] += out[i - 1];
    }
    return out;

  } else {

    // From constrained to unconstrained
    out = X;
    for (unsigned int i = 1; i < n_thr; i++) {
      out[i] = log(X[i] - X[i - 1]);
    }

    return out;

  }

 }







#endif
