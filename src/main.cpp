#include <Rcpp.h>
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
#include "extractParams.h"
#include "latent.h"
#include "irtMod.h"

//' GRTCM computations with Gaussian quadrature
//'
//' Compute log-likelihood an gradient of the Graded response time-censored model
//'
//' @param THETA Parameter vector
//' @param EXAMS_GRADES Matrix of exam grades
//' @param EXAMS_DAYS Matrix of exam days
//' @param EXAMS_SET Matrix of binary values. 1 to include an exam in the study-plan, 0 otherwise
//' @param EXAMS_OBSFLAG Matrix of binary values. 1 if the exam is observed, 0 otherwise
//' @param MAX_DAY Vector with maximum day per observation
//' @param GRID Grid of quadrature points
//' @param WEIGHTS Weights of quadrature points
//' @param N_GRADES Number of possible grades modelled
//' @param N_EXAMS Number of possible exams modelled
//' @param ROTGRID TRUE to rotate the quadrature grid using the latent covariance matrix
//' @param GRFLAG Set to true to return the gradient along the log-likelihood
//'
//' @export
// [[Rcpp::export]]
Rcpp::List GRTCM_GH(
    Eigen::VectorXd& THETA,
    Eigen::MatrixXd& EXAMS_GRADES,
    Eigen::MatrixXd& EXAMS_DAYS,
    Eigen::MatrixXd& EXAMS_SET,
    Eigen::MatrixXd& EXAMS_OBSFLAG,
    Eigen::VectorXd& MAX_DAY,
    Eigen::MatrixXd GRID,
    Eigen::VectorXd WEIGHTS,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const bool GRFLAG = true,
    const bool ROTGRID = true
){
  double ll = 0;
  Eigen::VectorXd grll = Eigen::VectorXd::Zero(THETA.size());

  const unsigned int n = EXAMS_GRADES.rows();
  const unsigned int nq = GRID.rows();
  const unsigned int dim_irt = N_EXAMS*(N_GRADES+3);


  if(ROTGRID){
    Eigen::MatrixXd L{{1,0},{THETA(dim_irt), THETA(dim_irt+1)}};
    GRID = GRID * L.transpose();
  }

  Eigen::MatrixXd llMat(n,nq);
  for(unsigned int i = 0; i < n; i++){

    // Initialize conditional IRT model
    GRTC_MOD irt_mod(THETA,
                     EXAMS_GRADES.row(i),
                     EXAMS_DAYS.row(i),
                     EXAMS_SET.row(i),
                     EXAMS_OBSFLAG.row(i),
                     MAX_DAY(i),
                     N_GRADES,
                     N_EXAMS,
                     ROTGRID);

    Eigen::VectorXd f(nq);
    Eigen::MatrixXd gr = Eigen::MatrixXd::Zero(THETA.size(), nq);

    for(unsigned int point = 0; point < nq; point++){
      f(point) = exp(irt_mod.ll(GRID(point, 0), GRID(point, 1)));
      llMat(i,point)=(irt_mod.ll(GRID(point, 0), GRID(point, 1)));

      if(GRFLAG){
        Eigen::VectorXd gr_point = irt_mod.grll(GRID(point, 0), GRID(point, 1));
        gr_point *= f(point);
        gr.col(point) = gr_point;
      }

    }
    double lli = std::max(-10000.0, log(f.dot(WEIGHTS)));

    ll += lli;
    if(GRFLAG){
      grll += gr*WEIGHTS/exp(lli);
    }
  }

  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("llMat") = llMat,
      Rcpp::Named("gr") = grll,
      Rcpp::Named("ll") = ll
    );

  return output;

}
