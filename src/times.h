#ifndef times_H
#define times_H
#include "extractParams.h"

//' Evaluate the c.d.f or p.d.f of the last attempt to an exam
//'
//' @param EXAM Exam of interest. Integers from 0 to N_EXAMS -1.
//' @param DAY Day of interest.
//' @param THETA Portion of the parameter vector related to the IRT model
//' @param N_GRADES Number of grades modelled.
//' @param N_EXAMS Number of exams.
//' @param SPEED speed value.
//' @param CDFFLAG `TRUE` for c.d.f. of time. `FALSE` for p.d.f.
//' @param LOGFLAG Set TRUE to return log value.
//'
//' @export
// [[Rcpp::export]]
double pTimeExam(
    const unsigned int EXAM,
    const double DAY,
    Eigen::VectorXd& THETA,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const double SPEED,
    const bool CDFFLAG,
    const bool LOGFLAG = false
){
  std::vector<double> pars(2);
  pars[0] = extract_params_irt(THETA, N_GRADES, N_EXAMS, 3, EXAM)(0);
  pars[1] = extract_params_irt(THETA, N_GRADES, N_EXAMS, 4, EXAM)(0);
  const double mean = pars[0] - SPEED;
  const double sd = 1/pars[1];
  double out;
  if(CDFFLAG){
    out = R::plnorm(DAY, mean, sd, true, LOGFLAG);
  }else{
    out = R::dlnorm(DAY, mean, sd, LOGFLAG);
  }

  return(out);
}


// //' @export
// // [[Rcpp::export]]
// Rcpp::List internal_pTimeExam(
//     const unsigned int EXAM,
//     const double DAY,
//     Eigen::VectorXd& THETA,
//     const unsigned int N_GRADES,
//     const unsigned int N_EXAMS,
//     const double SPEED,
//     const bool CDFFLAG,
//     const bool LOGFLAG = false
// ){
//   std::vector<double> pars(2);
//   pars[0] = extract_params_irt(THETA, N_GRADES, N_EXAMS, 3, EXAM)(0);
//   pars[1] = extract_params_irt(THETA, N_GRADES, N_EXAMS, 4, EXAM)(0);
//   const double mean = pars[0] - SPEED;
//   const double sd = 1/pars[1];
//   double out;
//   if(CDFFLAG){
//     out = R::plnorm(DAY, mean, sd, true, LOGFLAG);
//   }else{
//     out = R::dlnorm(DAY, mean, sd, LOGFLAG);
//   }
//
//   Rcpp::List output = Rcpp::List::create(
//     Rcpp::Named("loc") = pars[0],
//     Rcpp::Named("invsd") = pars[1],
//     Rcpp::Named("mean") = mean,
//     Rcpp::Named("sd") = sd,
//     Rcpp::Named("out") = out
//   );
//   return(output);
// }
//' Evaluate the gradient of the c.d.f or p.d.f of the last attempt to an exam
//'
//' @param EXAM Exam of interest. Integers from 0 to N_EXAMS -1.
//' @param DAY Day of interest.
//' @param THETA Portion of the parameter vector related to the IRT model
//' @param N_GRADES Number of grades modelled.
//' @param N_EXAMS Number of exams.
//' @param SPEED speed value.
//' @param ABILITY ability value.
//' @param CDFFLAG `TRUE` for c.d.f. of time. `FALSE` for p.d.f.
//' @param ROTATED Have latent scores been rotated using their variance?
//' @param LOGFLAG TRUE to compute the gradient of the log density.(not available for cdf)
//'
//' @export
// [[Rcpp::export]]
Eigen::VectorXd gr_pTimeExam(
     const unsigned int EXAM,
     const double DAY,
     Eigen::VectorXd& THETA,
     const unsigned int N_GRADES,
     const unsigned int N_EXAMS,
     const double SPEED,
     const double ABILITY,
     const bool CDFFLAG,
     const bool ROTATED,
     const bool LOGFLAG = false
){
  const unsigned int dim_irt = N_EXAMS*(N_GRADES+3);

  Eigen::VectorXd gr = Eigen::VectorXd::Zero(dim_irt+2);
  std::vector<double> pars(2);
  pars[0] = extract_params_irt(THETA, N_GRADES, N_EXAMS, 3, EXAM)(0);
  pars[1] = extract_params_irt(THETA, N_GRADES, N_EXAMS, 4, EXAM)(0);
  const double mean = pars[0] - SPEED;
  const double sd = 1/pars[1];


  std::vector<unsigned int> idx_m = extract_params_idx_irt(THETA, N_GRADES, N_EXAMS, 3, EXAM);
  std::vector<unsigned int> idx_v = extract_params_idx_irt(THETA, N_GRADES, N_EXAMS, 4, EXAM);

  if(CDFFLAG){
    const double tmp = R::dnorm(pars[1]*(log(DAY)-mean), 0, 1, false);
    gr(idx_m[0]) = -tmp*pars[1];
    gr(idx_v[0]) = tmp*(log(DAY)-mean)*pars[1];
    if(ROTATED){
      double speed_raw = (SPEED -THETA(dim_irt)*ABILITY)/THETA(dim_irt+1);
      gr(dim_irt) = tmp*pars[1]*ABILITY;
      gr(dim_irt+1) = tmp*pars[1]*speed_raw;
    }
  }else{
    if(!LOGFLAG){
      const double tmp = R::dlnorm(DAY, mean, sd, false);
      gr(idx_m[0]) = tmp*pow(pars[1],2)*(log(DAY)-mean);
      gr(idx_v[0]) = (tmp/pars[1] - tmp * pars[1] * pow(log(DAY)-mean, 2))*pars[1];
      if(ROTATED){
        double speed_raw = (SPEED -THETA(dim_irt)*ABILITY)/THETA(dim_irt+1);
        gr(dim_irt) =  -tmp * pow(pars[1],2)*(log(DAY)-mean)*ABILITY;
        gr(dim_irt+1) =  -tmp * pow(pars[1],2)*(log(DAY)-mean)*speed_raw;
      }
    }else{
      gr(idx_v[0]) = 1 - pow(pars[1],2) * pow(log(DAY)-mean, 2);
      gr(idx_m[0]) = pow(pars[1],2)*(log(DAY)-mean);
      if(ROTATED){
        double speed_raw = (SPEED -THETA(dim_irt)*ABILITY)/THETA(dim_irt+1);
        gr(dim_irt) =  -pow(pars[1],2)*(log(DAY)-mean)*ABILITY;
        gr(dim_irt+1) =  -pow(pars[1],2)*(log(DAY)-mean)*speed_raw;
      }
    }

  }

  return(gr);
}

#endif
