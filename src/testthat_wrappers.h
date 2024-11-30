#ifndef testthat_wrappers_H
#define testthat_wrappers_H
#include "grades.h"
#include "times.h"
#include "exams.h"
#include "extractParams.h"
#include "latent.h"
#include "cr.h"
#include "conditional_models.h"

// EXAMS //
//' @export
// [[Rcpp::export]]
Rcpp::List cpp_pGreaterGrades(
    const unsigned int GRADE,
    const unsigned int EXAM,
    Eigen::Map<Eigen::VectorXd> THETA,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const double ABILITY,
    const bool LOGFLAG){

  double prob = exams::pGreaterGrades(GRADE, EXAM, THETA, N_GRADES, N_EXAMS, ABILITY, LOGFLAG);
  Eigen::VectorXd gr=exams::grad::gr_pGreaterGrades(GRADE, EXAM, THETA, N_GRADES, N_EXAMS, ABILITY);

  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("prob")=prob,
    Rcpp::Named("gr")=gr
  );

  return output;
}

//' @export
// [[Rcpp::export]]
Rcpp::List cpp_pGrade(
    const unsigned int GRADE,
    const unsigned int EXAM,
    Eigen::Map<Eigen::VectorXd> THETA,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const double ABILITY,
    const bool LOGFLAG){

  double prob = exams::pGrade(GRADE, EXAM, THETA, N_GRADES, N_EXAMS, ABILITY, LOGFLAG);
  Eigen::VectorXd gr=exams::grad::gr_pGrade(GRADE, EXAM, THETA, N_GRADES, N_EXAMS, ABILITY);

  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("prob")=prob,
    Rcpp::Named("gr")=gr
  );

  return output;
}

//' @export
// [[Rcpp::export]]
Rcpp::List cpp_pTimeExam(
    const unsigned int EXAM,
    const double DAY,
    Eigen::Map<Eigen::VectorXd> THETA,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const double SPEED,
    const double ABILITY,
    const bool CDFFLAG,
    const bool ROTATED,
    const bool LOGFLAG){

  double prob=exams::pTimeExam(EXAM, DAY, THETA, N_GRADES, N_EXAMS, SPEED, CDFFLAG, LOGFLAG);
  Eigen::VectorXd gr=exams::grad::gr_pTimeExam(EXAM, DAY, THETA, N_GRADES, N_EXAMS, SPEED, ABILITY, CDFFLAG, ROTATED, LOGFLAG);

  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("prob")=prob,
    Rcpp::Named("gr")=gr
  );

  return output;
}

//' @export
// [[Rcpp::export]]
Rcpp::List cpp_examLik(
    const unsigned int EXAM,
    const unsigned int GRADE,
    const double DAY,
    const double MAX_DAY,
    const bool OBSFLAG,
    Eigen::Map<Eigen::VectorXd> THETA,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const double ABILITY,
    const double SPEED,
    const bool ROTATED
){
  double ll=exams::examLik(EXAM, GRADE, DAY, MAX_DAY, OBSFLAG, THETA, N_GRADES, N_EXAMS, ABILITY, SPEED, true);
  Eigen::VectorXd grll=exams::grad::grl_examLik(EXAM, GRADE, DAY, MAX_DAY, OBSFLAG, THETA, N_GRADES, N_EXAMS, ABILITY, SPEED, ROTATED);

  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("ll")=ll,
    Rcpp::Named("grll")=grll
  );

  return output;
}



// COMPETING RISK //
//' @export
// [[Rcpp::export]]
Rcpp::List cpp_survival(
    const unsigned int YEAR_FIRST,
    const unsigned int YEAR_LAST,
    Eigen::Map<Eigen::VectorXd>  THETA,
    Eigen::Map<Eigen::VectorXd>  COVARIATES,
    const unsigned int YB,
    const unsigned int YEAR_LAST_EXAM,
    const bool LOGFLAG = false
){
  double prob = cr::survival(YEAR_FIRST,
                             YEAR_LAST,
                         THETA,
                         COVARIATES,
                         YB,
                         YEAR_LAST_EXAM,
                         LOGFLAG);
  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("prob")=prob
  );

  return output;
}

//' @export
// [[Rcpp::export]]
Rcpp::List cpp_outcome(
    const unsigned int OUTCOME,
    const unsigned int YEAR_FIRST,
    const unsigned int YEAR_LAST,
    Eigen::Map<Eigen::VectorXd>  THETA,
    Eigen::Map<Eigen::VectorXd>  COVARIATES,
    const unsigned int YB,
    const unsigned int YEAR_LAST_EXAM,
    const bool LOGFLAG = false
){
  double prob = cr::outcomeLik(OUTCOME,
                               YEAR_FIRST,
                               YEAR_LAST,
                               THETA,
                               COVARIATES,
                               YB,
                               YEAR_LAST_EXAM,
                               LOGFLAG);

  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("prob")=prob
  );

  return output;
}
//' @export
// [[Rcpp::export]]
Rcpp::List cpp_crmod(
    const unsigned int OUTCOME,
    const unsigned int YEAR_FIRST,
    const unsigned int YEAR_LAST,
    Eigen::Map<Eigen::VectorXd> THETA,
    Eigen::Map<Eigen::VectorXd> EXT_COVARIATES,
    const unsigned int YB,
    const unsigned int YEAR_LAST_EXAM,
    Eigen::Map<Eigen::VectorXd> LAT_POINTS,
    const bool ROTATE
){

  const int q = EXT_COVARIATES.size();
  const int dim_cr = 2*(YB+q+2)+1;
  const int dim_irt = THETA.size()-dim_cr-2;

  CR_MOD cr(THETA, OUTCOME, EXT_COVARIATES, YB, YEAR_FIRST, YEAR_LAST, YEAR_LAST_EXAM, ROTATE);
  Eigen::VectorXd lat(2);
  if(ROTATE){
    Eigen::MatrixXd L{{1,0},{THETA(dim_irt), THETA(dim_irt+1)}};
    lat = L*LAT_POINTS;
  }else{
    lat = LAT_POINTS;
  }

  double ll = cr.ll(lat(0), lat(1));
  Eigen::VectorXd grll=cr.grll(lat(0), lat(1));

  Rcpp::List output = Rcpp::List::create(
    // Rcpp::Named("L")=L,
    Rcpp::Named("ll")=ll,
    Rcpp::Named("grll")=grll
  );

  return output;
}


#endif
