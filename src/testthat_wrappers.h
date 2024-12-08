#ifndef testthat_wrappers_H
#define testthat_wrappers_H
#include "grades.h"
#include "times.h"
#include "exams.h"
#include "extractParams.h"
#include "latent.h"
#include "cr.h"
#include "conditional_models.h"
#include "em.h"

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
Rcpp::List cpp_hazard(
    const unsigned int OUTCOME,
    const unsigned int YEAR,
    Eigen::Map<Eigen::VectorXd>  THETA,
    Eigen::Map<Eigen::VectorXd>  COVARIATES,
    const unsigned int YB,
    const bool LOGFLAG = false
){
  double prob = cr::hazard(OUTCOME,
                           YEAR,
                           THETA,
                           COVARIATES,
                           YB,
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

// Latent
//' @export
// [[Rcpp::export]]
Rcpp::List cpp_lat(
    Eigen::Map<Eigen::VectorXd> THETA,
    const double ABILITY,
    const double SPEED,
    const unsigned int DIM_IRT,
    const bool ROTATED,
    const bool GRFLAG
){

  LAT lat(THETA, DIM_IRT, ROTATED);
  double ll = lat.ll(ABILITY, SPEED);

  Eigen::VectorXd grll = Eigen::VectorXd::Zero(THETA.size());

  if(GRFLAG){
    grll = lat.grll(ABILITY, SPEED);
  }

  Rcpp::List output = Rcpp::List::create(
    // Rcpp::Named("L")=L,
    Rcpp::Named("ll")=ll,
    Rcpp::Named("grll")=grll
  );

  return output;
}

// Expectation Maximization
//' @export
// [[Rcpp::export]]
Rcpp::List cpp_estep(
    Eigen::VectorXd THETA,
    Eigen::MatrixXd EXAMS_GRADES,
    Eigen::MatrixXd EXAMS_DAYS,
    Eigen::MatrixXd EXAMS_SET,
    Eigen::MatrixXd EXAMS_OBSFLAG,
    Eigen::VectorXd MAX_DAY,
    Eigen::VectorXd OUTCOME,
    Eigen::MatrixXd EXT_COVARIATES,
    Eigen::VectorXd YEAR_FIRST,
    Eigen::VectorXd YEAR_LAST,
    Eigen::VectorXd YEAR_LAST_EXAM,
    Eigen::MatrixXd GRID,
    Eigen::VectorXd WEIGHTS,
    const unsigned int YB,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const std::string MOD
){

  const int m = EXT_COVARIATES.cols();
  const int dim_cr = 2*(YB+m+2)+1;
  const int dim_irt = THETA.size()-dim_cr-2;
  Eigen::MatrixXd L{{1,0},{THETA(dim_irt), THETA(dim_irt+1)}};
  GRID = GRID * L.transpose();

  Eigen::MatrixXd Ew = EM::Estep(THETA,
                                 EXAMS_GRADES,
                                 EXAMS_DAYS,
                                 EXAMS_SET,
                                 EXAMS_OBSFLAG,
                                 MAX_DAY,
                                 OUTCOME,
                                 EXT_COVARIATES,
                                 YEAR_FIRST,
                                 YEAR_LAST,
                                 YEAR_LAST_EXAM,
                                 GRID,
                                 WEIGHTS,
                                 YB,
                                 N_GRADES,
                                 N_EXAMS,
                                 MOD);




  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("Ew")=Ew,
    Rcpp::Named("nodes")=GRID

  );
  return output;

}

// Expectation Maximization
//' @export
// [[Rcpp::export]]
Rcpp::List cpp_mstep(
    Eigen::VectorXd THETA,
    Eigen::MatrixXd EXAMS_GRADES,
    Eigen::MatrixXd EXAMS_DAYS,
    Eigen::MatrixXd EXAMS_SET,
    Eigen::MatrixXd EXAMS_OBSFLAG,
    Eigen::VectorXd MAX_DAY,
    Eigen::VectorXd OUTCOME,
    Eigen::MatrixXd EXT_COVARIATES,
    Eigen::VectorXd YEAR_FIRST,
    Eigen::VectorXd YEAR_LAST,
    Eigen::VectorXd YEAR_LAST_EXAM,
    Eigen::MatrixXd GRID,
    Eigen::VectorXd EWEIGHTS,
    const unsigned int YB,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const std::string MOD
){

  const int m = EXT_COVARIATES.cols();
  const int dim_cr = 2*(YB+m+2)+1;
  const int dim_irt = THETA.size()-dim_cr-2;



  Eigen::VectorXd gr = Eigen::VectorXd::Zero(THETA.size());

  EM::Ejnll eclass(
      EXAMS_GRADES,
      EXAMS_DAYS,
      EXAMS_SET,
      EXAMS_OBSFLAG,
      MAX_DAY,
      OUTCOME,
      EXT_COVARIATES,
      YEAR_FIRST,
      YEAR_LAST,
      YEAR_LAST_EXAM,
      YB,
      N_GRADES,
      N_EXAMS,
      MOD
  );

  eclass.update_quadrature(GRID, EWEIGHTS);

  double nll = eclass.f_grad(THETA, gr);
  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("nll")=nll,
    Rcpp::Named("gr")=gr
  );
  return output;

}

// Expectation Maximization
//' @export
// [[Rcpp::export]]
Rcpp::List cpp_mstep2(
    Eigen::VectorXd THETA,
    Eigen::MatrixXd EXAMS_GRADES,
    Eigen::MatrixXd EXAMS_DAYS,
    Eigen::MatrixXd EXAMS_SET,
    Eigen::MatrixXd EXAMS_OBSFLAG,
    Eigen::VectorXd MAX_DAY,
    Eigen::VectorXd OUTCOME,
    Eigen::MatrixXd EXT_COVARIATES,
    Eigen::VectorXd YEAR_FIRST,
    Eigen::VectorXd YEAR_LAST,
    Eigen::VectorXd YEAR_LAST_EXAM,
    Eigen::MatrixXd GRID,
    Eigen::VectorXd EWEIGHTS,
    const unsigned int YB,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const std::string MOD
){

  const int m = EXT_COVARIATES.cols();
  const int dim_cr = 2*(YB+m+2)+1;
  const int dim_irt = THETA.size()-dim_cr-2;


  double irt_nll = 0;
  for(unsigned int i = 0; i < EXAMS_GRADES.rows(); i++){

    // Initialize conditional IRT model
    GRTC_MOD irt_mod(THETA,
                     EXAMS_GRADES.row(i),
                     EXAMS_DAYS.row(i),
                     EXAMS_SET.row(i),
                     EXAMS_OBSFLAG.row(i),
                     MAX_DAY(i),
                     N_GRADES,
                     N_EXAMS,
                     false);

    irt_nll-=irt_mod.cll(GRID(0,0), GRID(0,1));
  }


  Eigen::VectorXd gr = Eigen::VectorXd::Zero(THETA.size());

  EM::EAPLOGJ eclass(
      EXAMS_GRADES,
      EXAMS_DAYS,
      EXAMS_SET,
      EXAMS_OBSFLAG,
      MAX_DAY,
      OUTCOME,
      EXT_COVARIATES,
      YEAR_FIRST,
      YEAR_LAST,
      YEAR_LAST_EXAM,
      YB,
      N_GRADES,
      N_EXAMS,
      MOD
  );

  eclass.update_quadrature(GRID, EWEIGHTS);

  double nll = eclass.irt_cll(THETA, GRID(0,0), GRID(0,1));
  double nllf = eclass.f_grad(THETA, gr);

  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("gr")=gr,
    Rcpp::Named("grid")=eclass._grid,
    Rcpp::Named("nll")=nll,
    Rcpp::Named("nllf")=nllf,
    Rcpp::Named("irt_nll")=irt_nll

  );
  return output;

}
#endif
