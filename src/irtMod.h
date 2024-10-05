#ifndef irtMod_H
#define irtMod_H
#include "grades.h"
#include "times.h"
#include "exams.h"
#include "extractParams.h"
#include "latent.h"



class GRTC_MOD
{
private:
  Eigen::VectorXd _theta;
  Eigen::VectorXd _exams_grades;
  Eigen::VectorXd _exams_days;
  Eigen::VectorXd _exams_set;
  Eigen::VectorXd _exams_obsflag;
  unsigned int _max_day;
  unsigned int _n_grades;
  unsigned int _n_exams;
  unsigned int _dim_irt;
  bool _rotated;

public:
  GRTC_MOD(Eigen::VectorXd THETA,
          Eigen::VectorXd EXAMS_GRADES,
          Eigen::VectorXd EXAMS_DAYS,
          Eigen::VectorXd EXAMS_SET,
          Eigen::VectorXd EXAMS_OBSFLAG,
          const unsigned int MAX_DAY,
          const unsigned int N_GRADES,
          const unsigned int N_EXAMS,
          const bool ROTATED):
  _theta(THETA),
  _exams_grades(EXAMS_GRADES),
  _exams_days(EXAMS_DAYS),
  _exams_set(EXAMS_SET),
  _exams_obsflag(EXAMS_OBSFLAG),
  _max_day(MAX_DAY),
  _n_grades(N_GRADES),
  _n_exams(N_EXAMS),
  _rotated(ROTATED)
  {
    _dim_irt = N_EXAMS*(N_GRADES+3);
  }

  // conditional log-likelihood
  double ll(const double ABILITY, const double SPEED);

  //gradient of conditional log-likelihood
  Eigen::VectorXd grll(const double ABILITY, const double SPEED);

  // complete log-likelihood
  double cll(const double ABILITY, const double SPEED);

};

double GRTC_MOD::ll(const double ABILITY, const double SPEED) {

  double out = 0;

  for(unsigned int exam = 0; exam < _n_exams; exam++){

    if(_exams_set[exam]){
      out += examLik(exam,
                     _exams_grades(exam),
                     _exams_days(exam),
                     _max_day,
                     _exams_obsflag(exam),
                     _theta,
                     _n_grades,
                     _n_exams,
                     ABILITY, SPEED, 1);

    }
  }

  return out;
}
double GRTC_MOD::cll(const double ABILITY, const double SPEED) {

  LAT_DISTR lat(_theta.segment(_dim_irt,2));
  double out = lat.ll(ABILITY, SPEED);

  for(unsigned int exam = 0; exam < _n_exams; exam++){

    if(_exams_set[exam]){
      out += examLik(exam,
                     _exams_grades(exam),
                     _exams_days(exam),
                     _max_day,
                     _exams_obsflag(exam),
                     _theta,
                     _n_grades,
                     _n_exams,
                     ABILITY, SPEED, 1);

    }
  }

  return out;
}

Eigen::VectorXd GRTC_MOD::grll(const double ABILITY, const double SPEED){

  Eigen::VectorXd gr = Eigen::VectorXd::Zero(_theta.size());
  Eigen::VectorXd gr_irt = Eigen::VectorXd::Zero(_dim_irt+2);

  for(unsigned int exam = 0; exam < _n_exams; exam++){

    if(_exams_set[exam]){
      gr_irt += grl_examLik(exam,
                            _exams_grades(exam),
                             _exams_days(exam),
                             _max_day,
                             _exams_obsflag(exam),
                             _theta,
                             _n_grades,
                             _n_exams,
                             ABILITY, SPEED,
                             _rotated);
    }
  }

  gr.segment(0, _dim_irt+2) = gr_irt;
  return gr;
}

//' Log-likelihood and gradient of the conditional GRTC Model on one observation
//'
//' Used for internal testing
//'
//' @param THETA Parameter vector
//' @param EXAMS_GRADES Vector of exam grades
//' @param EXAMS_DAYS Vector of Eexam days
//' @param EXAMS_SET Vector of binary values. 1 to include an exam in the study-plan, 0 otherwise
//' @param EXAMS_OBSFLAG Vector of binary values. 1 if the exam is observed, 0 otherwise
//' @param MAX_DAY Maximum day of observation
//' @param N_GRADES Number of possible grades modelled
//' @param N_EXAMS Number of possible exams modelled
//' @param ABILITY Ability latent parameter
//' @param SPEED Speed latent parameter
//' @param ROTATE TRUE to rotate latent scores using their covariance matrix
//' @param GRFLAG Set to true to return the gradient along the log-likelihood
//'
//' @export
// [[Rcpp::export]]
Rcpp::List conditional_igrtcm(
    Eigen::VectorXd& THETA,
    Eigen::VectorXd& EXAMS_GRADES,
    Eigen::VectorXd& EXAMS_DAYS,
    Eigen::VectorXd& EXAMS_SET,
    Eigen::VectorXd& EXAMS_OBSFLAG,
    const unsigned int MAX_DAY,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const double ABILITY,
    double SPEED,
    const bool ROTATE,
    const bool GRFLAG = true
){

  const unsigned int dim_irt = N_EXAMS*(N_GRADES+3);
  double speed = SPEED;

  if(ROTATE){
    speed = THETA(dim_irt)*ABILITY + THETA(dim_irt+1)*SPEED;
  }

  GRTC_MOD irt_mod(THETA,
                   EXAMS_GRADES,
                   EXAMS_DAYS,
                   EXAMS_SET,
                   EXAMS_OBSFLAG,
                   MAX_DAY,
                   N_GRADES,
                   N_EXAMS,
                   ROTATE);

  double ll = irt_mod.ll(ABILITY, speed);
  Eigen::VectorXd gr = Eigen::VectorXd::Zero(dim_irt+2);
  if(GRFLAG){
    gr = irt_mod.grll(ABILITY, speed);
  }

  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("gr") = gr,
      Rcpp::Named("ll") = ll
    );

  return output;
}

//' Log-likelihood and gradient of the complete GRTC Model on one observation
//'
//' Used for internal testing
//'
//' @param THETA Parameter vector
//' @param EXAMS_GRADES Vector of exam grades
//' @param EXAMS_DAYS Vector of exam days
//' @param EXAMS_SET Vector of binary values. 1 to include an exam in the study-plan, 0 otherwise
//' @param EXAMS_OBSFLAG Vector of binary values. 1 if the exam is observed, 0 otherwise
//' @param MAX_DAY Maximum day of observation
//' @param N_GRADES Number of possible grades modelled
//' @param N_EXAMS Number of possible exams modelled
//' @param ABILITY Ability latent parameter
//' @param SPEED Speed latent parameter
//'
//' @export
// [[Rcpp::export]]
double complete_igrtcm(
     Eigen::VectorXd& THETA,
     Eigen::VectorXd& EXAMS_GRADES,
     Eigen::VectorXd& EXAMS_DAYS,
     Eigen::VectorXd& EXAMS_SET,
     Eigen::VectorXd& EXAMS_OBSFLAG,
     const unsigned int MAX_DAY,
     const unsigned int N_GRADES,
     const unsigned int N_EXAMS,
     const double ABILITY,
     double SPEED
 ){

   const unsigned int dim_irt = N_EXAMS*(N_GRADES+3);

   GRTC_MOD irt_mod(THETA,
                    EXAMS_GRADES,
                    EXAMS_DAYS,
                    EXAMS_SET,
                    EXAMS_OBSFLAG,
                    MAX_DAY,
                    N_GRADES,
                    N_EXAMS,
                    false);

   double cll = irt_mod.cll(ABILITY, SPEED);


   return cll;
 }











#endif
