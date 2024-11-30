#ifndef condMods_H
#define condMods_H
#include "grades.h"
#include "times.h"
#include "exams.h"
#include "extractParams.h"
#include "latent.h"
#include "cr.h"

// COMPETING RISK CONDITIONAL MODEL
class CR_MOD
{
private:
  Eigen::VectorXd _theta;
  unsigned int _outcome;
  Eigen::VectorXd _ext_covariates;
  unsigned int _yb;
  unsigned int _year_first;
  unsigned int _year_last;
  unsigned int _year_last_exam;
  bool _rotated;

public:
  CR_MOD(Eigen::VectorXd    THETA,
         const unsigned int OUTCOME,
         Eigen::VectorXd    EXT_COVARIATES,
         const unsigned int YB,
         const unsigned int YEAR_FIRST,
         const unsigned int YEAR_LAST,
         const unsigned int YEAR_LAST_EXAM,
         const bool         ROTATED):
  _theta(THETA),
  _outcome(OUTCOME),
  _ext_covariates(EXT_COVARIATES),
  _yb(YB),
  _year_first(YEAR_FIRST),
  _year_last(YEAR_LAST),
  _year_last_exam(YEAR_LAST_EXAM),
  _rotated(ROTATED){}

  // conditional log-likelihood
  double ll(const double ABILITY, const double SPEED);

  //gradient of conditional log-likelihood
  Eigen::VectorXd grll(const double ABILITY, const double SPEED);

};
double CR_MOD::ll(const double ABILITY, const double SPEED){

  Eigen::VectorXd covariates(_ext_covariates.size()+2);
  covariates << _ext_covariates, ABILITY, SPEED;
  double out = cr::outcomeLik(
    _outcome,
    _year_first,
    _year_last,
    _theta,
    covariates,
    _yb,
    _year_last_exam,
    true);

  return out;
}
Eigen::VectorXd CR_MOD::grll(const double ABILITY, const double SPEED){
  Eigen::VectorXd covariates(_ext_covariates.size()+2);
  covariates << _ext_covariates, ABILITY, SPEED;

  Eigen::VectorXd out=cr::grad::grl_outcomeLik(
    _outcome,
    _year_first,
    _year_last,
    _theta,
    covariates,
    _yb,
    _year_last_exam,
    _rotated
  );

  return(out);
}

// GRADED RESPONSE TIME CENSORED MODEL
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
      out += exams::examLik(exam,
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
      out += exams::examLik(exam,
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
      gr_irt += exams::grad::grl_examLik(exam,
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


// joint COMPETING RISK GRADED RESPONSE TIME CENSORED MODEL

#endif
