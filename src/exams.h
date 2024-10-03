#ifndef exams_H
#define exams_H
#include "grades.h"
#include "times.h"
#include "extractParams.h"

//' Evaluate exam specific likelihood
//'
//' @param EXAM Exam of interest. Integers from 0 to N_EXAMS -1.
//' @param GRADE Grade used as reference. Integers from 1 to N_GRADES.
//' @param DAY Day of interest.
//' @param MAX_DAY Last day of observation.
//' @param OBSFLAG TRUE for observed, FALSE for not-observed.
//' @param THETA_IRT Portion of the parameter vector related to the IRT model
//' @param N_GRADES Number of grades modelled.
//' @param N_EXAMS Number of exams.
//' @param ABILITY ability value.
//' @param SPEED speed value.
//' @param LOGFLAG Set TRUE to return log value.
//'
//' @returns It returns the probability of observing or not a specific
//' grade on a given exam before a given day conditioned on ability and speed.
//'
//' @export
// [[Rcpp::export]]
double examLik(
    const unsigned int EXAM,
    const unsigned int GRADE,
    const double DAY,
    const double MAX_DAY,
    const bool OBSFLAG,
    Eigen::VectorXd& THETA_IRT,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const double ABILITY,
    const double SPEED,
    const bool LOGFLAG = false
){
  double out, logpExam, logpTime;

  if(OBSFLAG){
    logpExam = pGrade(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true);
    logpTime = pTimeExam(EXAM, DAY, THETA_IRT, N_GRADES,  N_EXAMS, SPEED, false, true);
    out = logpExam+logpTime;
  }else{
    logpExam = pGreaterGrades(1, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true);
    logpTime = pTimeExam(EXAM, MAX_DAY, THETA_IRT, N_GRADES,  N_EXAMS, SPEED, true, true);
    out = log1mexp(-logpExam-logpTime);
    out = std::max(-10000.0, out);  // avoid log(0)
  }

  if(LOGFLAG){
    return(out);
  }else{
    return(exp(out));
  }


  return(out);

}


//' Evaluate the gradient of exam-specific log-likelihood
//'
//' @param EXAM Exam of interest. Integers from 0 to N_EXAMS -1.
//' @param GRADE Grade used as reference. Integers from 1 to N_GRADES.
//' @param DAY Day of interest.
//' @param MAX_DAY Last day of observation.
//' @param OBSFLAG TRUE for observed, FALSE for not-observed.
//' @param THETA_IRT Portion of the parameter vector related to the IRT model
//' @param N_GRADES Number of grades modelled.
//' @param N_EXAMS Number of exams.
//' @param ABILITY ability value.
//' @param SPEED speed value.
//'
//' @returns It returns the probability of observing or not a specific
//' grade on a given exam before a given day conditioned on ability and speed.
//'
//' @export
// [[Rcpp::export]]
Eigen::VectorXd grl_examLik(
    const unsigned int EXAM,
    const unsigned int GRADE,
    const double DAY,
    const double MAX_DAY,
    const bool OBSFLAG,
    Eigen::VectorXd& THETA_IRT,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const double ABILITY,
    const double SPEED,
    const bool ROTATED
){
  double logp, logpExam, logpTime;
  Eigen::VectorXd gr = Eigen::VectorXd::Zero(THETA_IRT.size());


  if(OBSFLAG){
     logpExam = pGrade(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true);
     logpTime = pTimeExam(EXAM, DAY, THETA_IRT, N_GRADES,  N_EXAMS, SPEED, false, true);
     Eigen::VectorXd grTime = gr_pTimeExam(EXAM, DAY, THETA_IRT, N_GRADES, N_EXAMS, SPEED, ABILITY, false, ROTATED);
     Eigen::VectorXd grExam = gr_pGrade(GRADE,EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY);

     gr = grTime*exp(logpExam) + exp(logpTime)*grExam;
     logp = logpExam+logpTime;
  }else{
     logpExam = pGreaterGrades(1, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true);
     logpTime = pTimeExam(EXAM, MAX_DAY, THETA_IRT, N_GRADES,  N_EXAMS, SPEED, true, true);
     Eigen::VectorXd grExam = gr_pGreaterGrades(1,EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY);
     Eigen::VectorXd grTime = gr_pTimeExam(EXAM, MAX_DAY, THETA_IRT, N_GRADES, N_EXAMS, SPEED, ABILITY, true, ROTATED);

     gr = -(grTime*exp(logpExam) + exp(logpTime)*grExam);
     logp = log1mexp(-logpExam-logpTime);
     logp = std::max(-10000.0, logp);  // avoid log(0)
   }

   gr/=std::max(1e-16, exp(logp));


   return(gr);

 }
#endif
