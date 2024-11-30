#ifndef exams_H
#define exams_H
#include "grades.h"
#include "times.h"
#include "extractParams.h"

namespace exams{
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
  double examLik(
      const unsigned int EXAM,
      const unsigned int GRADE,
      const double DAY,
      const double MAX_DAY,
      const bool OBSFLAG,
      const Eigen::Ref<const Eigen::VectorXd> THETA_IRT,
      const unsigned int N_GRADES,
      const unsigned int N_EXAMS,
      const double ABILITY,
      const double SPEED,
      const bool LOGFLAG = false
  ){
    double out, logpExam, logpTime;

    if(OBSFLAG){
      logpExam = exams::pGrade(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true);
      logpTime = exams::pTimeExam(EXAM, DAY, THETA_IRT, N_GRADES,  N_EXAMS, SPEED, false, true);
      out = logpExam+logpTime;
    }else{
      logpExam = exams::pGreaterGrades(1, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true);
      logpTime = exams::pTimeExam(EXAM, MAX_DAY, THETA_IRT, N_GRADES,  N_EXAMS, SPEED, true, true);
      out = log1mexp(-logpExam-logpTime);
      out = std::max(-10000.0, out);  // avoid log(0)
    }

    if(LOGFLAG){
      return(out);
    }else{
      return(exp(out));
    }

  }
}

namespace exams::grad{
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
 //' @param ROTATED Have the latent scores been rotated using their covariance matrix?
 //'
 //' @returns It returns the probability of observing or not a specific
 //' grade on a given exam before a given day conditioned on ability and speed.
 //'
 Eigen::VectorXd grl_examLik(
     const unsigned int EXAM,
     const unsigned int GRADE,
     const double DAY,
     const double MAX_DAY,
     const bool OBSFLAG,
     const Eigen::Ref<const Eigen::VectorXd> THETA_IRT,
     const unsigned int N_GRADES,
     const unsigned int N_EXAMS,
     const double ABILITY,
     const double SPEED,
     const bool ROTATED
 ){
   double logp, logpGrade, logpTime;
   const unsigned int dim_irt = N_EXAMS*(N_GRADES+3);

   Eigen::VectorXd gr = Eigen::VectorXd::Zero(dim_irt+2);


   if(OBSFLAG){
     logpGrade = exams::pGrade(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true);
     logpTime = exams::pTimeExam(EXAM, DAY, THETA_IRT, N_GRADES,  N_EXAMS, SPEED, false, true);
     Eigen::VectorXd grlTime = exams::grad::gr_pTimeExam(EXAM, DAY, THETA_IRT, N_GRADES, N_EXAMS, SPEED, ABILITY, false, ROTATED, true);
     Eigen::VectorXd grGrade = exams::grad::gr_pGrade(GRADE,EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY);
     gr = grlTime + grGrade /exp(logpGrade);
     logp = logpGrade+logpTime;
   }else{
     logpGrade = exams::pGreaterGrades(1, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true);
     logpTime = exams::pTimeExam(EXAM, MAX_DAY, THETA_IRT, N_GRADES,  N_EXAMS, SPEED, true, true);
     Eigen::VectorXd grGrade = exams::grad::gr_pGreaterGrades(1,EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY);
     Eigen::VectorXd grTime = exams::grad::gr_pTimeExam(EXAM, MAX_DAY, THETA_IRT, N_GRADES, N_EXAMS, SPEED, ABILITY, true, ROTATED, false);
     gr = -(grTime*exp(logpGrade) + exp(logpTime)*grGrade);
     logp = log1mexp(-logpGrade-logpTime);
     logp = std::max(-1000.0, logp);  // avoid log(0)
     gr/= exp(logp);
   }

   // gr/=std::max(1e-20, exp(logp));
   // logp/= 10;

   // Eigen::VectorXi is_selected = (gr.array() != 0).cast<int>();
   // gr/= exp(logp);


   return(gr);

 }
}


#endif
