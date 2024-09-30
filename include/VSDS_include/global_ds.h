#ifndef GLOBALDS_H_
#define GLOBALDS_H_

#include "utility.h"
#include "cmath"
#include <iostream>
#include "cmath"
#include "eigen3/Eigen/Dense"
using namespace special_math_functions;
using namespace GeneralFunction;



class GlobalDS{

private:

int n_DOF_ ;
Mat Mu_ ;
Vec Priors_ ;
Mat A_g_;
Mat Sigma_;
Vec att_;


public:

    GlobalDS(int n_dof) ;
    GlobalDS() ;

    Vec GlobalDS_eval(Vec x) ;
    realtype gaussPDF(Vec Data_in, Vec Mu , Mat Sigma) ;
    void Read_SEDS_FROMFile(string MuFile, string SigmaFile, string PriorsFile, string A_gFile,string att_File,Mat &Mu, Mat &Sigma, Vec &Priors, Mat &A_g, int M, int K_SEDS,Vec &att) ;
    Vec GetAttractor()  ;
    Mat FindDampingBasis(Vec x) ;
} ;






#endif
