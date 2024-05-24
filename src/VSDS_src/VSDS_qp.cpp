
/*
 * Author: Youssef Michel
 * Description: This scipt implements the VSDS approach based on QP-optimization. It overrides the base VSDS function of
 *              setting the desired force
 *              For more details, please check https://ieeexplore.ieee.org/abstract/document/10288346
 */


#include "VSDS_base.h"
#include "cmath"
#include "eigen3/Eigen/Dense"

using namespace special_math_functions;
using namespace GeneralFunction;
using Eigen::MatrixXd;
using namespace Eigen;

#include <iostream>

Vec VSDS_qp::GetDesiredForce(Vec x,Vec x_dot,Vec q)
{
    Vec x_t(n_DOF_);
    Vec x_dot_t(n_DOF_) ;
    ROS_INFO("QP TEST ") ;


    if (x.size() == 3 && n_DOF_ == 2){
        x_t(0) = x(1);
        x_t(1) = x(2);
        x_dot_t(0)=x_dot(1) ;
        x_dot_t(1)=x_dot(2) ;

    }
    else{
        x_t = x;
        x_dot_t=x_dot ;

    };

    Mat fl = Mat::Zero(n_DOF_,K_);
    Mat fd = Mat::Zero(n_DOF_,K_);

    Vec g = Omega(x_t);
    realtype act = StartActivation(x_t);
    Mat Q = GlobalDS_.FindDampingBasis(GlobalDS_ .global_ds_eval(x_t));
    Vec F_control(2*n_DOF_) ;

        for (int i=0; i<K_; i++)
        {

                Mat D_k(2,2) ;
                Mat L= Mat::Identity(2,2) ;
                L(0,0)=0.8 ;
                L(1,1)=0.8;
                D_k=-A_.block(0,n_DOF_*i,n_DOF_,n_DOF_)*6 ;
                D_k= 2*L*D_k.pow(0.5) ;
                D_k = Q * D_k * Q.transpose();
                fl.col(i) =(g(i) * A_.block(0,n_DOF_*i,n_DOF_,n_DOF_) * (x_t - x_rec_.block(0,i+1,n_DOF_,1)) ) + g(i)*Force_fields_.block(0,i,n_DOF_,1)   ; // Add local force fields to improve tracking
                fd.col(i)=- g(i)*D_k*(x_dot_t) ;

        }

        Vec fdmp_sum= fd.rowwise().sum();
        Vec fvs_sum =  act *fl.rowwise().sum();


        F_control<<fvs_sum,fdmp_sum ;




    return F_control;
}
