/*
 * Author: Youssef Michel
 * Description: This scipt implements the damping controller from https://ieeexplore.ieee.org/document/7358081.
 *              It overrides the base VSDS function of setting the desired force
 *              For more details, please check https://ieeexplore.ieee.org/abstract/document/10288346
 *
 *
 */


#include "vsds_base.h"
#include "cmath"
#include "eigen3/Eigen/Dense"

using namespace special_math_functions;
using namespace GeneralFunction;
using Eigen::MatrixXd;
using namespace Eigen;

#include <iostream>

namespace vsds_transl_control {

Vec VSDS_picds::GetDesiredForce(Vec x,Vec x_dot,Vec q)
{
    Vec x_t(n_DOF_);
    Vec x_dot_t(n_DOF_) ;

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

    Vec g = ComputeOmega(x_t);
    Mat Q = GlobalDS_.FindDampingBasis(GlobalDS_ .GlobalDS_eval(x_t));
    Vec F_control(2*n_DOF_) ;

    std::vector<double> d_picds;
    Mat L=Mat::Identity(2,2)*300 ;
    if (ros::param::get("/D_PICDS", d_picds)) {
        Eigen::VectorXd eigenVector = Eigen::Map<Eigen::VectorXd>(d_picds.data(), d_picds.size());
        L = eigenVector.asDiagonal();
    }
    else {
        ROS_WARN("D_PICDS Param Not Found !!") ;
    }


    Mat D_PICDS=Q * L * Q.transpose();
    Vec x_dot_d=pre_scale_*GlobalDS_ .GlobalDS_eval(x_t) ;
    Vec F_PICDS= D_PICDS*(x_dot_d - x_dot_t) ;
    F_control << F_PICDS , Vec::Zero(n_DOF_) ;


    return F_control;
}
}
