#include "PassiveControl.h"
#include "iostream"
#include "eigen3/Eigen/Dense"
#include "cmath"


using Eigen::MatrixXd;
using namespace Eigen;

using namespace std;
using namespace special_math_functions;
using namespace StiffnessProfiles;

// Implements the passification of the VSDS control law (energy tank +conservative potential),
// and also the damping injection
// For more details, please check https://ieeexplore.ieee.org/abstract/document/10288346

ClosedLoopControl::ClosedLoopControl( int M):
    n_DOF_(M)

{
    K_o_ = 40.0;
    D_o_ = 15;
    K_1_ = 2000;
    D_1_ = 60;
    D_eig_=Mat::Zero(2,2) ;
    D_eig_ << 64, 0,
            0, 72;

}

Vec ClosedLoopControl::GetControlTorque(Vec x_dot_d, Vec x_dot, Mat Jac, realtype x_c, Vec x)
{
    Vec F(3);
    Mat Q = FindDampingBasis(x_dot_d);
    Mat D = Q * D_eig_ * Q.transpose();
    if (n_DOF_ == 3){
        F = - D * x_dot + x_dot_d;
    }
    else{
        F(0) = K_1_ * (x_c - x(0)) - D_1_ * x_dot(0);
        Vec F_temp = - D * x_dot.block(1,0,2,1) + D_eig_(0,0) * x_dot_d;
        F(1) = F_temp(0);
        F(2) = F_temp(1);
    }


    Mat Jac_t = Jac.block(0, 0, 3, 7);
    Vec tau_x = Jac_t.transpose() * F;
    return F ;
}

Vec ClosedLoopControl::ComputeOrientationForce(Mat  R_a, Mat R_d,Mat Jac,Vec q_dot) {


    // Implementation Based on Franka Emika
    Eigen::Matrix3d R_aa=R_a;
    Eigen::Matrix3d R_dd = R_d;
    Vec error=Vec::Zero(3) ;

    Eigen::Quaterniond orientation_d(R_dd);
    Eigen::Quaterniond orientation(R_aa);
    
    if (orientation_d.coeffs().dot(orientation.coeffs()) < 0.0) {
        orientation.coeffs() << -orientation.coeffs();
    }

    Eigen::Quaterniond error_quaternion(orientation.inverse() * orientation_d);
    error << error_quaternion.x(), error_quaternion.y(), error_quaternion.z();
    error<< R_a * error;

    Eigen::Matrix3d R_ad=R_aa.transpose()*R_dd ;
    Eigen::Quaterniond orientation_ad(R_ad);
    error << orientation_ad.x(), orientation_ad.y(), orientation_ad.z();
    error<< R_a * error;


    Mat Jac_o = Jac.block(3, 0, 3, 7);

    Vec F_o = K_o_*error -D_o_*( (Jac_o*q_dot))   ;

    Vec tau_o=(Jac_o.transpose() * F_o) ;

    return  tau_o;

}

Mat ClosedLoopControl::FindDampingBasis(Vec x)
{
    if (n_DOF_ == 2){
        Vec y(2);
        y << 1,
                -x(0)/(x(1)+eps());
        Mat B = Mat::Zero(2,2);
        B.col(0) = x/x.norm();
        B.col(1) = y/y.norm();
        return B;
    } else {
        Vector3d x_ = x;
        Vector3d y(1,1,(-x(0)-x(1))/(x(2)+eps()));
        Vector3d z;
        z = x_.cross(y);
        Mat B = Mat::Zero(3,3);
        B.col(0) = x/x.norm();
        B.col(1) = y/y.norm();
        B.col(2) = z/z.norm();
        return B;
    }
}

//----------------------------------Energy Tank and Stablization Implementation----------------------
PassiveClosedLoopControl::PassiveClosedLoopControl(int M,Vec x_0,Vec x_att_):
    ClosedLoopControl(M){


    D_eig_ << 64, 0,
            0, 72;
    

    wind_act_=4800 ; //4800
    x_att_=x_att_ ;
    realtype act=Smooth_Activation_Exp(x_0-x_att_) ;

    stepness_ = x_att_ ; //0.0006
    stepness_ <<0.0006,0.004 ;


    ko_=0*x_att_ ;
    ko_<< 0.50, 0.55;

    beta_=1 ;
    gamma_=beta_ ;
    s_tank_init_=30 ;
    s_tank_max_=act*s_tank_init_ ;
    s_tank_=0.79*s_tank_max_ ;
    ds_=0.1*s_tank_max_ ;
    Pow_VSDS_=0 ;

}


realtype PassiveClosedLoopControl::Smooth_Activation_Exp(Vec x){

    return (1-exp(-wind_act_*x.dot(x))) ;
}

realtype PassiveClosedLoopControl::Get_Pow_VSDS_aka_z(){
    return Pow_VSDS_ ;
}




realtype PassiveClosedLoopControl::Get_Tank_gamma(){
    return gamma_ ;
}

realtype PassiveClosedLoopControl::Get_Tank_beta(){
    return beta_ ;
}

realtype PassiveClosedLoopControl::Get_Tank_Energy(){
    return s_tank_ ;
}


Vec PassiveClosedLoopControl::Stiff_Potential_Consv(Vec x_err) {

    int size=max(x_err.rows(),x_err.cols()) ;
    Vec center = Vec::Zero(size) ;
    Vec F = Vec::Zero(size) ;

    realtype tau_min = 0.0;

    for (int i=0 ; i< size ; i++){

        F(i)=-(ko_(i)/stepness_(i))*x_err(i)*exp(-pow(x_err(i)-center(i),2) / (2*stepness_(i))) - (2*tau_min*x_err(i)) ;

    }

    return F ;
}


Vec PassiveClosedLoopControl::GetControlForce_Passive(Vec x_dot_d, Vec x_dot, realtype x_c, Vec x,realtype dt)
{
    Vec F=Vec::Zero(3);
    Mat Q = FindDampingBasis(x_dot_d.head(2));
    Mat D = Q * D_eig_ * Q.transpose();
    string VS_type="qp" ; // org , qp , fd
    ros::param::get("/VS_type",VS_type) ;

    if (n_DOF_ == 3){

        Vec F_VSDS= x_dot_d ;
        Vec F_Damp= - D * x_dot ;

        realtype Pow_Damp= x_dot.transpose()*F_Damp  ;
        realtype Pow_VSDS= x_dot.transpose() *F_VSDS;

        Vec F_constPot=Stiff_Potential_Consv(x.tail(2)-x_att_) ;
        realtype beta=Energy_Tank_smooth(Pow_VSDS,Pow_Damp, x,  dt) ;
        Vec F_tmp= beta*F_VSDS +F_Damp +F_constPot ;
        F=F_tmp ;

    }
    else{

        F(0) = K_1_ * (x_c - x(0)) - D_1_ * x_dot(0);

        if(VS_type=="PI") {
            Vec F_tmp= x_dot_d.head(2) ;
            F(1)=F_tmp(0) ; F(2) =F_tmp(1) ;

        }
        else {

            Vec F_constPot=Stiff_Potential_Consv(x.tail(2)-x_att_) ;
            Vec F_VSDS=  1.0*Smooth_Activation_Exp( x_att_-x.tail(2))*x_dot_d.head(2) ;
            Vec F_Damp= x_dot_d.tail(2) ;

            realtype Pow_Damp= -x_dot.tail(2).transpose()  *   F_Damp  ; // Added - since should be positive power
            Pow_VSDS_=  x_dot.tail(2).transpose()   *  F_VSDS ;

            Energy_Tank_smooth(Pow_VSDS_,Pow_Damp, x.tail(2),  dt) ; // gamma will get updated

            Vec  F_temp(2) ;
            F_temp= gamma_*F_VSDS +F_Damp + 1*F_constPot ;

            F(1) = F_temp(0);
            F(2) = F_temp(1);
        }

    }

    return F ;
}

realtype PassiveClosedLoopControl::Energy_Tank_smooth(realtype Pow_VSDS,realtype Pow_Damp, Vec x, realtype dt){

    realtype act=Smooth_Activation_Exp(x-x_att_) ;
    s_tank_max_= act*s_tank_init_ ;
    ds_=0.1*s_tank_max_ ;
    ds_=2 ;
    realtype s_tank_min =0 ;
    realtype z=Pow_VSDS ;
    realtype dz=0.1 ;


    if(abs(z)<0.001){
        z=0 ;
    }

    beta_  = smooth_rise_fall_2d(z,s_tank_, 0.000, 0.000, dz, s_tank_min, s_tank_max_, ds_);
    beta_=max(beta_,0.0) ;
    gamma_ = smooth_fall_gamma(z,-0.1, 0, beta_ ) ;
    realtype alpha = smooth_rise_fall(s_tank_, 0, ds_, s_tank_max_-ds_, s_tank_max_);

    //New Formulation
    alpha = min(0.99, alpha);
    realtype sdot = alpha*Pow_Damp - beta_*z  -(1.05-act)*s_tank_ ;
    s_tank_ = s_tank_ + sdot*dt ;

    return gamma_ ;

}


realtype PassiveClosedLoopControl::smooth_fall(realtype val,realtype hi, realtype lo){

    realtype smooth_value=1 ;

    if(hi>lo) {
        cout<<"Error in Smooth Fall"<<endl ;
        smooth_value = -1;
    }

    if(val>=lo){
        smooth_value = 0; }
    else if(val<hi){
        smooth_value = 1;}
    else{
        realtype T = 2*(lo-hi);
        smooth_value = 0.5+0.5*sin(2*pi*(val-lo)/T - pi*0.5);
    }
    return smooth_value ;

}

realtype PassiveClosedLoopControl::smooth_rise(realtype val, realtype lo,realtype hi){
    realtype smooth_value=1 ;
    if(lo>hi){
        cout<<"Error in Smooth Rise"<<endl ;
        smooth_value = -1;

    }

    if(val>=hi){
        smooth_value = 1; }
    else if(val<lo){
        smooth_value = 0;}
    else{
        realtype T = 2*(hi-lo);
        smooth_value = 0.5+0.5*sin(2*pi*(val-lo)/T - pi*0.5);
    }
    return smooth_value ;
}

realtype PassiveClosedLoopControl::smooth_rise_2d(realtype x, realtype y, realtype xlo, realtype dx, realtype ylo, realtype dy){
    realtype h1 = smooth_rise(x, xlo, xlo+dx);
    realtype h2 = smooth_fall(y, ylo, ylo+dy);

    realtype h = 1 - h1*h2;
    return h ;
}

realtype PassiveClosedLoopControl::smooth_rise_fall_2d( realtype x, realtype y, realtype xlo, realtype xhi,realtype dx, realtype ylo, realtype yhi, realtype dy){
    realtype h1_1 = smooth_rise(x, xlo-dx, xlo);
    realtype h1_2 = smooth_fall(y, ylo, ylo+dy);
    realtype h1 = h1_1*h1_2;

    realtype h2_1 = smooth_fall(x, xhi, xhi+dx);
    realtype h2_2 = smooth_rise(y, yhi-dy, yhi);
    realtype h2 = h2_1*h2_2;

    realtype h = 1 - h1 - h2;
    return h ;
}

realtype PassiveClosedLoopControl::smooth_rise_fall(realtype val, realtype l0,realtype l1,realtype h1, realtype h0){
    realtype hr = smooth_rise( val,  l0,  l1);
    realtype hf = smooth_fall( val, h1, h0);
    realtype h = hr*hf;
    return h ;
}


realtype PassiveClosedLoopControl::smooth_fall_gamma(realtype val,realtype hi, realtype lo,realtype beta ){

    realtype smooth_value=1 ;

    if(hi>lo) {
        cout<<"Error in Smooth Fall"<<endl ;
        smooth_value = -1;
    }
    realtype scale=(1-beta) ;

    if(val>=lo){
        smooth_value = beta; }
    else if(val<hi){
        smooth_value = 1;}
    else{
        realtype T = 2*(lo-hi);
        smooth_value = (0.5+0.5*sin(2*pi*(val-lo)/T - pi*0.5))*scale + beta;
    }
    return smooth_value ;

}


