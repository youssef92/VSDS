#ifndef PASSIVECONTROL_H
#define PASSIVECONTROL_H


#include "Utility.h"
#include "MotionGeneration.h"


class ClosedLoopControl {
protected:
    Mat D_eig_;
    int n_DOF_;
    realtype K_o_;
    realtype D_o_;
    realtype K_1_;
    realtype D_1_;

public:
    ClosedLoopControl( int M);
    Vec GetControlTorque(Vec x_dot_d, Vec x_dot, Mat Jac, realtype x_c, Vec x);
    Mat FindDampingBasis(Vec x);
    Vec ComputeOrientationForce(Mat  R_a, Mat R_d,Mat Jac,Vec q_dot);

};

class PassiveClosedLoopControl:public ClosedLoopControl{

private:
    realtype wind_act_ ;
    Vec x_att_ ;
    realtype s_tank_init_;
    realtype s_tank_max_ ;
    realtype s_tank_ ;
    realtype beta_ ;
    realtype ds_ ;
    realtype gamma_ ;
    Vec  stepness_ ;
    Vec  ko_;
    realtype Pow_VSDS_ ;


public:
    PassiveClosedLoopControl( int M,Vec x_0, Vec x_att);
    realtype Smooth_Activation_Exp(Vec x) ;
    realtype Energy_Tank_smooth(realtype Pow_VSDS,realtype Pow_Damp, Vec x, realtype dt) ;


    Vec GetControlForce_Passive(Vec x_dot_d, Vec x_dot, realtype x_c, Vec x,realtype dt) ;
    realtype Get_Tank_gamma() ;
    realtype Get_Tank_beta() ;
    realtype Get_Tank_Energy() ;


    realtype smooth_fall(realtype val,realtype hi, realtype lo) ;
    realtype smooth_rise(realtype val, realtype lo,realtype hi) ;
    realtype smooth_rise_2d(realtype x, realtype y, realtype xlo, realtype dx, realtype ylo, realtype dy) ;
    realtype smooth_rise_fall_2d(realtype x, realtype y, realtype xlo, realtype xhi,realtype dx, realtype ylo, realtype yhi, realtype dy) ;
    realtype smooth_rise_fall(realtype val, realtype l0,realtype l1,realtype h1, realtype h0) ;
    realtype smooth_fall_gamma(realtype val,realtype hi, realtype lo,realtype beta ) ;

    Vec Stiff_Potential_Consv(Vec x_err) ;
    realtype Get_Pow_VSDS_aka_z() ;

} ;

#endif

