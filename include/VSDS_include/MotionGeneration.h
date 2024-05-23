#ifndef MOTIONGENERATION_H_
#define MOTIONGENERATION_H_

#include "Utility.h"

#include "global_ds.h"

//enum Tracking_Terms{
//    VarDamping_Flag=0 ,
//    ForceField_Flag=1 ,
//} ;



class VSDS {

private:
    int K_; // number of local DS
    int n_DOF_; // number of data dimension
    int N_; // number of data points
    int   VarDamping_Flag_ ;
    int   ForceField_Flag_ ;
    int   PICDS_Flag_ ;

    int N_init_;
    Mat x_cen_;
    Vec x_len_;
    realtype pre_scale_;
    Vec att_;
    Vec x0_;
    realtype th_begin_;
    realtype sigmascale_;
    realtype dt_;
    realtype threshold_;
    realtype traj_len_;
    realtype F_exp_;
    global_ds GlobalDS_ ;
    Mat Force_fields_ ;
    Mat Damping_Fields_ ;
    ros::NodeHandle nh;
    void Read_ForceFields_FromFile(string file_forcefields, Mat &Force_fields, int M, int N_viapoints) ;


public:


    VSDS(Vec x_0) ;
    Mat x_rec_;
    Mat A_;
    Mat B_;
    Vec Omega(Vec x);
    Mat GetTempPoints(Vec x0, realtype dt, realtype threshold, bool first_gen);
    void GetABMatrix();
    Vec GetDesiredForce(Vec x,Vec x_dot,Vec q);
    realtype StartActivation(Vec x);
    Mat GetStiffness(Vec x);
    void GetViaPointsReduced(Mat x_temp, Vec x0);
    realtype PosCheck(Vec x);
    void MotionRegenerate(Vec x);
    Mat GetX_rec() ;
    VSDS() ;
    Vec GetAttractor() ;
    

};



#endif
