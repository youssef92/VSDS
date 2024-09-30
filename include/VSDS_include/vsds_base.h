#ifndef VSDS_BASE_H_
#define VSDS_BASE_H_

#include "utility.h"
#include "global_ds.h"

/*
 * Author: Youssef Michel
 * Description: Base class implementation for VSDS, which is adapted depending
 *              on the chosen VSDS formulation.
 *              For more details, please check https://ieeexplore.ieee.org/abstract/document/10288346
 *
 * Notes:
 *  - Note that VSDS implementations only differ in the setting the desire force function, hence
 *    this is the function to be overriden
 *  - We assume that the first order DS has been already learnt, and its relevant paramaters are provided in .txt files
 *  - You can either use one of the provided demos in the /config directory, or use your own DS, for which you have to provide the parameters
 */



namespace vsds_transl_control {

class VSDS_base {

protected:
    int n_viapoints_; // number of via points
    int n_DOF_; // number of data dimension
    int n_VSDS_; // number of VSDS models
    int N_init_;
    Mat x_cen_;
    Vec x_len_;
    realtype pre_scale_;
    Vec att_;
    Vec x0_;
    realtype th_begin_;
    realtype sigmascale_;
    realtype dt_;
    realtype threshold_tube_;
    realtype traj_len_;
    GlobalDS GlobalDS_ ;
    Mat Force_fields_ ;
    Mat Damping_Fields_ ;
    ros::NodeHandle node_handle_;
    void Read_ForceFields_FromFile(string file_forcefields, Mat &Force_fields, int M, int N_viapoints) ;
    Mat x_rec_;
    Mat A_;
    Mat B_;
    Vec ComputeOmega(Vec x);
    Mat GetTempPoints(Vec x0, realtype dt, realtype threshold, bool first_gen);
    void Get_A_matrix(); // System Matrix (Matrix A in the paper)


public:


    VSDS_base(Vec x_0) ;
    virtual Vec GetDesiredForce(Vec x,Vec x_dot,Vec q);
    realtype StartActivation(Vec x);
    Mat GetStiffness(Vec x);
    void GetViaPointsReduced(Mat& x_temp, Vec x0);
    realtype PosCheck(Vec x);
    void MotionRegenerate(Vec x);
    Mat GetX_rec() ;
    VSDS_base() ;
    Vec GetAttractor() ;


};

class VSDS_org: public VSDS_base{

public:

    using VSDS_base::VSDS_base ;
    Vec GetDesiredForce(Vec x,Vec x_dot,Vec q) override;

} ;

class VSDS_qp: public VSDS_base{

public:

    using VSDS_base::VSDS_base ;
    Vec GetDesiredForce(Vec x,Vec x_dot,Vec q) override;

} ;

class VSDS_vf: public VSDS_base{

public:
    using VSDS_base::VSDS_base ;
    Vec GetDesiredForce(Vec x,Vec x_dot,Vec q) override;

} ;

class VSDS_picds: public VSDS_base{

public:
    using VSDS_base::VSDS_base ;
    Vec GetDesiredForce(Vec x,Vec x_dot,Vec q) override;

} ;

}


#endif
