#include "MotionGeneration.h"
#include "cmath"
#include "eigen3/Eigen/Dense"

using namespace special_math_functions;
using namespace GeneralFunction;
using Eigen::MatrixXd;
using namespace Eigen;

#include <iostream>


void VSDS:: Read_ForceFields_FromFile(string file_forcefields, Mat &Force_fields, int M, int N_viapoints){

    std::ifstream ForceFileRead;

    ForceFileRead.open(file_forcefields.c_str());
    if (ForceFileRead){
        for(int i = 0; i < M; i++){
            for(int j=0; j<N_viapoints; j++){
                ForceFileRead >> Force_fields(i,j);
            }
        }
        ForceFileRead.close();
    }
    else {
        std::cout << "no MU file" << std::endl;
    }


}

VSDS::VSDS(Vec x0):
    x0_(x0){

    ROS_INFO("--------------Initializing VSDS-----------") ;
    int K_SEDS=3 ;
    nh.getParam("/K_SEDS", K_SEDS);
    realtype prescale = 2.0; // Scales the Global DS during the via points generation 2.0 for force fields
    nh.getParam("/pre_scale",prescale) ;


    N_= 20 ;
    x_rec_= x0 ;
    if ( ! ros::param::get("/N_VSDS",N_) ){
    ROS_WARN("Number of Via Points !!");
    }

    N_init_= N_ ;
    n_DOF_=2 ;
    sigmascale_=1.0 ;

    std::string DS_ModelName ;
    nh.getParam("/VSDS_name",DS_ModelName);
    std::string packPath = ros::package::getPath("VSDS");

    string file_forcefields= packPath + "/config/"+ DS_ModelName+ "/force_fields.txt";

    GlobalDS_  = global_ds() ;

    Force_fields_=Mat::Ones(n_DOF_,N_) ;
    att_=Vec::Zero(n_DOF_) ;
    att_= GlobalDS_.GetAttractor() ;
    Read_ForceFields_FromFile( file_forcefields,   Force_fields_, n_DOF_, N_) ;

    realtype  damping_vf=230  ;
    if ( ! nh.getParam("/damping_vf",damping_vf))
        ROS_WARN("Damping Gain for VF condition Not Found !! ");

    Mat Damping_fields=Mat::Ones(n_DOF_,N_) ;
    Damping_fields.row(0)=Damping_fields.row(0)*damping_vf ;
    Damping_fields.row(1)=Damping_fields.row(1)*damping_vf;

    threshold_ = 0.0001 ;
    dt_ = 0.01;


    pre_scale_ = 2.0; // Scales the Global DS during the via points generation 2.0 for force fields
    if ( ! nh.getParam("/pre_scale",pre_scale_))
        ROS_WARN("PreScale param Not found !!");

    Mat xtemp;
    cout <<"x0: "<<x0 <<endl ;
    xtemp = GetTempPoints(x0, dt_, threshold_, true);
    GetViaPointsReduced(xtemp,x0);
    GetABMatrix();

    string VS_type="qp" ; // org , qp , fd
    nh.getParam("/VS_type",VS_type) ;


    if(VS_type=="qp"){
        ForceField_Flag_=1 ;
    }
    else if(VS_type=="fd"){
        VarDamping_Flag_=1 ;

    }
    else if(VS_type=="PI"){
        VarDamping_Flag_=0 ;
        ForceField_Flag_=0 ;
        PICDS_Flag_=1 ;
    }

    else{
        VarDamping_Flag_=0 ;
        ForceField_Flag_=0 ;
        PICDS_Flag_=0 ;
    }

   ROS_INFO("Initialized VSDS !!! ");
}


Vec VSDS::Omega(Vec x)
{
    Vec omega = Vec::Zero(K_);
    Vec delta = sigmascale_ * x_len_;
    Vec x_t(n_DOF_);
    if (x.size() == 3 && n_DOF_ == 2){
        x_t(0) = x(1);
        x_t(1) = x(2);
    }
    else{
        x_t = x;
    };

    for(int i=0;i <K_;i++){
        omega(i) = exp(-(1/(2*delta(i)*delta(i)))*(x_t-x_cen_.block(0,i,n_DOF_,1)).transpose()*(x_t-x_cen_.block(0,i,n_DOF_,1)));
    }

    realtype omega_sum = omega.sum();
    omega = omega / omega_sum;
    return omega;
}

Mat VSDS::GetTempPoints(Vec x0, realtype dt, realtype threshold, bool first_gen)
{
    Vec xnow = x0;
    Vec xd = Vec::Zero(n_DOF_);
    Mat x_temp;
    realtype l_sum;

    while ((xnow - att_).norm() > threshold) {
        xd = pre_scale_* GlobalDS_.global_ds_eval(xnow);

        x_temp = PushBackColumn(x_temp, xnow);
        xnow = xnow + dt*xd;
    }
    x_temp = PushBackColumn(x_temp, att_);
    l_sum = (x_temp.block(0,0,n_DOF_,x_temp.cols()-1) - x_temp.block(0,1,n_DOF_,x_temp.cols()-1)).colwise().norm().sum();
    std::cout << "l_sum:" << l_sum << std::endl;
    th_begin_ = 0.1 * l_sum;

    if (first_gen) {
        traj_len_ = l_sum;
    }
    else{
        N_ = (int)ceil(N_init_ * (l_sum / traj_len_));
    }

    return x_temp;
}

void VSDS::GetViaPointsReduced(Mat x_temp, Vec x0)
{
    int temp_size = x_temp.cols();
    realtype dis = 0;
    int i = 1;
    int j = 0;
    Vec xnow = x0;
    Vec len(N_);
    Vec x_len_temp = (x_temp.block(0,0,n_DOF_,x_temp.cols()-1) - x_temp.block(0,1,n_DOF_,x_temp.cols()-1)).colwise().norm();
    realtype l_sum = x_len_temp.sum();
    realtype len_sub = l_sum/(N_);
    for (int j=0; j<N_-1; j++){
        len(j) = (j+1)*len_sub;
    }
    len(N_-1) = l_sum;
    j = 0;
    while (i < temp_size) {
        if (dis < len(j)) {
            dis = dis + x_len_temp(i-1);
            i++;
        }
        else {
            xnow = x_temp.block(0,i-1,2,1);
            x_rec_ = PushBackColumn(x_rec_,xnow);
            dis = dis + x_len_temp(i-1);
            i++;
            j++;
        }
    }
    x_rec_ = PushBackColumn(x_rec_, att_);
    K_ = x_rec_.cols()-1;
    x_cen_ = (x_rec_.block(0,0,n_DOF_,x_rec_.cols()-1) + x_rec_.block(0,1,n_DOF_,x_rec_.cols()-1))/2;
    x_len_ = (x_rec_.block(0,0,n_DOF_,x_rec_.cols()-1) - x_rec_.block(0,1,n_DOF_,x_rec_.cols()-1)).colwise().norm();
    std::cout << "N:" << x_rec_.cols() << std::endl;
    std::cout << "X_rec:" << x_rec_ << std::endl;
    ROS_INFO("Computed Via Points") ;
}

Mat VSDS::GetX_rec(){
    return x_rec_ ;
}

void VSDS::GetABMatrix()
{
    Mat B_temp;
    Mat A_temp;
    for (int i=0; i<K_; i++)
    {
        B_temp = GlobalDS_.FindDampingBasis(x_rec_.block(0,i+1,n_DOF_,1) - x_rec_.block(0,i,n_DOF_,1));
        A_temp = -B_temp * GetStiffness(x_rec_.block(0,i,n_DOF_,1)) * B_temp.transpose();
        B_ = PushBackColumnMat(B_,B_temp);
        A_ = PushBackColumnMat(A_,A_temp);
    }
    return;

}



Vec VSDS::GetDesiredForce(Vec x,Vec x_dot,Vec q)
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

    Mat fl = Mat::Zero(n_DOF_,K_);
    Mat fd = Mat::Zero(n_DOF_,K_);

    Vec g = Omega(x_t);
    realtype act = StartActivation(x_t);
    Vec x_dot_d=pre_scale_*GlobalDS_ .global_ds_eval(x_t) ;

    Mat Q = GlobalDS_.FindDampingBasis(GlobalDS_ .global_ds_eval(x_t));

    Vec F_control(2*n_DOF_) ;

    if(PICDS_Flag_==0){

        for (int i=0; i<K_; i++)
        {

            if(VarDamping_Flag_ && ForceField_Flag_){

                fl.col(i) =(g(i) * A_.block(0,n_DOF_*i,n_DOF_,n_DOF_) * (x_t - x_rec_.block(0,i+1,n_DOF_,1)) ) ;
            }

            else if(VarDamping_Flag_){
                Mat D_k(2,2) ;
                D_k.diagonal()=Damping_Fields_.block(0,i,n_DOF_,1) ;
                fl.col(i) =(g(i) * A_.block(0,n_DOF_*i,n_DOF_,n_DOF_) * (x_t - x_rec_.block(0,i+1,n_DOF_,1)) ) + 1*g(i)*D_k*x_dot_d  ;
                fd.col(i)=- g(i)*D_k*(x_dot_t) ;

            }
            else if (ForceField_Flag_){
                Mat D_k(2,2) ;
                Mat L= Mat::Identity(2,2) ;
                L(0,0)=0.79 ; // 0.55
                L(1,1)=0.79; //0.37 for LASA DataSet
                D_k=-A_.block(0,n_DOF_*i,n_DOF_,n_DOF_)*6 ;
                D_k= 2*L*D_k.pow(0.5) ;
                D_k = Q * D_k * Q.transpose();

                fl.col(i) =(g(i) * A_.block(0,n_DOF_*i,n_DOF_,n_DOF_) * (x_t - x_rec_.block(0,i+1,n_DOF_,1)) ) + g(i)*Force_fields_.block(0,i,n_DOF_,1)   ; // Add local force fields to improve tracking
                fd.col(i)=- g(i)*D_k*(x_dot_t) ;

            }
            else{
                Mat D_k(2,2) ;
                D_k=-A_.block(0,n_DOF_*i,n_DOF_,n_DOF_)*4 ; //6
                D_k= 2*0.33*special_math_functions::Eig_Decomp_Sqrt(D_k) ; //0.47 or 0.59
                fl.col(i) =(g(i) * A_.block(0,n_DOF_*i,n_DOF_,n_DOF_) * (x_t - x_rec_.block(0,i+1,n_DOF_,1)) ) ;
                fd.col(i)=- g(i)*D_k*(x_dot_t) ;
            }


        }

        Vec fdmp_sum= fd.rowwise().sum();
        Vec fvs_sum =  act *fl.rowwise().sum();

        if( ForceField_Flag_==1 ){
            fvs_sum = 1 * fl.rowwise().sum();
        }

        F_control<<fvs_sum,fdmp_sum ;

    }
    else {

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
        Vec F_PICDS= D_PICDS*(x_dot_d - x_dot_t) ;
        F_control << F_PICDS , Vec::Zero(n_DOF_) ;

    }

    return F_control;
}

realtype VSDS::StartActivation(Vec x)
{
    //    realtype b = 0.01;
    realtype b = 0.6;
    realtype temp;
    Vec x_t(n_DOF_);
    if (x.size() == 3 && n_DOF_ == 2){
        x_t(0) = x(1);
        x_t(1) = x(2);
    }
    else{
        x_t = x;
    };

    if ((x_t-x0_).norm() > th_begin_)
    {
        return 1;
    } else {
        temp = asin(1-b) * ((x_t-x0_).norm()/th_begin_);
        return sin(temp) + b;
    }
}

Mat VSDS::GetStiffness(Vec x)
{
    Mat K_des(n_DOF_,n_DOF_);
    if (n_DOF_ == 2){

        realtype k11;
        realtype k22;

        realtype t_max = 0.05;
        realtype t_min = 0.0;
        realtype k22_min = 100;
        realtype k22_diff = 600;
        k22 = k22_min + smooth_transition_rising(x(1),t_max,t_min,k22_diff);

        realtype e_max = 0.05;//0.2;
        realtype e_min = 0.0;//0.17;
        realtype k11_min = 700;
        realtype k11_diff = 500;
        k11 = k11_min + k11_diff*smooth_transition_fall(x(1),e_max,e_min);

        k11=(800 + 500*(sin(10*x(0)+0.8)+1)/2) ;
        k22=(700 + 400*(sin(10*x(1)+0.8)+1)/2) ;

        k11=(800 +  300*(sin(15*x(0)+0.8)+1)/2) ; // Used For Lasa Data Set
        k22=(1000 + 400*(sin(15*x(0)+0.8)+1)/2) ; // Used For Lasa Data Set

   //    k11= 1000 + 800*smooth_transition_fall(x(1),0.16,0.04); // Used For Drilling
  //    k22= 1000 + 800*smooth_transition_fall(x(1),0.16,0.04); // Used For Drilling

        string Stiff_type ;
        nh.getParam("/Stiff_type",Stiff_type) ;

        if(Stiff_type=="v"){
            K_des << k11, 0,
                    0, k22;
        }
        else{
            K_des << 1200, 0,
                    0, 1500 ;

        }


    } else {
        K_des << 20, 0, 0,
                0, 50, 0,
                0, 0, 100;
    }
    return K_des;
}


realtype VSDS::PosCheck(Vec x)
{
    Vec omega = Vec::Zero(K_);
    Vec delta = sigmascale_ * x_len_;
    Vec x_t(n_DOF_);
    if (x.size() == 3 && n_DOF_ == 2){
        x_t(0) = x(1);
        x_t(1) = x(2);
    }
    else{
        x_t = x;
    };

    for(int i=0;i <K_;i++){
        omega(i) = exp(-(1/(2*delta(i)*delta(i)))*(x_t-x_cen_.block(0,i,n_DOF_,1)).transpose()*(x_t-x_cen_.block(0,i,n_DOF_,1)));
    }
    realtype omega_sum = omega.sum();
    return omega.maxCoeff();
}

void VSDS::MotionRegenerate(Vec x)
{
    Vec x0_new(2);
    x0_new(0) = x(1);
    x0_new(1) = x(2);
    x0_ = x0_new;
    Mat xtemp;
    x_rec_ = x0_;
    Mat A_new;
    Mat B_new;
    A_ = A_new;
    B_ = B_new;
    xtemp = GetTempPoints(x0_, dt_, threshold_, false);
    GetViaPointsReduced(xtemp,x0_);
    GetABMatrix();
}


Vec VSDS::GetAttractor() {
    return att_ ;
}


