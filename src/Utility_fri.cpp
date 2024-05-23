


#include "Utility_fri.h"

int startCartImpedanceCtrl(FastResearchInterface *fri, float *commCartPose){
    unsigned int controlScheme = FastResearchInterface::CART_IMPEDANCE_CONTROL;
    int resultValue;
    if(fri->GetCurrentControlScheme() != controlScheme || !fri->IsMachineOK()){
        // Stop

        fri->StopRobot();

        fri->GetMeasuredCartPose(commCartPose);
        fri->SetCommandedCartPose(commCartPose);

        // Restart
        resultValue	= fri->StartRobot(controlScheme);

        if (resultValue != EOK){
            std::cout << "An error occurred during starting up the robot..." << std::endl;
            return -1;
        }
    }
    return 0;
}

Vec ApplyDisturbance(realtype F_distMag,realtype t,realtype t_init,realtype t_final, Vec Des_vel) {

    Vec dir_vel= Des_vel/Des_vel.norm() ;
    Vec dir_force= Vec::Zero(2) ;
    Vec dist_force= Vec::Zero(2) ;
    Vec dist_force_tot= Vec::Zero(3) ;

    dir_force(0)= dir_vel(1) ;
    dir_force(1)= -dir_vel(0) ;

    if(t>t_init && t<t_final) {
        dist_force=-F_distMag * dir_force ;
    }

    dist_force_tot(1)=dist_force(0) ;
    dist_force_tot(2)=dist_force(1) ;

    return dist_force_tot ;
}


int JointGravityCompensation(FastResearchInterface *FRI, float tot_time) {

    printf("Gravity compensation... (please wait)\n");

    float currentCartPose[FRI_CART_FRM_DIM] ;
    float currentJointPosition[LWR_JNT_NUM];
    float  JointStiffnessValues[LBR_MNJ] ;
    std::string packPath = ros::package::getPath("VSDS");
    std::string demoCartTrajFile = packPath + "/data/DemonstratedTrajectory.txt";
    std::string demoJointTrajFile= packPath+ "/data/DemonstratedJointTrajectory.txt" ;
    std::vector< std::vector <float> > demo_;
    std::vector<std::vector<float>>demo_Joints ;

    if(startJointImpedanceCtrl(FRI, currentJointPosition)==0){

        printf("Gravity compensation Started \n");

        // Set stiffness
        for(int i = 0; i < LWR_JNT_NUM; i++){
            JointStiffnessValues[i] = (float)MIN_STIFFNESS; // max stiffness 0-2 -> 2000.0, max 3-5 200.0
        }

        FRI->SetCommandedJointStiffness(JointStiffnessValues);
        int it_= 0 ;

        int LoopValue = int(tot_time / FRI->GetFRICycleTime());
        demo_.clear();
        demo_Joints.clear();

        for(int i = 0; i<LoopValue; ++i){
            demo_.push_back( std::vector <float>() );
            demo_Joints.push_back(std::vector<float>());
            for(int j=0; j<FRI_CART_FRM_DIM; ++j){
                demo_[i].push_back(0.0);

            }
            for (int j=0; j<7;++j){
                demo_Joints[i].push_back(0.0);
            }
        }
        int done_g=1 ;

        while (FRI->IsMachineOK()  && (it_<LoopValue)&& (done_g == 1)  ) {

            FRI->WaitForKRCTick();

           if (dhdKbHit() && dhdKbGet()=='q') done_g = -1;
            FRI->GetMeasuredJointPositions(currentJointPosition);
            FRI->SetCommandedJointPositions(currentJointPosition);
            FRI->GetMeasuredCartPose(currentCartPose);
            FRI->SetCommandedCartPose(currentCartPose);

            for(int i=0; i<FRI_CART_FRM_DIM; ++i){
                demo_[it_][i] = currentCartPose[i];
            }
            for (int  i=0;i<NUMBER_OF_JOINTS;++i){
                demo_Joints[it_][i]=currentJointPosition[i] ;
            }

            it_++;


        }
        cout<<"Joint Angles: \n "  << float_2_Vec(currentJointPosition,7)*(180/PI)	<<endl ;
        cout<<"Current Pose: \n "  << GetTranslation(currentCartPose)	<<endl ;

        saveVectorMatrixToFile(demoCartTrajFile, demo_);
        saveVectorMatrixToFile(demoJointTrajFile,demo_Joints);

    }

    return 1;

}


int CartesianGravityCompensation_YZ(FastResearchInterface *FRI, float tot_time) {
    float   CartStiffnessValues[FRI_CART_VEC],
            CartDampingValues[FRI_CART_VEC];
    float currentCartPose[FRI_CART_FRM_DIM] ;
    float CommandedForcesAndTorques [NUMBER_OF_CART_DOFS] ;
    float x_d_demonstration [FRI_CART_FRM_DIM] ;
    float currentJointPosition[LWR_JNT_NUM];

    FRI->GetMeasuredCartPose(currentCartPose);

    for(int i=0; i<FRI_CART_FRM_DIM; i++){
        x_d_demonstration[i]=currentCartPose[i] ;
    }
    for (int i = 0; i < NUMBER_OF_CART_DOFS; i++)
    {
        if (i<3){
            CartStiffnessValues[i] =(float)1500;
            CartDampingValues[i]=(float)0.7;
        }
        else{
            CartStiffnessValues[i] =(float)200.0 ;
            CartDampingValues[i]=(float)0.7;
        }
        CommandedForcesAndTorques	[i]	=	(float)0.0;
    }
    CartStiffnessValues[1] =(float)0.01;
    CartStiffnessValues[2] =(float)0.01;

     FRI->SetCommandedCartPose(x_d_demonstration);

    int done=1 ;
    int CycleCounter=0 ;

    printf("       please enter Demo Number \n");
    char n ;
    n	=	WaitForKBCharacter(NULL);
    printf("\n\n\n");
    std::string ss; ss.push_back(n);
    std::string packPath = ros::package::getPath("VSDS");
    std::string x_file_demo = packPath + "/data/DemoDrill_1"+ss+  "_X.txt" ;
    std::string q_file_demo = packPath + "/data/DemoDrill_1"+ss+  "_Q.txt" ;
    std::vector< std::vector <float> > x_rob_Vector ;
    std::vector< std::vector <float> > q_rob_Vector ;

    if(startCartImpedanceCtrl(FRI, currentJointPosition)==0){

         printf(" Cartesian Compensation Started-For Saving \n");
    while ((FRI->IsMachineOK()) && ((float)CycleCounter * FRI->GetFRICycleTime()<60.0 ) && done==1  ){

        FRI->WaitForKRCTick();
        if (dhdKbHit() && dhdKbGet()=='q') done = -1;

        FRI->GetMeasuredCartPose(currentCartPose);
        Vec x=GetTranslation(currentCartPose) ;

        FRI->GetMeasuredJointPositions(currentJointPosition);


        x_d_demonstration[7]=currentCartPose[7]   ;
        x_d_demonstration[11]=currentCartPose[11] ;

        FRI->SetCommandedCartPose(x_d_demonstration);
        x_rob_Vector.push_back( std::vector <float>() );
        q_rob_Vector.push_back(std::vector <float>());

        x_rob_Vector[CycleCounter].push_back(x(0));
        x_rob_Vector[CycleCounter].push_back(x(1));
        x_rob_Vector[CycleCounter].push_back(x(2));
        for (int i=0 ;i <7;i++){
            q_rob_Vector[CycleCounter].push_back(currentJointPosition[i]);
        }

        CycleCounter++ ;
    }
    saveVectorMatrixToFile(x_file_demo,x_rob_Vector);
    saveVectorMatrixToFile(q_file_demo,q_rob_Vector);
    printf(" Saved Demo !!\n");


}




}


int CartesianGravityCompensation(FastResearchInterface *FRI, float tot_time) {

    float   CartStiffnessValues[FRI_CART_VEC],
            CartDampingValues[FRI_CART_VEC];
    float currentCartPose[FRI_CART_FRM_DIM] ;
    float CommandedForcesAndTorques [NUMBER_OF_CART_DOFS] ;
    float x_d_demonstration [FRI_CART_FRM_DIM] ;
    float currentJointPosition[LWR_JNT_NUM];

    FRI->GetMeasuredCartPose(currentCartPose);

    for(int i=0; i<FRI_CART_FRM_DIM; i++){
        x_d_demonstration[i]=currentCartPose[i] ;
    }

    for (int i = 0; i < NUMBER_OF_CART_DOFS; i++)
    {
        if (i<3){
            CartStiffnessValues[i] =(float)1500;
            CartDampingValues[i]=(float)0.7;
        }
        else{
            CartStiffnessValues[i] =(float)200.0 ;
            CartDampingValues[i]=(float)0.7;

        }
        CommandedForcesAndTorques	[i]	=	(float)0.0;
    }
    FRI->SetCommandedCartPose(x_d_demonstration);

    int done=1 ;
    int CycleCounter=0 ;
    std::vector< std::vector <float> > x_robgrav_Vector;
    std::vector< std::vector <float> > q_robgrav_Vector;


    if(startCartImpedanceCtrl(FRI, currentJointPosition)==0){
        CartStiffnessValues[1] =(float)0.01;
        CartStiffnessValues[2] =(float)0.01;
        printf(" Cartesian Compensation Started -Free \n");
        while ((FRI->IsMachineOK()) && ((float)CycleCounter * FRI->GetFRICycleTime()<tot_time) && done==1  ){

            FRI->WaitForKRCTick();
            if (dhdKbHit() && dhdKbGet()=='q') done = -1;

            FRI->GetMeasuredCartPose(currentCartPose);
            Vec x=GetTranslation(currentCartPose) ;

            FRI->GetMeasuredJointPositions(currentJointPosition);


            x_d_demonstration[7]=currentCartPose[7]   ;
            x_d_demonstration[11]=currentCartPose[11] ;

            FRI->SetCommandedCartPose(x_d_demonstration);
            x_robgrav_Vector.push_back( std::vector <float>() );
            q_robgrav_Vector.push_back(std::vector <float>());

            x_robgrav_Vector[CycleCounter].push_back(x(0));
            x_robgrav_Vector[CycleCounter].push_back(x(1));
            x_robgrav_Vector[CycleCounter].push_back(x(2));
            for (int i=0 ;i <7;i++){
                q_robgrav_Vector[CycleCounter].push_back(currentJointPosition[i]);
            }

            CycleCounter++ ;
        }
    }

    return 1 ;

}

vector< std::vector <float> > Kuka_MoveCartesian_MinJerk(FastResearchInterface *FRI, float tot_time, float dt, Vec x_df) {
    printf("Moving to Desired Cartesian Position with a minimum Jerk Trajectorty  \n");
    float CartStiffnessValues[FRI_CART_VEC],  CartDampingValues[FRI_CART_VEC], CommandedForcesAndTorques[FRI_CART_VEC] ;
    float currentCartPose[FRI_CART_FRM_DIM], commCartPose[FRI_CART_FRM_DIM],currentjointpos[7];
    FRI->GetMeasuredCartPose(currentCartPose);
    for(int i=0; i<FRI_CART_FRM_DIM; i++){
        commCartPose[i]=currentCartPose[i] ;
    }
    std::vector< std::vector <float> > X_pos_Vec;

    MinJerk *min_jerk_profile;
    min_jerk_profile = new MinJerk(dt,tot_time, GetTranslation(currentCartPose),x_df );


    for (int i = 0; i < 6; i++)
    {
        if (i==0 || i==1 || i==2){
            CartStiffnessValues[i] = (float)2500;
            CartDampingValues[i]=(float)0.7;
        }
        else{
            CartStiffnessValues[i] =(float)200.0 ;
            CartDampingValues[i]=(float)0.7;
        }
        CommandedForcesAndTorques	[i]	=	(float)0.0;
    }

    CartStiffnessValues[1] = (float)2000.0;

    FRI->SetCommandedCartStiffness(CartStiffnessValues);
    FRI->SetCommandedCartDamping(CartDampingValues);
    FRI->SetCommandedCartForcesAndTorques( CommandedForcesAndTorques);
    FRI->SetCommandedCartPose(commCartPose);

    Vec err_x = Vec::Ones(3) ; realtype tol=0.01 ; realtype Max_time=55 ;  realtype t=0 ; Vec x=GetTranslation(currentCartPose);
    Vec x_d=GetTranslation(currentCartPose) ; int wv_index=0 ; Mat R=GetRotationMatrix(currentCartPose) ;
    Mat R_d=R ;
    R_d << -1,0,0,
            0,1,0,
            0,0,-1;
    R_d=R ;

    //
    int done_g=1 ;
    while ( (FRI->IsMachineOK()) && err_x.norm()>tol && t<Max_time  ){
        // while ( (FRI->IsMachineOK()) &&  t<Max_time  && done_g==1 ){

        FRI->WaitForKRCTick();
        X_pos_Vec.push_back( std::vector <float>() );
        min_jerk_profile->Update(t);
        x_d=min_jerk_profile->GetDesiredPos() ;
        commCartPose[3]=x_d(0);commCartPose[7]=x_d(1);commCartPose[11]=x_d(2);
        if (dhdKbHit() && dhdKbGet()=='q') done_g = -1;


        FRI->GetMeasuredCartPose(currentCartPose);
        for (int i = 0; i < FRI_CART_VEC; i++){
            if(i<3){
                //   CartStiffnessValues[i] =  GeneralFunction::smooth_transition_rising(t, 6.0, 0, 300);    ;
            }

            else{
                //    CartStiffnessValues[i] = GeneralFunction::smooth_transition_rising(t, 5.0, 0, 200);
                CartStiffnessValues[i] =200 ;
            }

        }

        Set_Desired_Pose_FRI(R_d,x_d,commCartPose);
        FRI->SetCommandedCartStiffness(CartStiffnessValues);


        FRI->GetMeasuredCartPose(currentCartPose);
        FRI->GetMeasuredJointPositions(currentjointpos) ;
        Vec q=float_2_Vec(currentjointpos,7) ;
        x=GetTranslation(currentCartPose) ; R=GetRotationMatrix(currentCartPose) ;


        err_x=x-x_df ;
        FRI->SetCommandedCartPose(commCartPose);
        for(int j=0; j<3; ++j){
            X_pos_Vec[wv_index].push_back(x(j));
        }

        for(int j=0; j<3; ++j){
            X_pos_Vec[wv_index].push_back(x_d(j));
        }
        for(int j=0; j<7; ++j){
            X_pos_Vec[wv_index].push_back(q(j));
        }


        wv_index++ ;
        t+=dt ;

    }

    FRI->GetMeasuredCartPose(currentCartPose);

    FRI->SetCommandedCartPose(currentCartPose);

    std::cout << "Cartesian Position reached..." << std::endl;
    cout<<"Current Cartesisian position: \n "  << GetTranslation(currentCartPose)	<<endl ;
    cout<<"Desired Final GOal: \n "  << x_df	<<endl ;

    return X_pos_Vec ;

}


void Set_Desired_Pose_FRI(Mat Rot_d,Vec x_d,float *CartPose_d){

    int tot_count=0 ;

    for (int i=0;i<3;i++){

        int j=0 ;

        while(j<3){

            if ( !(tot_count==3 || tot_count==7|| tot_count==11 ) ){

                CartPose_d[tot_count]=Rot_d(i,j) ;

                j=j+1 ;

            }
            tot_count++ ;

        }



    }
    CartPose_d[3]=x_d(0) ;
    CartPose_d[7]=x_d(1) ;
    CartPose_d[11]=x_d(2) ;



}



int HD_gravityCompensation()
{


    if (dhdSetForceAndTorqueAndGripperForce (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) < DHD_NO_ERROR) {
        //   printf ("error: cannot set force (%s)\n", dhdErrorGetLastStr());
        return -1 ;
    }
    return 1 ;

}


int startJointImpedanceCtrl(FastResearchInterface *fri, float *commJointPosition){
    unsigned int controlScheme = FastResearchInterface::JOINT_IMPEDANCE_CONTROL;
    int resultValue;
    if(fri->GetCurrentControlScheme() != controlScheme || !fri->IsMachineOK()){
        // Stop
        //cout << fri->WaitForKRCTick()<< endl;
        fri->StopRobot();


        fri->GetMeasuredJointPositions(commJointPosition);
        fri->SetCommandedJointPositions(commJointPosition);

        // Restart
        resultValue	= fri->StartRobot(controlScheme);
        if (resultValue != EOK){
            std::cout << "An error occurred during starting up the robot..." << std::endl;
            return -1;
        }
    }
    return 0;
}







Mat GetRotationMatrix(float *CartPose){

    int tot_count=0 ;
    Mat Rot_Mat(3,3) ;
    for (int i=0;i<3;i++){

        int j=0 ;

        while(j<3){

            if ( !(tot_count==3 || tot_count==7|| tot_count==11 ) ){

                Rot_Mat(i,j)=CartPose[tot_count] ;
                j=j+1 ;

            }
            tot_count++ ;

        }



    }

    return Rot_Mat ;


}

Vec GetTranslation(float *CartPose){

    int tot_count=0 ;
    Vec transl(3) ;
    int i=0 ;

    while(i<3){

        if  (tot_count==3 || tot_count==7|| tot_count==11 ) {

            transl(i)=CartPose[tot_count] ;
            i++ ;


        }
        tot_count++ ;

    }

    return transl ;
}

Mat Convert_Jacobian_2Mat(float ** Jacobian){

    Mat output(6,7) ;

    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 7; j++)
        {
            //JacobianMatrix[i][j]  =   this->ReadData.data.FRIJacobianMatrix[(i==3)?(5):((i==5)?(3):(i))*NUMBER_OF_JOINTS+j];
            output(i,j)    =   Jacobian[i][j];
        }
    }

    return output ;

}

Mat Convert_MassMat_2Mat(float ** Mass){

    Mat output(7,7) ;

    for (int i = 0; i < 7; i++)
    {
        for (int j = 0; j < 7; j++)
        {
            //JacobianMatrix[i][j]  =   this->ReadData.data.FRIJacobianMatrix[(i==3)?(5):((i==5)?(3):(i))*NUMBER_OF_JOINTS+j];
            output(i,j)    =   Mass[i][j];
        }
    }

    return output ;

}



Vec float_2_Vec(float *inp,int size){

    Vec out(size) ;

    for (int i=0;i<size;i++){

        out(i)=inp[i] ;
    }
    return out ;

}

Mat Tool_2_World_Jacobian(Mat Jac_tool,Mat Rot){
    Mat Jac_world=Jac_tool ;
    Mat zeros=MatrixXd::Zero(3,3);
    Mat Trans(6,6) ;
    Trans<<Rot,zeros ,
            zeros,Rot ;
    Jac_world=Trans *Jac_tool ;

    return Jac_world ;



}

Vec Quart_Orient_err(Mat R_a, Mat R_d){

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
    return error ;
}
