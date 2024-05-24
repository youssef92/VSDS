/*
 * Author: Youssef Michel
 * Description: This scipt implements the main control loop of the VSDS control law,
 *              handling the communication with a Kuka Robot via FRI, recieving sensor
 *              readings and sending force commands to the robot
 *
 * Usage: This code can be compiled and run using a standard C++ compiler.
 *        Example compilation command: g++ example.cpp -o example
 *        Example execution command: ./example
 *
 * Notes:
 *  - This should be adapted to integrate VSDS with your own robot interface
 */




#include "boost/filesystem.hpp"

#include "MotionGeneration.h"
#include "PassiveControl.h"
#include "Utility_fri.h"
#include "VSDS_base.h"

using namespace boost::filesystem;



int main(int argc, char *argv[])
{

    ros::init(argc, argv, "VSDS_main");
    ros::NodeHandle nh;
    string DS_ModelName ;

    string VS_type="qp" ; // org , qp , fd
    nh.getParam("/VS_type",VS_type) ;

    FastResearchInterface	*FRI;
    std::string packPath = ros::package::getPath("VSDS");
    std::cout << packPath << "\n";

    nh.getParam("/VSDS_name",DS_ModelName);
    cout <<"----------------------Started VSDS ROS Node with DS: "<<packPath + "/config/"+ DS_ModelName+ "/"<<" ---------------- "<<endl ;
    fprintf(stdout, "You may need superuser permission to run this program.\n");
    fflush(stdout);

    bool					Run							=	true;
    char					c							=	0 ;
    int i = 0;
    int						ResultValue					=	0;
    float    CartStiffnessValues[FRI_CART_VEC],
            CartDampingValues[FRI_CART_VEC];
    float currentCartPose[FRI_CART_FRM_DIM] ;
    float CartPose_init[FRI_CART_FRM_DIM] ;
    float currentJointPosition[LWR_JNT_NUM];
    float desiredCartPose [FRI_CART_FRM_DIM] ;
    float CommandedForcesAndTorques [NUMBER_OF_CART_DOFS] ;

    FRI = new FastResearchInterface((packPath + "/data/Control-FRI-Driver_2ms.init").c_str());
    fprintf(stdout, "Check .\n");
    fprintf(stdout, "OK-OK \n");
    fflush(stdout);
    std::string initJointFile    = packPath + "/data/InitialJointPosition_LASA.txt";

    while (Run)
    {

        std::string DataPathEnding;
        cout<<"\n"<<"Which test you would like to do?\n"<<"\n";
        cin>>DataPathEnding;
        cout<<"\n\n";
        printf("---------------------------------------------------------------------------------------\n");
        printf("Press     q  for exit thisss program\n");
        printf("          h  for start the joint position controller and go to home position\n");
        printf("          g  for the gravity compensation mode\n");
        printf("          w  for VSDS \n");
        printf("---------------------------------------------------------------------------------------\n\n");
        printf("Please press any key...\n");

        c	=	WaitForKBCharacter(NULL);

        printf("\n\n\n");

        switch (c)
        {
        case 'h':
        case 'H':{

            printf("Going to home position... (please wait)\n");
            RunJointTrajectory(FRI, initJointFile);
            break ;
        }

        case 'g':
        case 'G':{

            JointGravityCompensation(FRI, 60) ;
            break ;}


        case 'w':
        case 'W':{
            printf("VSDS control loop \n");
            drdSleep(1);

            if(startCartImpedanceCtrl(FRI, currentJointPosition)==0){


                float **ptr_jacobian;
                ptr_jacobian = new float *[7];
                for(int i = 0; i <7; i++){
                    ptr_jacobian[i] = new float[7];
                }

                float  dt=FRI->GetFRICycleTime() ;
                FRI->GetMeasuredJointPositions(currentJointPosition);
                FRI->GetMeasuredCartPose(currentCartPose);
                FRI->GetMeasuredCartPose( CartPose_init);
                Vec x0=GetTranslation(currentCartPose) ; // 3*1 Vec
                Vec q_prev=float_2_Vec(currentJointPosition,7)  ;
                Vec x = x0;
                int M = 2;
                realtype x_c = x0(0);

                string VS_type="qp" ; // org , qp , fd
                nh.getParam("/VS_type",VS_type) ;

                std::unique_ptr<VSDS_base> MyVSDS;

                if (VS_type == "qp") {
                    MyVSDS = std::make_unique<VSDS_qp>(x0.tail(2));
                }
                else if (VS_type == "fd") {
                    MyVSDS = std::make_unique<VSDS_qp>(x0.tail(2));
                }
                else if (VS_type == "PI") {
                    MyVSDS = std::make_unique<VSDS_picds>(x0.tail(2));
                }

                else {
                   MyVSDS = std::make_unique<VSDS_org>(x0.tail(2));
                }

                global_ds MyGlobalDS =global_ds(M) ;
                PassiveClosedLoopControl *my_control_passive  ;
                my_control_passive = new PassiveClosedLoopControl( M,x0.tail(2),MyVSDS->GetAttractor() );


                for (i = 0; i < NUMBER_OF_CART_DOFS; i++)
                {
                    if (i<3){
                        CartStiffnessValues[i] =(float)0.01;
                        CartDampingValues[i]=(float)0.0;
                    }
                    else{
                        CartStiffnessValues[i] =(float)200.0 ;
                        CartDampingValues[i]=(float)0.7;

                    }
                    CommandedForcesAndTorques	[i]	=	(float)0.0;
                }
                FRI->SetCommandedCartStiffness(  CartStiffnessValues);
                FRI->SetCommandedCartDamping(CartDampingValues);
                FRI->SetCommandedCartForcesAndTorques( CommandedForcesAndTorques);

                for(int i=0;i<12;i++){
                    desiredCartPose[i]=CartPose_init[i] ;
                }
                FRI->SetCommandedCartPose(desiredCartPose);


                double t = 0;
                unsigned CycleCounter=0 ;
                Vec q_dot_filt_prev= Vec::Zero(7) ;

                int done=1 ;

                int DISTURBANCE_FLAG= 0 ;
                ros::param::get("/Disturbance_Flag", DISTURBANCE_FLAG) ;
                Vec F_dist=Vec::Zero(3) ;
                realtype F_distMag, t_init, t_final ;

                if(DISTURBANCE_FLAG) {
                    ROS_WARN("Will Apply Disturbance !!") ;

                    if ( ! ros::param::get("/F_disturbance", F_distMag)) {
                      ROS_WARN("Could not read Disturbance Force Magn. !!") ;
                    }

                    if ( ! ros::param::get("/t_init_dist", t_init)) {
                      ROS_WARN("Could not read Disturbance T_init !!") ;
                    }

                    if ( ! ros::param::get("/t_final_dist", t_final)) {
                      ROS_WARN("Could not read Disturbance T_final !!") ;
                    }
                }


                while ((FRI->IsMachineOK()) && ((float)CycleCounter * FRI->GetFRICycleTime()<20.0 ) && done==1  ){
                    FRI->WaitForKRCTick();
                    if (dhdKbHit() && dhdKbGet()=='q') done = -1;

                    //--------------------------------------------------Get Current robot Data----------------------------------------------------------
                    FRI->GetMeasuredCartPose(currentCartPose);
                    FRI->GetMeasuredJointPositions(currentJointPosition);
                    FRI->GetCurrentJacobianMatrix (ptr_jacobian);

                    x=GetTranslation(currentCartPose) ;
                    Vec q=float_2_Vec(currentJointPosition,7) ;

                    Mat Jac_temp=Convert_Jacobian_2Mat(ptr_jacobian) ;
                    Mat Jacobian_Matrix_tool=Jac_temp ;
                    Jacobian_Matrix_tool.row(3)=Jac_temp.row(5) ;
                    Jacobian_Matrix_tool.row(5)=Jac_temp.row(3) ;
                    Mat Rot_mat=GetRotationMatrix(currentCartPose) ;
                    Mat Jacobian_Matrix_world=Tool_2_World_Jacobian(Jacobian_Matrix_tool,Rot_mat) ;

                    // Get Robot Velocity
                    Vec q_dot=(q-q_prev)/dt ;
                    q_prev=q ;
                    Vec q_dot_filt=low_pass( q_dot, q_dot_filt_prev,10,FRI->GetFRICycleTime() ) ;
                    q_dot_filt_prev=q_dot_filt ;
                    Vec x_dot_estimated=Jacobian_Matrix_world.block(0,0,3,7)*q_dot_filt ;

                    // Get Desired Force
                    Vec Force_d =  MyVSDS->GetDesiredForce(x,x_dot_estimated,q);
                    Vec x_dot_gds = MyGlobalDS.global_ds_eval(x.tail(2)) ;

                    // Get Control Force Passive
                    Vec F_x_passive=my_control_passive->GetControlForce_Passive(Force_d,x_dot_estimated, x_c, x,FRI->GetFRICycleTime()) ;
                    Vec Impedance_Force_tool=Rot_mat.transpose()*F_x_passive ;

                    if(DISTURBANCE_FLAG){

                            F_dist=ApplyDisturbance(F_distMag,t, t_init,t_final,x_dot_gds) ;
                            F_dist=Rot_mat.transpose()*F_dist ;
                    }

                    for(int i=0;i <3;i++){

                        CommandedForcesAndTorques[i]=Impedance_Force_tool[i] + F_dist[i] ;

                    }

                    desiredCartPose[3]=currentCartPose[3]   ;
                    desiredCartPose[7]=currentCartPose[7]   ;
                    desiredCartPose[11]=currentCartPose[11] ;
                    FRI->SetCommandedCartPose(desiredCartPose);
                    FRI->SetCommandedCartForcesAndTorques(CommandedForcesAndTorques) ;
                    t += dt;
                    CycleCounter++;
                }


                printf("Motion Completed.\n");

                if (ResultValue != EOK)
                {
                    fprintf(stderr, "An error occurred during stopping the robot...\n");
                }
                else
                {
                    fprintf(stdout, "Robot successfully stopped.\n");
                }

                sleep(2);


            }

            break;
        }


        default:
            printf("This key is not supported yet...\n");
            break;

        case 'q':
        case 'Q':

            ResultValue=FRI->StopRobot() ;
            delete FRI;
            sleep(3);
            return(EXIT_SUCCESS) ;
        }

    }

}




