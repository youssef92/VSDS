#ifndef VSDSTASKPLANNER_H
#define VSDSTASKPLANNER_H

#include "Utility.h"
#include "MotionGeneration.h"

#include "PassiveControl.h"
#include "min_jerk.h"
#include "Utility_fri.h"

class vsds_task_planner{

private:

    FastResearchInterface	*FRI_;

    float **ptr_jacobian_;

    float  JointStiffnessValues_[LBR_MNJ],
            CommandedJointTorques_[LBR_MNJ],
            JointDampingValues_[LBR_MNJ],
            EstimatedExternalCartForcesTorques_[FRI_CART_VEC],
            EstimatedExternalJointTorques_[LBR_MNJ];
    float currentCartPose_[FRI_CART_FRM_DIM] ;
    float currentJointPosition_[LWR_JNT_NUM];
    void init_datalogging()  ;

public:

    vsds_task_planner() ;
    int init(FastResearchInterface	*FRI,ros::NodeHandle nh,string DSname)   ;

    void  run()  ;
    void push_data_toVector() ;
    void save_data_toFile()   ;





};


#endif
