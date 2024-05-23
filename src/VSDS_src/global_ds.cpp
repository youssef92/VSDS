

#include "global_ds.h"

void global_ds::Read_SEDS_FROMFile(string MuFile, string SigmaFile, string PriorsFile, string A_gFile,string att_File,Mat &Mu, Mat &Sigma, Vec &Priors, Mat &A_g, int M, int K_SEDS,Vec &att){


    std::ifstream MuFileRead;
    std::ifstream SigmaFileRead;
    std::ifstream PriorsFileRead;
    std::ifstream A_gFileRead;
    std::ifstream b_gFileRead;
    std::ifstream attFileRead;

    MuFileRead.open(MuFile.c_str());
    if (MuFileRead){
        for(int i = 0; i < M; i++){
            for(int j=0; j<K_SEDS; j++){
                MuFileRead >> Mu(i,j);
            }
        }
        MuFileRead.close();
    }
    else {
        std::cout << "no MU file" << std::endl;
    }

    SigmaFileRead.open(SigmaFile.c_str());
    if (SigmaFileRead){
        for(int i = 0; i < M; i++){
            for(int j=0; j<M*K_SEDS; j++){
                SigmaFileRead >> Sigma(i,j);
            }
        }
        SigmaFileRead.close();
    }
    else {
        std::cout << "no Sigma file" << std::endl;
    }

    PriorsFileRead.open(PriorsFile.c_str());
    if (PriorsFileRead){
        for(int i = 0; i < K_SEDS; i++){
            PriorsFileRead >> Priors(i);
        }
        PriorsFileRead.close();
    }
    else {
        std::cout << "no Priors file" << std::endl;
    }

    A_gFileRead.open(A_gFile.c_str());
    if (A_gFileRead){
        for(int i = 0; i < M; i++){
            for(int j=0; j<M*K_SEDS; j++){
                A_gFileRead >> A_g(i,j);
            }
        }
        A_gFileRead.close();
    }
    else {
        std::cout << "no A_g file" << std::endl;
    }


    attFileRead.open(att_File.c_str());
    if (attFileRead){
        for(int i = 0; i < M; i++){
            attFileRead >> att(i);

        }
        attFileRead.close();
    }
    else {
        std::cout << "no Attractor file" << std::endl;
    }


}



global_ds::global_ds()  {

ROS_INFO("--------------Initializing VSDS-----------") ;

int K_SEDS=3 ;
ros::param::get("/K_SEDS", K_SEDS);
n_DOF_=2 ;
std::string DS_ModelName ;
ros::param::get("/VSDS_name",DS_ModelName);
std::string packPath = ros::package::getPath("VSDS");

std::string MuFile = packPath + "/config/"+ DS_ModelName+ "/Mu.txt";
std::string SigmaFile = packPath + "/config/"+ DS_ModelName+ "/Sigma.txt";
std::string PriorsFile = packPath + "/config/"+ DS_ModelName+ "/Priors.txt";
std::string A_gFile = packPath + "/config/"+ DS_ModelName+ "/A_g.txt";
std::string attFile = packPath + "/config/"+ DS_ModelName+ "/att.txt";



Mu_=Mat::Ones(n_DOF_,K_SEDS);
Priors_=Vec::Ones(K_SEDS);
A_g_=Mat::Ones(n_DOF_,n_DOF_*K_SEDS);
Sigma_=Mat::Ones(n_DOF_,n_DOF_*K_SEDS);
att_=Vec::Zero(n_DOF_) ;

Read_SEDS_FROMFile(MuFile, SigmaFile,  PriorsFile, A_gFile,attFile, Mu_, Sigma_, Priors_, A_g_,  n_DOF_,  K_SEDS,att_) ;


ROS_INFO("Initialized Global DS !!! ");

}

Mat global_ds::FindDampingBasis(Vec x)
{
    if (n_DOF_ == 2){
        Vec y(n_DOF_);
        y << 1,   -x(0)/(x(1)+eps());
        Mat B = Mat::Zero(2,2);
        B.col(0) = x/x.norm();
        B.col(1) = y/y.norm();
        return B;
    } else {
        Vector3d x_tmp = x;
        Vector3d y(1,1,(-x(0)-x(1))/(x(2)+eps()));
        Vector3d z;
        z = x_tmp.cross(y);
        Mat B = Mat::Zero(3,3);
        B.col(0) = x/x.norm();
        B.col(1) = y/y.norm();
        B.col(2) = z/z.norm();
        return B;
    }
}




Vec global_ds::global_ds_eval(Vec x) {

    int K_g = Priors_.size();
    Mat fg = Mat::Zero(n_DOF_,K_g);
    Vec fg_sum(n_DOF_);
    Vec beta = Vec::Zero(K_g);

    for (int i=0; i<K_g; i++)
    {
        Vec Mu_temp = Mu_.block(0,i,n_DOF_,1);
        Mat Sigma_temp = Sigma_.block(0,n_DOF_*i,n_DOF_,n_DOF_);
        beta(i) = exp(-0.5 * (x-Mu_temp).transpose() * Sigma_temp.inverse() * (x-Mu_temp)) / sqrt(pow(2*pi,n_DOF_) * abs(Sigma_temp.determinant()));
        beta(i) = beta(i) * Priors_(i);
    }
    beta = beta / beta.sum();

    for (int i=0; i<K_g; i++)
    {

        fg.col(i) = beta(i) * (A_g_.block(0,n_DOF_*i,n_DOF_,n_DOF_) * (x - att_));
    }
    fg_sum = fg.rowwise().sum();

    return fg_sum;

}

Vec global_ds::GetAttractor() {
    return att_ ;
}
