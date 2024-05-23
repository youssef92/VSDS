#pragma once
#include "eigen3/Eigen/Dense"
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include "vector"
#include <fstream>
#include "iostream"
#include <ros/ros.h>
#include <ros/package.h>
using namespace std;
using namespace Eigen ;
#ifndef PASSIVE_DS_TYPEDEFS_H
#define PASSIVE_DS_TYPEDEFS_H


	//#ifdef USE_DOUBLE_PRECISION
	typedef double realtype;
	/*#else
	typedef float realtype;
	#endif*/
	typedef Eigen::Matrix<realtype, Eigen::Dynamic, Eigen::Dynamic> Mat;
	typedef Eigen::Matrix<realtype, Eigen::Dynamic, 1> Vec;

#endif
         const realtype pi=         3.141592653589793238462643383279502884 /* pi */;

         template <typename T> int sgn(T val) {
             return (T(0) < val) - (val < T(0));
         }

namespace special_math_functions
{

	Mat Skew(const Vec &v);
    Quaterniond quatMult(Eigen::Quaterniond q1, Eigen::Quaterniond q2) ;
    Vec quatlog(Eigen::Quaterniond q1) ;
    Mat PushBackColumn(Mat m, Vec x);
    Mat PushBackColumnMat(Mat m, Mat x);
    Mat Eig_Decomp_Sqrt(Mat M) ;
    realtype eps();
}

namespace StiffnessProfiles
{
	Mat StiffnessProfile_MainTask(realtype t);
	Mat StiffnessProfile_NullSpace(realtype t);
}
namespace GeneralFunction {
	void Write_To_File(string file_name, vector<vector<float> > Input_Data);
	realtype smooth_transition_rising(realtype t, realtype t_max, realtype t_min, realtype scale);
	realtype smooth_transition_fall(realtype E, realtype E_max, realtype E_min);
    void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove) ;
    realtype mod(realtype x, realtype y) ;
    void  Vec2double(Vec a,double q[]) ;
    Vec double2Vec(double q[]) ;
    int loadVectorMatrixFromFile (std::string fileName, int cols, vector<vector<float>> &outMat) ;
    Vec SaturationFunc(Vec inp,float max) ;
    void saveVectorMatrixToFile (string fileName, vector < vector <float> > outMat) ;
    float getSquaredDistance(float a[3], float b[3]) ;
    float low_pass(float signal, float prev_filt, float cutt_off, float cycle_time) ;
    Vec low_pass(Vec signal, Vec prev_filt, float cutt_off, float cycle_time) ;

    template <size_t size_x, size_t size_y>

    Mat double2Mat(double(&T)[size_x][size_y])
    {
        int rows = size_x;
        int cols = size_y;
        Mat out = Mat::Zero(rows, cols);
        for (int i = 0; i<rows; i++)
        {
            for (int j = 0; j<cols; j++)
            {
                out(i, j) = T[i][j];
            }

        }
        return out;
    }



}
