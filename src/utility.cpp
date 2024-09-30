#include "iostream"
#include "utility.h"

namespace special_math_functions {
	Mat Skew(const Vec &v) {
		if (v.size() == 1) {
			Mat out = Mat::Zero(2, 2);

			out << 0, -v(0),
				v(0), 0;
			return out;
        
		}
		else
		{
			Mat out = Mat::Zero(3, 3);
			out << 0, -v(2), v(1),
				v(2), 0, -v(0),
				-v(1), v(0), 0;
			return out;
		}


	}


    Eigen::Quaterniond quatMult(Eigen::Quaterniond q1, Eigen::Quaterniond q2) {
        Eigen::Quaterniond resultQ;
        resultQ.setIdentity();
        

        resultQ.w() = q1.w() * q2.w() - q1.vec().dot(q2.vec());
        resultQ.vec() = q1.w() * q2.vec() + q2.w() * q1.vec() + q1.vec().cross(q2.vec());
      
        

        return resultQ;
    }


    Vec quatlog(Eigen::Quaterniond q1){

        if(q1.vec() != Vec::Zero(3)){

            return acos(q1.w()) * (q1.vec()/q1.vec().norm()  ) ;

        }

       return  Vec::Zero(3) ;

        }

    Mat PushBackColumn(Mat m, Vec x)
    {
        if (m.cols() == 0){
            m = x;
        }
        else{
            m.conservativeResize(Eigen::NoChange, m.cols()+1);
            m.col(m.cols()-1) = x;
        }
        return m;
    }

    Mat PushBackColumnMat(Mat m, Mat x)
    {
        if (m.cols() == 0){
            m = x;
        }
        else{
            int n = x.cols();
            m.conservativeResize(Eigen::NoChange, m.cols()+n);
            for (int i=0; i<n; i++)
            {
                m.col(m.cols()-n+i) = x.col(i);
            }
        }
        return m;
    }

    realtype eps()
    {
        return 0.00000000000000001;
    }

    Mat Eig_Decomp_Sqrt(Mat M){

        EigenSolver<MatrixXd> ces;
        ces.compute(M);
        Mat V=  ces.eigenvectors().real()  ;
        Vec eig_val_sqrt= ces.eigenvalues().real() ;

        for (int i=0;i<eig_val_sqrt.size();i++){
            eig_val_sqrt(i)=sqrt(eig_val_sqrt(i)) ;
        }



        return V * eig_val_sqrt.asDiagonal() * V.inverse()  ;


    }

}

namespace StiffnessProfiles {
	 

	Mat StiffnessProfile_MainTask(realtype t) {

		//temp 
		realtype k = 1*GeneralFunction::smooth_transition_rising(t, 4.0, 0, 600);

		return  Mat::Identity(3, 3)*k;


	}

	Mat StiffnessProfile_NullSpace(realtype t) {

		//	realtype k = 10* t;
        realtype k = 1 * GeneralFunction::smooth_transition_rising(t, 4.0, 0, 16);

		return  Mat::Identity(7, 7)*k;


	}

}

namespace GeneralFunction {

using namespace std;

	void Write_To_File(string file_name, vector<vector<float> > Input_Data) {
		ofstream myfile;
		myfile.open(file_name);

        if (!myfile)
        {
            cout << "No file found: " << file_name << endl;
            return;
        }

		for (int i = 0; i < Input_Data.size(); i++)
		{
			for (int j = 0; j < Input_Data[i].size(); j++)
			{

				if (j == Input_Data[i].size() - 1) {
					myfile << Input_Data[i][j] << endl;



				}
				else {

                    myfile << Input_Data[i][j] << " ";
				}


			}
		}
		myfile.close();
	}

	realtype smooth_transition_rising(realtype t, realtype t_max, realtype t_min, realtype scale){

        realtype alpha_min = 0.0;
        realtype alpha;

        if (t>t_max)
            alpha = 0;
        else if(t<t_min)
            alpha = 1;
        else
            alpha = (0.5*(1 + cos(pi / (t_max - t_min)*(t - t_min)))*(1 - alpha_min) + alpha_min);


         alpha = 1 - alpha;

        return  scale*alpha;

	}



    realtype smooth_transition_fall(realtype E, realtype E_max, realtype E_min){
        realtype alpha_min = 0.0;
        realtype alpha;

        if (E>E_max)
            alpha = 0;
        else if(E<E_min)
            alpha = 1;
        else
            alpha = (0.5*(1 + cos(pi / (E_max - E_min)*(E - E_min)))*(1 - alpha_min) + alpha_min);


            return  alpha;
	}
	
    void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
    {
        unsigned int numRows = matrix.rows();
        unsigned int numCols = matrix.cols()-1;

        if( colToRemove < numCols )
            matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

        matrix.conservativeResize(numRows,numCols);
    }

    realtype mod(realtype x, realtype y){
        return x - floor(x/y)*y ;
    }

    void  Vec2double(Vec a,double q[]) {

        int s = a.size();
        for (int i = 0; i < s; i++) {

            q[i] = a(i);
        }

    }

    int loadVectorMatrixFromFile (std::string fileName, int cols, vector<vector<float>> &outMat)
    {
        ifstream in(fileName.data());
        if (!in)
        {
            cout << "No file found: " << fileName << endl;
            return -1;
        }
        int counter = 0;
        while (!in.eof())
        {
            outMat.push_back( vector <float>() );
            for (int j = 0; j < cols; ++j)
            {
                double readf;
                in >> readf;
                outMat[counter].push_back(readf);
            }
            counter++;
        }
        outMat.pop_back();
        in.close();
        return 0;
    }

    Vec SaturationFunc(Vec inp,float max){

        int size=inp.rows() ;
        Vec out(size) ;

        for(int i=0;i<size;i++){

            if(inp[i]>max){
                out[i]=max ;

            }
            else if(inp[i]< -max){
                out[i]=-max ;

            }
            else{
                out[i]=inp[i] ;

            }

        }
        return out ;

    }

    void saveVectorMatrixToFile (string fileName, vector < vector <float> > outMat)
    {
        ofstream out(fileName.data());
        if (!out)
        {
            cout << "No file found: " << fileName << endl;
            return;
        }
        int rows = (int)outMat.size();
        int cols = (int)outMat[0].size();
        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                out << outMat[i][j] << "\t";
            }
            out << endl;
        }
        out.close();
        return;
    }

    float getSquaredDistance(float a[3], float b[3]){
        return (a[0]-b[0])*(a[0]-b[0]) +
                (a[1]-b[1])*(a[1]-b[1]) +
                (a[2]-b[2])*(a[2]-b[2]) ;
    }

    float low_pass(float signal, float prev_filt, float cutt_off, float cycle_time){

        return   ( ( signal*cutt_off*cycle_time)+ prev_filt )/(1+ cutt_off*cycle_time) ;

    }

    Vec low_pass(Vec signal, Vec prev_filt, float cutt_off, float cycle_time){

        Vec out=Vec::Zero(signal.size()) ;
        for (int i=0;i<signal.size();i++){

            out(i)=low_pass(signal(i),prev_filt(i),cutt_off,cycle_time) ;
        }
        return out ;


    }




}
