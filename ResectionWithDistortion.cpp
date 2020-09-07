#include "ResectionWithDistortion.h"
#include "ResectionResidualWithDistortion.h"
#include "ceres\ceres.h"
#include "Eigen\Eigen"

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
using namespace Eigen;

ResectionWithDistortion::ResectionWithDistortion()
{
}


ResectionWithDistortion::~ResectionWithDistortion()
{
}

void ResectionWithDistortion::doResection()
{
	//double u[8] = { 1717.034, 2295.722, 2405.384, 2827.728, 2648.367, 875.615, 1641.026, 2966.603 };
	//double v[8] = { 1527.403, 1526.682, 1035.059,  891.408,  616.730, 713.137,  920.560,  639.781 };
	//double X[8] = { -13.41,     -14.29,   -76.96,   -39.37,   -60.12,  -50.43,   -44.72,   -18.09 };
	//double Y[8] = {  -6.54,      -3.40,   -13.74,    -0.38,    -4.41,  -50.31,   -23.70,     0.61 };
	//double Z[8] = {  -1.91,      -1.97,     3.62,     3.83,    12.75,   12.92,      5.1,     3.46 };
	//int count = 8;
	//Matrix<double, 3, 1> Tci;
	//Tci(0, 0) = 0.0;		
	//Tci(1, 0) = 1.22957;	
	//Tci(2, 0) = 0.0917;//
	//double yaw_iw = -108.1216*EIGEN_PI / 180.0;
	//double pitch_iw = 2.1562*EIGEN_PI / 180.0;
	//double roll_iw = 1.4324*EIGEN_PI / 180.0;//

	double u[7] = { 1393.380, 2689.992, 2495.512, 2192.429, 2162.663, 1752.160, 1168.072 };
	double v[7] = { 1640.526, 1733.708, 1329.823,  960.686,  757.620,  804.560,  240.564 };
	double X[7] = { 10.355,    9.851,   29.828,   59.028,   43.487,   29.995,   47.240 };
	double Y[7] = { 5.839,    0.386,    2.825,   12.698,    9.822,   12.228,   34.291 };
	double Z[7] = { -1.850,   -1.890,   -1.823,    5.230,    7.590,    4.950,   22.260 };
	int count = 7;

	// camera position in imu coordinate system
	Matrix<double, 3, 1> Tci;

	Tci(0, 0) = 0.0566;
	Tci(1, 0) = 1.2348;
	Tci(2, 0) = 0.088;

	// camera pose in imu coordinate system
	Matrix<double, 3, 3> Rci;
	Rci.setIdentity();
	double y_ci = 0;
	double p_ci = -90;
	double r_ci = 0;
	angle2matrix<double>(Rci, y_ci* EIGEN_PI / 180.0, p_ci* EIGEN_PI / 180.0, r_ci* EIGEN_PI / 180.0);

	// imu pose in world coordinate system	
	double yaw_iw = 74.5924*EIGEN_PI / 180.0;
	double pitch_iw = 2.4326*EIGEN_PI / 180.0;
	double roll_iw = 0.8178*EIGEN_PI / 180.0;
	Matrix<double, 3, 3>  Riw;
	angle2matrix<double>(Riw, yaw_iw, pitch_iw, roll_iw);


	// coordinates of point cloud have been translated to the center of imu
	Matrix<double, 3, 1> Tiw;
	Tiw.setZero();

	// camera rotation matrix in world coordinate system
	Matrix<double, 3, 3> Rcw = Riw*Rci;
	Matrix<double, 3, 1> Tcw = Riw*Tci + Tiw;

	double yaw_cw = 0, pitch_cw = 0, roll_cw = 0;
	matrix2angle<double>(Rcw, yaw_cw, pitch_cw, roll_cw);

	double yaw0 = yaw_cw / EIGEN_PI * 180;
	double pitch0 = pitch_cw / EIGEN_PI * 180;
	double roll0 = roll_cw / EIGEN_PI * 180;

	// variables to solve for£ºcamera position and pose in world coordinate system
	double resection_params[6] = { 0 };   // initial value 
	resection_params[0] = yaw_cw;// yaw
	resection_params[1] = pitch_cw;// pitch
	resection_params[2] = roll_cw; // roll

	resection_params[3] = Tcw(0, 0);
	resection_params[4] = Tcw(1, 0);
	resection_params[5] = Tcw(2, 0);

	double fx = 2393.6262497151; //
	double fy = 2393.7278226826; //focal length in pixel unit

	double cx = 2067.7557568795, cy = 1084.3810612941;// principle point 

													  // distortion parameters: k1,k2,p1,p2,k3,k4,k5,k6
	double distortion[8] = { -0.0701450958, 0.3824863463, -0.0001965670, -0.0004805648, -0.9935583385, 0.1023169459, 0.3095196006, -1.0635275856 };

	Problem problem;
	for (int i = 0; i < count; ++i)
	{
		ResectionResidualWithDistortion* residual =
			new ResectionResidualWithDistortion(u[i], v[i], X[i], Y[i], Z[i],
				fx, fy, cx, cy, distortion,
				Tcw(0, 0), Tcw(1, 0), Tcw(2, 0));

		problem.AddResidualBlock(new AutoDiffCostFunction<ResectionResidualWithDistortion, 2, 3>(residual), NULL, resection_params);
	}

	Solver::Options m_options;
	Solver::Summary m_summary;
	m_options.minimizer_progress_to_stdout = false;

	// solve the problem
	ceres::Solve(m_options, &problem, &m_summary);

	// output the computed results
	FILE* fp = stdout;
	fprintf(fp, "\nresection with distortion\n");
	fprintf(fp, "camera position and pose in world coordinate system\n");
	fprintf(fp, "initial value\n");
	fprintf(fp, "Yaw0 = %.3lf¡ã\tPitch0 = %.3lf¡ã\tRoll0 = %.3lf¡ã\n", yaw0, pitch0, roll0);
	fprintf(fp, "Xs0 = %.3lf\tYs0 = %.3lf\tZs0 = %.3lf\n", Tcw(0, 0), Tcw(1, 0), Tcw(2, 0));
	fprintf(fp, "optimized value\n");
	fprintf(fp, "Yaw1 = %.3lf¡ã\tPitch1 = %.3lf¡ã\tRoll1 = %.3lf¡ã\n",
		resection_params[0] / EIGEN_PI * 180,
		resection_params[1] / EIGEN_PI * 180,
		resection_params[2] / EIGEN_PI * 180);
	fprintf(fp, "Xs1 = %.3lf\tYs1 = %.3lf\tZs1 = %.3lf\n",
		resection_params[3], resection_params[4], resection_params[5]);

	fprintf(fp, "residual errors in pixel unit\n");
	double errors[2] = { 0 };
	for (size_t i = 0; i < count; i++)
	{
		ResectionResidualWithDistortion res(u[i], v[i], X[i], Y[i], Z[i], fx, fy, cx, cy, distortion, Tcw(0, 0), Tcw(1, 0), Tcw(2, 0));
		res(resection_params, errors);
		fprintf(fp, "%d\t%lf\t%lf\n", i + 1, errors[0], errors[1]);
	}

	Matrix3d Rcw2;
	angle2matrix<double>(Rcw2, resection_params[0], resection_params[1], resection_params[2]);

	Matrix3d Rci2 = Riw.transpose()*Rcw2;
	double yaw_ci = 0, pitch_ci = 0, roll_ci = 0;
	matrix2angle<double>(Rci2, yaw_ci, pitch_ci, roll_ci);

	yaw_ci = yaw_ci / EIGEN_PI * 180;
	pitch_ci = pitch_ci / EIGEN_PI * 180;
	roll_ci = roll_ci / EIGEN_PI * 180;

	fprintf(fp, "camera pose in imu coordinate system\n");
	fprintf(fp, "initial value\n");
	fprintf(fp, "Yaw0 = %.3lf¡ã\tPitch0 = %.3lf¡ã\tRoll0 = %.3lf¡ã\n", y_ci, p_ci, r_ci);
	fprintf(fp, "optimized value\n");
	fprintf(fp, "Yaw  = %.3lf¡ã\tPitch  = %.3lf¡ã\tRoll  = %.3lf¡ã\n", yaw_ci, pitch_ci, roll_ci);

}
