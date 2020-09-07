#include "ResectionResidual.h"
#include "Resection.h"
#include "ceres\ceres.h"
using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

Resection::Resection()
{
}


Resection::~Resection()
{
}

void Resection::doResection()
{
	double x[4] = { -86.15, -53.4, -14.78, 10.46 };
	double y[4] = { -68.99, 82.21, -76.63, 64.43 };
	double X[4] = { 36589.41, 37631.08, 39100.97, 40426.54 };
	double Y[4] = { 25273.32, 31324.51, 24934.98, 30319.81 };
	double Z[4] = { 2195.17, 728.69, 2386.5, 757.31 };

	int count = 4;
	double resection_params[6] = { 0 };   // initial value
	for (int i = 0; i < count; i++)
	{
		resection_params[0] += X[i]; //X
		resection_params[1] += Y[i]; //Y
	}

	resection_params[0] /= count; // set the initial value of X 
	resection_params[1] /= count; // set the initial value of Y

	double f = 153.24 / 1000.0; //focal length, unit: meter
	double m = 5.0e4;// map scale is 1£º50000
	resection_params[2] = m * f;// estimate the height

	Problem problem;
	for (int i = 0; i < count; ++i)
	{
		ResectionResidual* residual = new ResectionResidual(x[i] / 1000.0, y[i] / 1000, f, X[i], Y[i], Z[i]);
		problem.AddResidualBlock(new AutoDiffCostFunction<ResectionResidual, 2, 6>(residual), NULL, resection_params);
	}

	Solver::Options m_options;
	Solver::Summary m_summary;
	m_options.max_num_iterations = 25;
	m_options.linear_solver_type = ceres::DENSE_QR;
	m_options.minimizer_progress_to_stdout = false;

	// solve the problem
	Solve(m_options, &problem, &m_summary);

	// output the computed results
	fprintf(stdout, "resection\n");
	fprintf(stdout, "Xcw=%.3lfm Ycw=%.3lfm Zcw=%.3lfm\n", resection_params[0], resection_params[1], resection_params[2]);
	fprintf(stdout, "phi=%.3lf¡ã omega=%.3lf¡ã kappa=%.3lf¡ã\n",
		resection_params[3] / EIGEN_PI * 180, resection_params[4] / EIGEN_PI * 180, resection_params[5] / EIGEN_PI * 180);

}

