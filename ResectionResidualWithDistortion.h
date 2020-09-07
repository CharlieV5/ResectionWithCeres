#pragma once
#include "Eigen\Eigen"
using namespace Eigen;

template <typename T>
void angle2matrix(Matrix<T, 3, 3>& m, const T& yaw, const T& pitch, const T& roll)
{
	T sinYaw = sin(yaw);
	T cosYaw = cos(yaw);
	T sinPitch = sin(pitch);
	T cosPitch = cos(pitch);
	T sinRoll = sin(roll);
	T cosRoll = cos(roll);

	m(0, 0) =  cosYaw*cosRoll + sinYaw*sinPitch*sinRoll;    m(0, 1) = sinYaw*cosPitch;    m(0, 2) =  cosYaw*sinRoll - sinYaw*sinPitch*cosRoll;
	m(1, 0) = -sinYaw*cosRoll + cosYaw*sinPitch*sinRoll;    m(1, 1) = cosYaw*cosPitch;    m(1, 2) = -sinYaw*sinRoll - cosYaw*sinPitch*cosRoll;
	m(2, 0) = -cosPitch*sinRoll;		                    m(2, 1) = sinPitch;           m(2, 2) =  cosPitch*cosRoll;
}

template <typename T>
void matrix2angle(const Matrix<T, 3, 3>& m, T& yaw, T& pitch, T& roll)
{
	yaw   = atan2(m(0, 1), m(1, 1));
	pitch = atan2(m(2, 1), sqrt(m(2, 0)*m(2, 0) + m(2, 2)*m(2, 2)));
	//pitch = asin(m(2, 1));
	roll  = atan2(-m(2, 0), m(2, 2));
}

struct ResectionResidualWithDistortion
{
	ResectionResidualWithDistortion(double _u, double _v, 
		double _X, double _Y, double _Z, 
		double _fx, double _fy, double _cx, double _cy,
		double* _distortion, double _Xcw, double _Ycw, double _Zcw)
		: u(_u), v(_v), X(_X), Y(_Y), Z(_Z), fx(_fx), fy(_fy), cx(_cx), cy(_cy), Xcw(_Xcw), Ycw(_Ycw), Zcw(_Zcw)
	{
		for (size_t i = 0; i < 8; i++)
		{
			distortion_params[i] = _distortion[i];// k1, k2, p2, p2, k3, k4, k5, k6
		}
	}

	template <typename T>
	bool operator () (const T * const resection_params, T* residual) const
	{
		// elements of exterior orientation
		// rotational angles 
		T Yaw   = resection_params[0];
		T Pitch = resection_params[1];
		T Roll  = resection_params[2];
		//// translation parameters
		//T Xcw = resection_params[3];
		//T Ycw = resection_params[4];
		//T Zcw = resection_params[5];

		Matrix<T, 3, 3> Rcw;
		angle2matrix(Rcw, Yaw, Pitch, Roll);
		Matrix<T, 3, 1> Pw;
		Matrix<T, 3, 1> Pc;
		Matrix<T, 3, 1> Tcw;

		Pw(0, 0) = T(X);
		Pw(1, 0) = T(Y);
		Pw(2, 0) = T(Z);
		Tcw(0, 0) = T(Xcw);
		Tcw(1, 0) = T(Ycw);
		Tcw(2, 0) = T(Zcw);

		// Pw = Rcw*Pc + Tcw
		Pc = Rcw.transpose()*(Pw - Tcw);

		T x = Pc(0, 0) / Pc(2, 0);
		T y = Pc(1, 0) / Pc(2, 0);

		T r2 = x*x + y*y;
		T r4 = r2*r2;
		T r6 = r4*r2;

		// radial distortion, and tangential distortion
		T k1 = (T)distortion_params[0];
		T k2 = (T)distortion_params[1];
		T p1 = (T)distortion_params[2];
		T p2 = (T)distortion_params[3];
		T k3 = (T)distortion_params[4];
		T k4 = (T)distortion_params[5];
		T k5 = (T)distortion_params[6];
		T k6 = (T)distortion_params[7];

		T radial_distortion = (T(1) + k1*r2 + k2*r4 + k3*r6)/(T(1) + k4*r2 + k5*r4 + k6*r6);
		T tan_distortion_x = T(2) * p1*x*y + p2*(r2 + T(2)*x*x);
		T tan_distortion_y = T(2) * p2*x*y + p1*(r2 + T(2)*y*y);
		T x_ = x*radial_distortion + tan_distortion_x;
		T y_ = y*radial_distortion + tan_distortion_y;

		T u_ = fx*x_ + cx;
		T v_ = fy*y_ + cy;

		// residual errors
		residual[0] = u_ - T(u);
		residual[1] = v_ - T(v);

		return true;
	}

private:
	// control point in image plane coordinate system in pixel unit
	const double u;
	const double v;

	// focal length in pixel unit
	const double fx;
	const double fy;

	// principal point in pixel unit
	const double cx;
	const double cy;

	////  linear elements of exterior orientation 
	const double Xcw;
	const double Ycw;
	const double Zcw;

	double distortion_params[8];

	// control point in world coordinate system
	const double X;
	const double Y;
	const double Z;

};