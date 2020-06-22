#pragma once
struct ResectionResidual
{
	ResectionResidual(double _x, double _y, double _f, double _X, double _Y, double _Z)
		: x(_x), y(_y), f(_f), X(_X), Y(_Y), Z(_Z) {}

	template <typename T>
	bool operator () (const T * const resection_params, T* residual) const
	{
		// elements of exterior orientation
		// translation parameters
		T Xs = resection_params[0];
		T Ys = resection_params[1];
		T Zs = resection_params[2];
		// rotational angles 
		T phi = resection_params[3];
		T omg = resection_params[4];
		T kap = resection_params[5];

		// rotational matrix 
		T a1 = cos(phi)*cos(kap) - sin(phi)*sin(omg)*sin(kap);
		T a2 = -cos(phi)*sin(kap) - sin(phi)*sin(omg)*cos(kap);
		T a3 = -sin(phi)*cos(omg);
		T b1 = cos(omg)*sin(kap);
		T b2 = cos(omg)*cos(kap);
		T b3 = -sin(omg);
		T c1 = sin(phi)*cos(kap) + cos(phi)*sin(omg)*sin(kap);
		T c2 = -sin(phi)*sin(kap) + cos(phi)*sin(omg)*cos(kap);
		T c3 = cos(phi)*cos(omg);

		// residual errors
		T X_ = T(a1*(X - Xs) + b1*(Y - Ys) + c1*(Z - Zs));
		T Y_ = T(a2*(X - Xs) + b2*(Y - Ys) + c2*(Z - Zs));
		T Z_ = T(a3*(X - Xs) + b3*(Y - Ys) + c3*(Z - Zs));
		residual[0] = T(x) + T(f) * X_ / Z_;
		residual[1] = T(y) + T(f) * Y_ / Z_;

		return true;
	}

private:
	// control point in image plane coordinate system
	const double x;
	const double y;

	// focal length
	const double f;

	// control point in ground coordinate system
	const double X;
	const double Y;
	const double Z;

};
