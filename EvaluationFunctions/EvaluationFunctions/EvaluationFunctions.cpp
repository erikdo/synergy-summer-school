// EvaluationFunctions.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include <algorithm>

#define DLLEXPORT extern "C" __declspec(dllexport)

DLLEXPORT double* welded_beam_design(double* x)
{
	// problem encoding
	double b = 0.125 + 4.875 * x[0];
	double h = 0.125 + (b - 0.125) * x[1];   // encoding of g3
	double l = 0.1 + 9.9 * x[2];
	double t = 0.1 + 9.9 * x[3];

	// auxiliary terms
	double sigma = 504000.0 / (t*t * b);
	double tau1 = 6000.0 / (std::sqrt(2.0) * h * l);
	double tau2 = 6000.0 * (14.0 + 0.5*l) * std::sqrt(0.25 * (l*l + (h + t)*(h + t))) / (1.414 * h*l * (l*l / 12.0 + 0.25*(h + t)*(h + t)));
	double tau = std::sqrt(tau1*tau1 + tau2*tau2 + l*tau1*tau2 / std::sqrt(0.25 * (l*l + (h + t)*(h + t))));
	double P_c = 64746.022 * (1.0 - 0.0282346*t) * t * b*b*b;

	// fitness functions
	double f1 = 1.10471 * h*h * l + 0.04811 * t * b * (14.0 + l);
	double f2 = 2.1952 / (t*t*t * b);

	// scaling roughly to range [0, 1]
	f1 /= 40.0;
	f2 /= 0.006;

	// constraints
	double g1 = tau - 13600.0;
	double g2 = sigma - 30000.0;
	double g3 = 6000.0 - P_c;

	// penalty for constraint violations
	double penalty = ((std::max)(1e-6 * g1, 0.0) + (std::max)(1e-6 * g2, 0.0) + (std::max)(1e-7 * g3, 0.0));
	if (penalty > 0.0) penalty += 0.5;
	f1 += penalty;
	f2 += penalty;

	double *ret = (double *)malloc(2 * sizeof(double));
	ret[0] = std::tanh(f1);
	ret[1] = std::tanh(f2);
	return ret;
}
