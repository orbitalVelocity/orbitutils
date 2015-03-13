#include <cmath>
#include <vector>
#include "lambert.h"

std::vector<double> boundingVelocities(double grav_param, std::vector<double> r0,
	std::vector<double> rf, double dt, bool long_way)
{
	// algorithm constants
	const int MAXITER = 1000;
	const double TOLERANCE = 1E-10;

	// switch to canonical units
	// TODO

	// problem constants
	double r0_norm = sqrt(r0[0]*r0[0] + r0[1]*r0[1] + r0[2]*r0[2]);
	double rf_norm = sqrt(rf[0]*rf[0] + rf[1]*rf[1] + rf[2]*rf[2]);
	double df = acos((r0[0]*rf[0] + r0[1]*rf[1] + r0[2]*rf[2])/(r0_norm*rf_norm));
	if(long_way == 1) df = 2.0*M_PI - df;
	double k = r0_norm*rf_norm*(1.0-cos(df));
	double l = r0_norm + rf_norm;
	double m = r0_norm*rf_norm*(1.0+cos(df));
	double p_i = k/(1.0+sqrt(2.0*m));
	double p_ii = k/(1.0-sqrt(2.0*m));

	// initial estimate
	double p = (p_i + p_ii)/2.0;

	// f and g functions
	double f,g,fdot,gdot;

	// Newton-Raphson iteration
	for(int i = 0; i < MAXITER; ++i)
	{
		double t, dtdp;
		// deal with negative values of p
		if(p < 0.0) p = ((MAXITER-i)*p_i+i*p_ii)/MAXITER;
		// deal with out-of-bounds p
		if(!long_way && p < p_i) p = (i*p_i+(MAXITER-i)*p_ii)/MAXITER;
		if(long_way && p > p_ii) p = p = (i*p_i+(MAXITER-i)*p_ii)/MAXITER;
		// compute the semi-major axis using the new p
		double a = m*k*p/((2*m-l*l)*p*p + 2*k*l*p-k*k);
		// compute the values of the f and g functions
		f = 1-rf_norm/p*(1-cos(df));
		g = r0_norm*rf_norm*sin(df)/sqrt(p);
		fdot = sqrt(1.0/p)*tan(df/2.0)*((1.0-cos(df))/p-1/r0_norm-1/rf_norm);
		gdot = (1.0+fdot*g)/f;
		// compute the time-of-flight and its derivative with respect to p
		if(a > 0.0)
		{
			double cosE = 1.0-r0_norm/a*(1.0-f);
			double sinE = -r0_norm*rf_norm*fdot/sqrt(a);
			double E = atan2(sinE,cosE);
			if(E < 0.0) E += 2.0*M_PI;
			t = g + sqrt(a*a*a)*(E-sin(E));
			dtdp = g/(2*p)-3/2*a*(t-g)*((k*k+(2*m-l*l)*p*p)/(m*k*p*p))+sqrt(a*a*a)*(2*k*sinE)/(p*(k-l*p));
		}
		else
		{
			double coshF = 1.0-r0_norm/a*(1-f);
			double F = acosh(coshF);
			t = g + sqrt(-a*a*a)*(sinh(F)-F);
			dtdp = -g/(2*p)-3/2*a*(t-g)*((k*k+(2*m-l*l)*p*p)/(m*k*p*p))-sqrt(-a*a*a)*(2*k*sinh(F))/(p*(k-l*p));
		}
		// iterate
		p = p + (dt-t)/dtdp;
		if(std::abs(dt-t) < TOLERANCE) break;
	}

	// compute and return the initial and final velocities
	std::vector<double> v(6);
	v[0] = (rf[0]-f*r0[0])/g;
	v[1] = (rf[1]-f*r0[1])/g;
	v[2] = (rf[2]-f*r0[2])/g;
	v[3] = fdot*r0[0]+gdot*v[0];
	v[4] = fdot*r0[1]+gdot*v[1];
	v[5] = fdot*r0[2]+gdot*v[2];
	return v;
}
