#include <cmath>
#include <vector>
#include "orbitalelements.h"
#include "oeconvert.h"

// TODO: Handle corner cases.

/**
 * Converts orbital elements to position.
 * @param	grav_param	The gravitational parameter of the two-body system.
 * @param	oe	The orbital elements to be converted from.
 * @returns	A vector corresponding to the position.
 */
std::vector<double> oe2r(double grav_param, OrbitalElements oe)
{
	double R11 = cos(oe.aop)*cos(oe.lan)-cos(oe.inc)*sin(oe.aop)*sin(oe.lan);
	double R12 = -cos(oe.lan)*sin(oe.aop)-cos(oe.inc)*cos(oe.aop)*sin(oe.lan);
	//double R13 = sin(oe.inc)*sin(oe.lan);
	double R21 = cos(oe.inc)*cos(oe.lan)*sin(oe.aop)+cos(oe.aop)*sin(oe.lan);
	double R22 = cos(oe.inc)*cos(oe.aop)*cos(oe.lan)-sin(oe.aop)*sin(oe.lan);
	//double R23 = -cos(oe.lan)*sin(oe.inc);
	double R31 = sin(oe.inc)*sin(oe.aop);
	double R32 = cos(oe.aop)*sin(oe.inc);
	//double R33 = cos(oe.inc);

	std::vector<double> r_pf(3);
	double r_norm = oe.sma*(1-oe.ecc*oe.ecc)/(1+oe.ecc*cos(oe.tra));
	r_pf[0] = r_norm*cos(oe.tra);
	r_pf[1] = r_norm*sin(oe.tra);
	r_pf[2] = 0.0;

	std::vector<double> r(3);
	r[0] = R11*r_pf[0] + R12*r_pf[1] /*+R13*r_pf[2]*/;
	r[1] = R21*r_pf[0] + R22*r_pf[1] /*+R23*r_pf[2]*/;
	r[2] = R31*r_pf[0] + R32*r_pf[1] /*+R33*r_pf[2]*/;

	return r;
}

/**
 * Converts orbital elements to position and velocity.
 * @param	grav_param	The gravitational parameter of the two-body system.
 * @param	oe	The orbital elements to be converted from.
 * @returns	A vector corresponding to the concatenated position and velocity.
 */
std::vector<double> oe2rv(double grav_param, OrbitalElements oe)
{
	// rotation matrix
	double R11 = cos(oe.aop)*cos(oe.lan)-cos(oe.inc)*sin(oe.aop)*sin(oe.lan);
	double R12 = -cos(oe.lan)*sin(oe.aop)-cos(oe.inc)*cos(oe.aop)*sin(oe.lan);
	//double R13 = sin(oe.inc)*sin(oe.lan);
	double R21 = cos(oe.inc)*cos(oe.lan)*sin(oe.aop)+cos(oe.aop)*sin(oe.lan);
	double R22 = cos(oe.inc)*cos(oe.aop)*cos(oe.lan)-sin(oe.aop)*sin(oe.lan);        
	//double R23 = -cos(oe.lan)*sin(oe.inc);        
	double R31 = sin(oe.inc)*sin(oe.aop);        
	double R32 = cos(oe.aop)*sin(oe.inc);        
	//double R33 = cos(oe.inc);        

	// semi-latus rectum
	double p = oe.sma*(1-oe.ecc*oe.ecc);

	// position in the perifocal frame
	std::vector<double> r_pf(3);        
	double r_norm = p/(1+oe.ecc*cos(oe.tra));        
	r_pf[0] = r_norm*cos(oe.tra);        
	r_pf[1] = r_norm*sin(oe.tra);        
	r_pf[2] = 0.0;

	// velocity in the perifocal frame
	std::vector<double> v_pf(3);
	v_pf[0] = sqrt(grav_param/p)*-sin(oe.tra);
	v_pf[1] = sqrt(grav_param/p)*(oe.ecc+cos(oe.tra));
	v_pf[2] = 0.0;

	// rotate the position and velocity into the body-fixed inertial frame
	std::vector<double> rv(6);        
	rv[0] = R11*r_pf[0] + R12*r_pf[1] /*+R13*r_pf[2]*/;
	rv[1] = R21*r_pf[0] + R22*r_pf[1] /*+R23*r_pf[2]*/;
	rv[2] = R31*r_pf[0] + R32*r_pf[1] /*+R33*r_pf[2]*/;
	rv[3] = R11*v_pf[0] + R12*v_pf[1] /*+R13*v_pf[2]*/;
	rv[4] = R21*v_pf[0] + R22*v_pf[1] /*+R23*v_pf[2]*/;
	rv[5] = R31*v_pf[0] + R32*v_pf[1] /*+R33*v_pf[2]*/;

	return rv;
}

/**
 * Converts position and velocity to orbital elements
 * @param	grav_param	The gravitational parameter of the two-body system.
 * @param	rv	The concatenated position and velocity to convert from.
 * @returns	The orbital elements.
 */
OrbitalElements rv2oe(double grav_param, std::vector<double> rv)
{
	OrbitalElements oe;

	// Semi-major Axis : Vis-viva Equation
	oe.sma = 1.0 / (
		2.0 / sqrt(rv[0]*rv[0] + rv[1]*rv[1] + rv[2]*rv[2])
		- (rv[3]*rv[3] + rv[4]*rv[4] + rv[5]*rv[5]) / grav_param
	);

	// Angular Momentum
	std::vector<double> h(3);
	h[0] = rv[1]*rv[5]-rv[2]*rv[4];
	h[1] = rv[2]*rv[3]-rv[0]*rv[5];
	h[2] = rv[0]*rv[4]-rv[1]*rv[3];

	// Norm of position
	double r = sqrt(rv[0]*rv[0] + rv[1]*rv[1] + rv[2]*rv[2]);

	// Eccentricity Vector :  e = v x h / mu - r/|r|
	std::vector<double> e(3);
	e[0] = (rv[4]*h[2]-rv[5]*h[1])/grav_param - rv[0]/r;
	e[1] = (rv[5]*h[0]-rv[3]*h[2])/grav_param - rv[1]/r;
	e[2] = (rv[3]*h[1]-rv[4]*h[0])/grav_param - rv[2]/r;

	// Eccentricity
	oe.ecc = sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);

	// Inclination
	oe.inc = acos(h[2]/sqrt(h[0]*h[0]+h[1]*h[1]+h[2]*h[2]));

	// Ascending Node Direction (In x-y plane)
	std::vector<double> n(2);
	n[0] = -h[1]/sqrt(h[0]*h[0]+h[1]*h[1]);
	n[1] = h[0]/sqrt(h[0]*h[0]+h[1]*h[1]);
	double n_norm = sqrt(n[0]*n[0]+n[1]*n[1]);

	// Longitude of the Ascending Node
	oe.lan = acos(n[0])/n_norm;
	if(n[1] < 0.0) oe.lan = 2*M_PI - oe.lan;

	// Argument of Periapsis
	oe.aop = acos((n[0]*e[0]+n[1]*e[1])/(n_norm*oe.ecc));
	if(e[2] < 0.0) oe.aop = 2*M_PI - oe.aop;

	// True Anomaly
	oe.tra = acos((rv[0]*e[0]+rv[1]*e[1]+rv[2]*e[2])/(r*oe.ecc));
	if(rv[0]*rv[3]+rv[1]*rv[4]+rv[2]*rv[5] < 0.0) oe.tra = 2*M_PI - oe.tra;

	return oe;
}
