#include <cmath>
#include "orbitalelements.h"
#include "kepler.h"

/**
 * Computes the delta-time until a particular true-anomaly is reached.
 * @param       grav_param      The gravitational parameter of the two-body system.
 * @param       oe      The orbital elements of the body of interest.
 * @param       true_anomaly    The true anomaly at the final time.
 * @returns     The amount of time until the true anomaly has been reached.
 */
double timeUntilAnomaly (double grav_param, OrbitalElements oe, double true_anomaly)
{
	if(oe.sma > 0.0)
	{
		if(true_anomaly < oe.tra) true_anomaly += 2*M_PI;

		double ecc_anomaly_i = 2*atan(sqrt((1-oe.ecc)/(1+oe.ecc))*tan(oe.tra/2));
		double ecc_anomaly_f = 2*atan(sqrt((1-oe.ecc)/(1+oe.ecc))*tan(true_anomaly/2));

		double mean_anomaly_i = ecc_anomaly_i-oe.ecc*sin(ecc_anomaly_i);
		double mean_anomaly_f = ecc_anomaly_f-oe.ecc*sin(ecc_anomaly_f);

		return sqrt(pow(oe.sma,3.0)/grav_param)*(mean_anomaly_f-mean_anomaly_i);
	}
	else
	{
		double hyp_anomaly_i = 2*atanh(sqrt((oe.ecc-1)/(1+oe.ecc))*tan(oe.tra/2));
		double hyp_anomaly_f = 2*atanh(sqrt((oe.ecc-1)/(1+oe.ecc))*tan(true_anomaly/2));

		double mean_anomaly_i = -hyp_anomaly_i+oe.ecc*sinh(hyp_anomaly_i);
		double mean_anomaly_f = -hyp_anomaly_f+oe.ecc*sinh(hyp_anomaly_f);

		return sqrt(pow(-oe.sma,3.0)/grav_param)*(mean_anomaly_f-mean_anomaly_i);
	}
}

/**
 * Computes the true anomaly after a delta-time has passed.
 * @param       grav_param      The gravitational parameter of the two-body system.
 * @param       oe      The orbital elements of the body of interest.
 * @param       delta_time      The elapsed time.
 * @returns     The true anomaly at the final time.
 */
double anomalyAfterTime (double grav_param, OrbitalElements oe, double delta_time)
{
	if(oe.sma > 0.0)
	{
		double ecc_anomaly_i = 2*atan(sqrt((1-oe.ecc)/(1+oe.ecc))*tan(oe.tra/2));
		double mean_anomaly_f = ecc_anomaly_i-oe.ecc*sin(ecc_anomaly_i)+sqrt(grav_param/pow(oe.sma,3.0))*delta_time;

		// perform Newton-Raphson iteration to determine the final eccentric anomaly
		double ecc_anomaly_f = mean_anomaly_f;
		double error = mean_anomaly_f-ecc_anomaly_f+oe.ecc*sin(ecc_anomaly_f);
		while (error > 1E-10)
		{
			ecc_anomaly_f = ecc_anomaly_f-error/(oe.ecc*cos(ecc_anomaly_f)-1);
			error = mean_anomaly_f-ecc_anomaly_f+oe.ecc*sin(ecc_anomaly_f);
		}

		return 2*atan(sqrt((1+oe.ecc)/(1-oe.ecc))*tan(ecc_anomaly_f/2));
	}
	else
	{
		double hyp_anomaly_i = 2*atanh(sqrt((oe.ecc-1)/(1+oe.ecc))*tan(oe.tra/2));    
		double mean_anomaly_f = -hyp_anomaly_i+oe.ecc*sinh(hyp_anomaly_i)+sqrt(grav_param/pow(-oe.sma,3.0))*delta_time;

		// perform Newton-Raphson iteration to determine the final eccentric anomaly
		double hyp_anomaly_f = mean_anomaly_f;
		double error = mean_anomaly_f+hyp_anomaly_f-oe.ecc*sinh(hyp_anomaly_f);
		while (error > 1E-10)
		{
			hyp_anomaly_f = hyp_anomaly_f-error/(oe.ecc*cosh(hyp_anomaly_f)-1);
			error = mean_anomaly_f+hyp_anomaly_f-oe.ecc*sinh(hyp_anomaly_f);
		}

		return 2*atan(sqrt((1+oe.ecc)/(oe.ecc-1))*tanh(hyp_anomaly_f/2));
	}
}
