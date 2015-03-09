#ifndef KEPLER
#define KEPLER

/**
 * Computes the delta-time until a particular true-anomaly is reached.
 * @param	grav_param	The gravitational parameter of the two-body system.
 * @param	oe	The orbital elements of the body of interest.
 * @param	true_anomaly	The true anomaly at the final time.
 * @returns	The amount of time until the true anomaly has been reached.
 */
double timeUntilAnomaly (double grav_param, OrbitalElements oe, double true_anomaly);

/**
 * Computes the true anomaly after a delta-time has passed.
 * @param	grav_param	The gravitational parameter of the two-body system.
 * @param	oe	The orbital elements of the body of interest.
 * @param	delta_time	The elapsed time.
 * @returns	The true anomaly at the final time.
 */
double anomalyAfterTime (double grav_param, OrbitalElements oe, double delta_time);

#endif
