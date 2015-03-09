#ifndef MATHUTILS
#define MATHUTILS

#include <vector>

/**
 * Computes the dot product of two length-3 vectors.
 * @param	vec1	One of the two vectors.
 * @param	vec2	The other vector.
 * @returns	The dot product of the two vectors.
 */
double dot3(std::vector<double> vec1, std::vector<double> vec2);

/**
 * Computes the cross product of two length-3 vectors.
 * @param	vec1	The first vector in the cross product.
 * @param	vec2	The second vector in the cross product.
 * @returns	The cross product of the two vectors.
 */
std::vector<double> cross3(std::vector<double> vec1, std::vector<double>);

/**
 * Computes the 2-norm of a length-3 vector.
 * @param	vec	The vector to take the norm of.
 * @returns	The 2-norm of the vector.
 */
double norm3(std::vector<double> vec);

/**
 * Computes a normalized length-3 vector.
 * @param	vec	The vector to be normalized.
 * @returns	The normalized vector.
 */
std::vector<double> normalize3(std:vector<double> vec);

/**
 * Rotates a length-3 vector counter-clockwise about the x-axis.
 * @param	angle	The angle to rotate by.
 * @param	vec	The vector to rotate.
 * @returns	The rotated vector.
 */
std::vector<double> rotx3(double angle, std::vector<double> vec);

/**
 * Rotates a length-3 vector counter-clockwise about the y-axis.
 * @param	angle	The angle to rotate by.
 * @param	vec	The vector to rotate.
 * @returns	The rotated vector.
 */
std::vector<double> roty3(double angle, std::vector<double> vec);

/**
 * Rotates a length-3 vector counter-clockwise about the z-axis.
 * @param	angle	The angle to rotate by.
 * @param	vec	The vector to rotate.
 * @returns	The rotated vector.
 */
std::vector<double> rotz3(double angle, std::vector<double> vec);

#endif
