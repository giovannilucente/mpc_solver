#ifndef UTILS_H
#define UTILS_H

#include <cmath> // for standard math functions

// Function to compute the Euclidean distance between two points (x1, y1) and (x2, y2)
double compute_distance(double x1, double y1, double x2, double y2);

// Function to compute the square of a number
double square(double value);

// Function to compute the dot product of two vectors (x1, y1) and (x2, y2)
double dot_product(double x1, double y1, double x2, double y2);

#endif // UTILS_H