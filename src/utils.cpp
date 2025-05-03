#include "utils.h"

// Function to compute the Euclidean distance between two points (x1, y1) and (x2, y2)
double compute_distance(double x1, double y1, double x2, double y2) {
    return std::sqrt(square(x2 - x1) + square(y2 - y1)); // Euclidean distance formula
}

// Function to compute the square of a number
double square(double value) {
    return value * value;
}

// Function to compute the dot product of two vectors (x1, y1) and (x2, y2)
double dot_product(double x1, double y1, double x2, double y2) {
    return x1 * x2 + y1 * y2; // Dot product formula
}