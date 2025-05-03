#include "vehicle_state.h"

void Lane::initialize_spline(const std::vector<double>& x, 
                            const std::vector<double>& y, 
                            const std::vector<double>& s) {
    if (s.size() > 2) {
        spline_x.set_points(s, x);  // Set x-coordinate spline
        spline_y.set_points(s, y);  // Set y-coordinate spline
        s_max = s.back();           // Store the last progress value
        present = true;
    } else {
        present = false;  // Not enough points to create a spline
    }
}

/** computes the heading on the spline x(s) and y(s) at parameter s*/
double Lane::compute_heading(double s)
{
    double psi;
    double dx = spline_x.deriv(1, s);
    double dy = spline_y.deriv(1, s);
    psi = atan2(dy, dx);
    if(psi < 0.0) {psi += 2*M_PI;}
    return psi;
}

/** computes the curvature on the spline x(t) and y(t) at time t*/
double Lane::compute_curvature(double s)
{
    double k;
    double dx = spline_x.deriv(1, s);
    double dy = spline_y.deriv(1, s);
    double ddx = spline_x.deriv(2, s);
    double ddy = spline_y.deriv(2, s);
    k = (ddy * dx - ddx * dy) / sqrt((dx * dx + dy * dy) * (dx * dx + dy * dy) * (dx * dx + dy * dy));
    return k;
}