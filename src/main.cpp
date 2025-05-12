#include <iostream>
#include <fstream>
#include <vector>
#include <chrono> 
#include "trust_region_mpc_solver.h"

void save_lanes_to_csv(const VehicleState& vehicle_state, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    // Write CSV header
    file << "lane_type,x,y,s\n";
    
    std::vector<std::pair<std::string, Lane>> lanes = {
        {"center", vehicle_state.centerlane},
        {"left", vehicle_state.leftlane},
        {"right", vehicle_state.rightlane}
    };

    for (const auto& [lane_type, lane] : lanes) {
        if (lane.present) {  // Only save if the lane exists
            int num_samples = 20;  // Number of points along the lane
            for (int i = 0; i < num_samples; i++) {
                double s = i * (lane.s_max / num_samples);
                double x = lane.spline_x(s);
                double y = lane.spline_y(s);

                file << lane_type << "," << x << "," << y << "," << s << "\n";
            }
        }
    }

    file.close();
    std::cout << "Lane saved to " << filename << std::endl;
}

void save_trajectories_to_csv(const VehicleState& vehicle_state, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    file << "vehicle_id,x,y,psi,s,time\n";

    for (const auto& point : vehicle_state.predicted_trajectory) {
        file << "0," << point.x << "," << point.y << "," << point.psi << "," << point.s << "," << point.t_start << "\n";
    }

    file.close();
    std::cout << "Trajectories saved to " << filename << std::endl;
}

int main() {
    
    //---------------------------------------- SMOOTH CURVE --------------------------------------------------------
    std::cerr<<"------------------------ Smooth Curve Scenario -----------------------------"<<"\n";

    // Define traffic participants (3 vehicles approaching an intersection)
    VehicleState vehicle_state_smooth_curve = 
        // x, y, v, psi, beta, a, v_target
        {0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 10.0};  
    
    int centerlane_length = 20;
    double step = 5.0; // Distance between points
    double curvature_radius = 100.0; // Larger = smoother curve
    
    std::vector<double> x_vals, y_vals, s_vals;
    
    for (int j = 0; j < centerlane_length; j++) {
        double s = j * step;
        double theta = s / curvature_radius; // Angle for arc (radians)
    
        double x = vehicle_state_smooth_curve.x + curvature_radius * sin(theta);
        double y = vehicle_state_smooth_curve.y + curvature_radius * (1 - cos(theta));
    
        s_vals.push_back(s);
        x_vals.push_back(x);
        y_vals.push_back(y);
    }

    vehicle_state_smooth_curve.centerlane.initialize_spline(x_vals, y_vals, s_vals);

    // Run the mpc_solver
    TrustRegionMPCSolver mpc_solver(vehicle_state_smooth_curve);

    auto start_time = std::chrono::high_resolution_clock::now();

    mpc_solver.run();

    auto end_time = std::chrono::high_resolution_clock::now();

    // Compute duration in milliseconds
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    std::cout << "Execution Time for run(): " << elapsed_time.count() << " ms" << std::endl;

    vehicle_state_smooth_curve = mpc_solver.vehicle_state;

    // Save trajectories to a CSV file
    save_trajectories_to_csv(vehicle_state_smooth_curve, "../trajectory_smooth_curve.csv");
    save_lanes_to_csv(vehicle_state_smooth_curve, "../lanes_smooth_curve.csv");

    //---------------------------------------- HARSH CURVE --------------------------------------------------------
    std::cerr << "------------------------ Harsh Curve Scenario -----------------------------" << "\n";

    // Define traffic participant (starting along X, turning to Y)
    VehicleState vehicle_state_harsh_curve = 
        // x, y, v, psi, beta, a, v_target
        {0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 5.0};  

    centerlane_length = 10;
    step = 5.0; // Distance between points
    curvature_radius = 10.0; // Larger = smoother curve

    std::vector<double> x_vals_harsh, y_vals_harsh, s_vals_harsh;

    for (int j = 0; j <= centerlane_length; j++) {
        double s = j * step;
        double theta = s / curvature_radius; // Angle for arc (radians)
    
        double x = vehicle_state_harsh_curve.x + curvature_radius * sin(theta);
        double y = vehicle_state_harsh_curve.y + curvature_radius * (1 - cos(theta));
    
        s_vals_harsh.push_back(s);
        x_vals_harsh.push_back(x);
        y_vals_harsh.push_back(y);
    }

    vehicle_state_harsh_curve.centerlane.initialize_spline(x_vals_harsh, y_vals_harsh, s_vals_harsh);

    // Run the solver
    TrustRegionMPCSolver mpc_solver_harsh(vehicle_state_harsh_curve);

    auto start_time_harsh = std::chrono::high_resolution_clock::now();
    mpc_solver_harsh.run();
    auto end_time_harsh = std::chrono::high_resolution_clock::now();

    auto elapsed_time_harsh = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_harsh - start_time_harsh);

    std::cout << "Execution Time for harsh run(): " << elapsed_time_harsh.count() << " ms" << std::endl;

    vehicle_state_harsh_curve = mpc_solver_harsh.vehicle_state;

    // Save results
    save_trajectories_to_csv(vehicle_state_harsh_curve, "../trajectory_harsh_curve.csv");
    save_lanes_to_csv(vehicle_state_harsh_curve, "../lanes_harsh_curve.csv");

    return 0;
}

