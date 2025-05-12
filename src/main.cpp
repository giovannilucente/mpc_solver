#include <iostream>
#include <fstream>
#include <vector>
#include <chrono> 
#include <random>
#include "trust_region_mpc_solver.h"


/** Dyanamic step */
TrajectoryPoint dynamic_step(TrajectoryPoint state, Input control)
{
    TrajectoryPoint next_state;
    double dt = state.t_end - state.t_start;
    
    /* Car dynamics parameters */
    double mass = 1500.0;  // kg
    double L = 2.5;
    double Lr = L / 2.0;
    double Lf = L - Lr;
    double drag_coeff = 0.3;
    double frontal_area = 2.2;
    double air_density = 1.225;
    double rolling_resistance_coeff = 0.015;
    double max_engine_force = 4000.0;
    double max_brake_force = 6000.0;
    double g = 9.81;
    double max_steer = 0.4;
    double Iz = 2250.0;

    // Example Pacejka parameters (for illustrative purposes)
    double B_f = 10;  // Stiffness factor for front tire
    double C_f = 1.9;  // Shape factor for front tire
    double D_f = 1500;  // Peak lateral force for front tire
    double E_f = 0.97;  // Curvature factor for front tire

    double B_r = 12;  // Stiffness factor for rear tire
    double C_r = 2.0;  // Shape factor for rear tire
    double D_r = 1600;  // Peak lateral force for rear tire
    double E_r = 0.98;  // Curvature factor for rear tire

    double throttle_val = 0.0;
    double brake_val = 0.0;
    
    double u = control.T;            // signed command
    double slope = 5.0;               // slope of the switch; bigger = sharper
    
    double sigmoid  = 0.5 * (1.0 + std::tanh(slope*u));   // ∈ (0,1)

    throttle_val =  sigmoid *  u;          // positive half
    brake_val    = (1.0 - sigmoid) * ( - u);   // negative half, always ≥ 0
    
    // -------- Vehicle longitudinal dynamics --------
    double engine_force = throttle_val * max_engine_force;
    double brake_force = brake_val * max_brake_force;
    double drag_force = 0.5 * air_density * drag_coeff * frontal_area * state.v * state.v;
    //double rolling_force = rolling_resistance_coeff * mass * g;
    double rolling_force = 0.0;

    double total_force = engine_force - brake_force - drag_force - rolling_force;
    
    // -------- Vehicle lateral/slip dynamics --------
    double yaw_rate = state.omega;
    double steering_angle = control.delta;
    double beta = steering_angle * Lr / L;

    if (state.v > 1.0)
    {
        double alpha_f = steering_angle - std::atan2((yaw_rate * Lf + state.v * std::sin(beta)), (state.v * std::cos(beta) + 1e-6));
        double alpha_r = - std::atan2((yaw_rate * Lr - state.v * std::sin(beta)), (state.v * std::cos(beta) + 1e-6));

        double Fyf = - D_f * std::sin(C_f * std::atan(B_f * alpha_f - E_f * (B_f * alpha_f - std::atan(B_f * alpha_f))));
        double Fyr = - D_r * std::sin(C_r * std::atan(B_r * alpha_r - E_r * (B_r * alpha_r - std::atan(B_r * alpha_r))));

        double yaw_acc = (Lf * Fyf - Lr * Fyr) / Iz;
        yaw_rate = yaw_rate + yaw_acc * dt;

        double lat_acceleration = (Fyf + Fyr) / mass;
        double v_lat = state.v * std::sin(beta) + lat_acceleration * dt;

        beta = std::atan2(v_lat, state.v + 1e-6);
    }else{
        yaw_rate = state.v * std::tan(steering_angle) * cos((Lr / L) * steering_angle)/ L;
    }

    // Derivatives computation:
    next_state.x = state.x + dt * state.v * cos(state.psi + beta);
    next_state.y = state.y + dt * state.v * sin(state.psi + beta);
    next_state.v = state.v + dt * total_force / mass;
    next_state.psi = state.psi + dt * yaw_rate;
    next_state.omega = yaw_rate;
    next_state.a = total_force / mass;
    next_state.s = state.s + dt * state.v;
    next_state.t_start = state.t_end;
    next_state.t_end = next_state.t_start + dt;

    return next_state;
}

void drive_trajectory(VehicleState& vehicle)
{
    int N = vehicle.predicted_control.size();

    static std::default_random_engine generator(42);  
    static std::normal_distribution<double> noise_dist(0.0, 0.02);  

    TrajectoryPoint state;
    TrajectoryPoint next_state;
    state.x = vehicle.x;
    state.y = vehicle.y;
    state.v = vehicle.v;
    state.psi = vehicle.psi;
    state.omega = vehicle.omega;
    state.a = vehicle.a;
    state.s = 0.0;
    state.t_start = 0.0;
    state.t_end = vehicle.predicted_trajectory[0].t_end;

    vehicle.predicted_trajectory.clear();

    for (int j = 0; j < N; j++){

        // Input control:
        Input control = vehicle.predicted_control[j];
        control.T = control.T + noise_dist(generator);
        control.delta = control.delta + noise_dist(generator);
        
        // Dynamic step: 
        next_state = dynamic_step(state, control);

        // Saving the new state and updating the state: 
        state = next_state;
        vehicle.predicted_trajectory.push_back(next_state);

    }
}

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
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0};  
    
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

    // Save planned trajectories to a CSV file
    save_trajectories_to_csv(vehicle_state_smooth_curve, "../trajectory_smooth_curve.csv");
    save_lanes_to_csv(vehicle_state_smooth_curve, "../lanes_smooth_curve.csv");

    // Drive trajectory
    drive_trajectory(vehicle_state_smooth_curve);

    // Save driven trajectory
    save_trajectories_to_csv(vehicle_state_smooth_curve, "../driven_trajectory_smooth_curve.csv");

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

    // Drive trajectory
    drive_trajectory(vehicle_state_harsh_curve);

    // Save driven trajectory
    save_trajectories_to_csv(vehicle_state_harsh_curve, "../driven_trajectory_harsh_curve.csv");

    return 0;
}

