#include <iostream>
#include <fstream>
#include <vector>
#include <chrono> 
#include "dynamic_game_planner.h"

void save_lanes_to_csv(const std::vector<VehicleState>& traffic, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    // Write CSV header
    file << "lane_type,x,y,s\n";

    for (const auto& vehicle : traffic) {
        std::vector<std::pair<std::string, Lane>> lanes = {
            {"center", vehicle.centerlane},
            {"left", vehicle.leftlane},
            {"right", vehicle.rightlane}
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
    }

    file.close();
    std::cout << "Lanes saved to " << filename << std::endl;
}

void save_trajectories_to_csv(const std::vector<VehicleState>& traffic, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    file << "vehicle_id,x,y,psi,s,time\n";

    for (size_t i = 0; i < traffic.size(); i++) {
        for (const auto& point : traffic[i].predicted_trajectory) {
            file << i << "," << point.x << "," << point.y << "," << point.psi << "," << point.s << "," << point.t_start << "\n";
        }
    }

    file.close();
    std::cout << "Trajectories saved to " << filename << std::endl;
}

int main() {
    
    //---------------------------------------- INTERSECTION --------------------------------------------------------
    std::cerr<<"------------------------ Intersection Scenario -----------------------------"<<"\n";

    // Define traffic participants (3 vehicles approaching an intersection)
    TrafficParticipants traffic_intersection = {
        // x, y, v, psi, beta, a, v_target
        {0.0, 0.0, 10.0, 0.0, 0.0, 1.0, 10.0},  // Vehicle 1 (moving along X)
        {10.0, -10.0, 10.0, 1.57, 0.0, 1.0, 10.0}, // Vehicle 2 (coming from bottom Y)
        {-10.0, 10.0, 10.0, -1.57, 0.0, 1.0, 10.0} // Vehicle 3 (coming from top Y)
    };

    // Generate center lanes for each vehicle
    int centerlane_length = 50;
    
    for (size_t i = 0; i < traffic_intersection.size(); i++) {
        std::vector<double> x_vals, y_vals, s_vals;

        for (int j = 0; j < centerlane_length; j++) {
            if (i == 0) { 
                x_vals.push_back(traffic_intersection[i].x + j * 5.0); // Move forward in X
                y_vals.push_back(traffic_intersection[i].y);
            } else if (i == 1) { 
                x_vals.push_back(traffic_intersection[i].x);
                y_vals.push_back(traffic_intersection[i].y + j * 5.0); // Move forward in Y
            } else if (i == 2) {
                x_vals.push_back(traffic_intersection[i].x);
                y_vals.push_back(traffic_intersection[i].y - j * 5.0); // Move downward in Y
            }
            s_vals.push_back(j * 5.0);
        }

        traffic_intersection[i].centerlane.initialize_spline(x_vals, y_vals, s_vals);
    }

    // Run the planner
    DynamicGamePlanner planner;

    auto start_time = std::chrono::high_resolution_clock::now();

    planner.run(traffic_intersection);

    auto end_time = std::chrono::high_resolution_clock::now();

    // Compute duration in milliseconds
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    std::cout << "Execution Time for run(): " << elapsed_time.count() << " ms" << std::endl;

    traffic_intersection = planner.traffic;

    // Save trajectories to a CSV file
    save_trajectories_to_csv(traffic_intersection, "../trajectories_intersection.csv");
    save_lanes_to_csv(traffic_intersection, "../lanes_intersection.csv");

    //-------------------------------------------- MERGING SCENARIO --------------------------------------------------
    std::cerr<<"------------------------ Merging Scenario -----------------------------"<<"\n";

    // Define traffic participants (2 vehicles with center lanes that unify)
    std::vector<VehicleState> traffic_merging = {
        {0.0, 0.0, 2.0, 0.0, 0.0, 1.0, 10.0},  // Vehicle 1 (moving along X)
        {10.0, -5.0, 10.0, 0.1, 0.0, 1.0, 10.0} // Vehicle 2 (coming from bottom Y)
    };

    // Length of the center lane (same for both vehicles)
    centerlane_length = 50;

    // Angle between the two center lanes in radians (10 degrees)
    double angle = 10.0 * M_PI / 180.0;  // 10 degrees to radians

    for (size_t i = 0; i < traffic_merging.size(); i++) {
        std::vector<double> x_vals, y_vals, s_vals;

        for (int j = 0; j < centerlane_length; j++) {
            if (i == 0) { 
                // Vehicle 1 moves straight along X-axis
                x_vals.push_back(traffic_merging[i].x + j * 5.0);
                y_vals.push_back(traffic_merging[i].y);
            } else if (i == 1) { 
                // Vehicle 2 moves from bottom, and its path should meet vehicle 1 at an angle of 10 degrees
                double progress = j / static_cast<double>(centerlane_length);
                
                // Using trigonometry to calculate the position of vehicle 2
                double x_offset = progress * 10.0 * cos(angle);
                double y_offset = progress * 10.0 * sin(angle);
                
                // Set the new x, y values
                x_vals.push_back(traffic_merging[i].x + x_offset);  // Moves leftward
                y_vals.push_back(traffic_merging[i].y + y_offset);  // Moves upward
            }
            
            s_vals.push_back(j * 5.0);  // Progression along the lane
        }

        // Initialize the lane with spline interpolation
        traffic_merging[i].centerlane.initialize_spline(x_vals, y_vals, s_vals);
    }

    // Run the planner
    start_time = std::chrono::high_resolution_clock::now();

    planner.run(traffic_merging);

    end_time = std::chrono::high_resolution_clock::now();

    // Compute duration in milliseconds
    elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    std::cout << "Execution Time for run(): " << elapsed_time.count() << " ms" << std::endl;

    traffic_merging = planner.traffic;

    // Save trajectories to a CSV file
    save_trajectories_to_csv(traffic_merging, "../trajectories_merging.csv");
    save_lanes_to_csv(traffic_merging, "../lanes_merging.csv");

    //---------------------------------------- OVERTAKING --------------------------------------------------------
    std::cerr<<"------------------------ Overtaking Scenario -----------------------------"<<"\n";

    // Define traffic participants (2 vehicles in column)
    TrafficParticipants traffic_overtaking = {
        // x, y, v, psi, beta, a
        {0.0, 0.0, 2.0, 0.0, 0.0, 1.0, 5.0},  // Vehicle 1 (moving along X)
        {-12.0, 0.0, 10.0, 0.0, 0.0, 1.0, 10.0}, // Vehicle 2 (moving along X)
    };

    // Generate center lanes for each vehicle
    centerlane_length = 50;
    
    for (size_t i = 0; i < traffic_overtaking.size(); i++) {
        std::vector<double> x_vals, y_vals, s_vals;

        for (int j = 0; j < centerlane_length; j++) {
            if (i == 0) { 
                x_vals.push_back(traffic_overtaking[i].x + j * 5.0); // Move forward in X
                y_vals.push_back(traffic_overtaking[i].y);
            } else if (i == 1) { 
                x_vals.push_back(traffic_overtaking[i].x + j * 5.0); // Move forward in X
                y_vals.push_back(traffic_overtaking[i].y);
            } 
            s_vals.push_back(j * 5.0);
        }

        traffic_overtaking[i].centerlane.initialize_spline(x_vals, y_vals, s_vals);

        x_vals.clear();
        y_vals.clear();
        s_vals.clear();

        for (int j = 0; j < centerlane_length; j++) {
            if (i == 0) { 
                x_vals.push_back(traffic_overtaking[i].x + j * 5.0); // Move forward in X
                y_vals.push_back(traffic_overtaking[i].y + 3.0);
            } else if (i == 1) { 
                x_vals.push_back(traffic_overtaking[i].x + j * 5.0); // Move forward in X
                y_vals.push_back(traffic_overtaking[i].y + 3.0);
            } 
            s_vals.push_back(j * 5.0);
        }

        traffic_overtaking[i].leftlane.initialize_spline(x_vals, y_vals, s_vals);
    }

    // Run the planner
    start_time = std::chrono::high_resolution_clock::now();

    planner.run(traffic_overtaking);

    end_time = std::chrono::high_resolution_clock::now();

    // Compute duration in milliseconds
    elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    std::cout << "Execution Time for run(): " << elapsed_time.count() << " ms" << std::endl;

    traffic_overtaking = planner.traffic;

    // Save trajectories to a CSV file
    save_trajectories_to_csv(traffic_overtaking, "../trajectories_overtaking.csv");
    save_lanes_to_csv(traffic_overtaking, "../lanes_overtaking.csv");

    return 0;
}

