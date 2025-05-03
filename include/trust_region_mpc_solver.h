#ifndef DYNAMIC_GAME_PLANNER_H
#define DYNAMIC_GAME_PLANNER_H

#include <vector>
#include <eigen3/Eigen/Dense>
#include <thread>
#include <iomanip>
#include <mutex>
#include "vehicle_state.h"
#include "utils.h"  // Utility functions

class TrustRegionMPCSolver {

private:
    const double tau = 2.0;
    const double k = 10.0;
    const double r_safe = 2.5;
    const double r_lane = 3.5;
    const double eps = 1e-6;
    const double length = 5.0;
    const double cg_ratio = 0.5;
    const double pi = 3.1415;
    const double v_max = 10.0;

    constexpr static const int N = 20;                                  /** number of integration nodes */
    constexpr static const int nX = 6;                                  /** <X, Y, V, PSI, S, L> */
    constexpr static const int nU = 2;                                  /** <d, F> */
    constexpr static const int N_interpolation = 60;
    constexpr static const double dt_interpolation = 0.1;

public:
    static const int nx = nX * (N + 1);                                 /** size of the state trajectory X */
    static const int nu = nU * (N + 1);                                 /** size of the input trajectory U */
    int nC;                                                             /** total number of inequality constraints */
    double dt = 0.3;                                                    /** integration time step */
    double d_up = 0.7;                                                  /** upper bound steering angle */
    double d_low = -0.7;                                                /** lower bound steering angle */
    double F_up = 2.0;                                                  /** upper bound force */
    double F_low = -3.0;                                                /** lower bound force */

    // Parameters:
    double qf = 1e-2;                                                   /** penalty for the final error in the lagrangian */
    double gamma = 1.3;                                                 /** increasing factor of the penalty weight */
    double rho = 1e-3;                                                  /** penalty weight */ 
    double weight_target_speed = 1e0;                                   /** weight for the maximum speed in the lagrangian */
    double weight_center_lane = 1e-1;                                   /** weight for the center lane in the lagrangian */
    double weight_heading = 1e2;                                        /** weight for the heading in the lagrangian */
    double weight_input = 0.0;                                          /** weight for the input in the lagrangian */
    
    std::vector<double> U_old;                                          /** solution in the previous iteration*/
    Eigen::MatrixXd ul;                                                 /** controls lower bound*/
    Eigen::MatrixXd uu;                                                 /** controls upper bound*/
    Eigen::MatrixXd time;                                               /** time vector */
    Eigen::MatrixXd lagrangian_multipliers;                             /** lagrangian multipliers*/
    
    enum STATES {x, y, v, psi, s, l};
    enum INPUTS {d, F};

    VehicleState vehicle_state;

    TrustRegionMPCSolver();  // Constructor
    ~TrustRegionMPCSolver(); // Destructor

    void run( VehicleState& vehicle_state );                                        /** Main method to execute the planner */
    void setup();                                                                   /** Setup function */
    void initial_guess( double* X, double* U );                                     /** Set the initial guess */
    void trust_region_solver( double* U );                                          /** solver of the mpc based on trust region */
    void integrate( double* X, const double* U );                                   /** Integration function */
    void dynamic_step( double* d_state, const double* state, const double* ref_state, 
                    const double* control );                                        /** Dynamic step function */
    void hessian_SR1_update( Eigen::MatrixXd & H_, const Eigen::MatrixXd & s_,            
                     const Eigen::MatrixXd & y_, const double r_ );                /** SR1 Hessian matrix update*/
    void increasing_schedule();                                                    /** function to increase rho = rho * gamma */
    void save_lagrangian_multipliers( double* lagrangian_multipliers_ );           /** function to save the lagrangian multipliers */
    void compute_lagrangian_multipliers( double* lagrangian_multipliers_, 
                                        const double* constraints_ );              /** computation of the lagrangian multipliers */
    
    void compute_constraints( double* constraints, const double* X_, 
                            const double* U_ );                                    /** computation of the inequality constraints */
    void compute_squared_lateral_distance_vector( double* squared_distances_, 
                            const double* X_);                                      /** computes a vector of the squared lateral distance 
                                                                                        between the trajectory and the allowed center 
                                                                                        lines at each time step*/
    double compute_cost_function(const double* X_, const double* U_);               /** computes the cost function for the vehicle */
    
    double compute_lagrangian(const double* X_, const double* U_);                    /** computes of the augmented lagrangian  
                                                                                    L = cost + lagrangian_multipliers * constraints */
    void compute_gradient(double* gradient, const double* U_);                      /** computes the gradient of the lagrangian with respect to 
                                                                                    U */
    void quadratic_problem_solver(Eigen::MatrixXd & s_, 
                                const Eigen::MatrixXd & G_, 
                                const Eigen::MatrixXd & H_, double Delta);          /** it solves the quadratic problem 
                                                                                        (GT * s + 0.5 * sT * H * s) with solution included in the 
                                                                                        trust region ||s|| < Delta */
    void constraints_diagnostic(const double* constraints, bool print);             /** shows violated constraints */
    void print_trajectories(const double* X, const double* U);                      /** prints trajectories */
    double compute_acceleration(const tk::spline & spline_v, double t);              /** computes the acceleration on the spline s(t) at time t*/
    VehicleState set_prediction(const double* X_, const double* U_);                /** sets the prediction to the traffic structure */
    double compute_heading(const tk::spline & spline_x, 
                           const tk::spline & spline_y, double s);                  /** computes the heading on the spline x(s) and y(s) at parameter s */
    double compute_curvature(const tk::spline & spline_x, 
                             const tk::spline & spline_y, double s);                 /** computes the curvature on the spline x(t) and y(t) at time t*/
    double gradient_norm(const double* gradient);                                               /** computes the norm of the gradient */
    void correctionU(double* U_);                                                    /** corrects U if outside the boundaries */
};

#endif // DYNAMIC_GAME_PLANNER_H
