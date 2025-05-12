#include "trust_region_mpc_solver.h"
#include <iostream>

TrustRegionMPCSolver::TrustRegionMPCSolver(VehicleState& vs) : vehicle_state(vs) 
{
}

TrustRegionMPCSolver::~TrustRegionMPCSolver() {
    // Destructor implementation
    std::cout << "TrustRegionMPCSolver destroyed." << std::endl;
}

void TrustRegionMPCSolver::run() {
    // Variables initialization and setup:
    setup();

    // definition of the control variable vector U and of the state vector X:
    double U[nu];
    double X[nx];
    double constraints[nC];

    initial_guess(X, U);
    trust_region_solver(U);
    integrate(X, U);
    print_trajectories(X, U);
    compute_constraints(constraints, X, U);
    constraints_diagnostic(constraints, false);
    set_prediction(X, U);
}

void TrustRegionMPCSolver::setup() {
    
    // Setup number of inequality constraints for the vehicle:
    // 2 * nU * (N + 1) inequality constraints for inputs 
    // (N + 1) constraints to remain in the lane
    nC = 2 * nU * (N + 1) + (N + 1);
    nG = nu;

    // resize and initialize limits for control input
    ul.resize(nU * (N + 1), 1);
    uu.resize(nU * (N + 1), 1);
    for (int j = 0; j < N + 1; j++){
        ul(nU * j + d, 0) = d_low;
        uu(nU * j + d, 0) = d_up;
        ul(nU * j + F, 0) = F_low;
        uu(nU * j + F, 0) = F_up;
    }

    // resize and initialize lagrangian multiplier vector
    lagrangian_multipliers.resize(nC, 1);
    lagrangian_multipliers = Eigen::MatrixXd::Zero(nC, 1);

}

/** Sets the intial guess of the game */
void TrustRegionMPCSolver::initial_guess(double* X_, double* U_)
{
    for (int j = 0; j < N + 1; j++){
        U_[nU * j + d] = 0.0;
        U_[nU * j + F] = 0.0;
    }
    integrate(X_, U_);
}

/** integrates the input U to get the state X */
void TrustRegionMPCSolver::integrate(double* X_, const double* U_)
{
    int tu;
    int td;
    double s_ref;
    double t;
    double s_t0[nX];
    double sr_t0[nX];
    double u_t0[nU];
    double ds_t0[nX];
    double ds_t1[nX];

    t = 0.0;

    // Initial state:
    s_t0[x] = vehicle_state.x;
    s_t0[y] = vehicle_state.y;
    s_t0[v] = vehicle_state.v;
    s_t0[psi] = vehicle_state.psi;
    s_t0[s] = 0.0;
    s_t0[l] = 0.0;

    ds_t0[x] = vehicle_state.v * std::cos(vehicle_state.psi + vehicle_state.beta);
    ds_t0[y] = vehicle_state.v * std::sin(vehicle_state.psi + vehicle_state.beta);
    ds_t0[v] = vehicle_state.a;
    ds_t0[psi] = vehicle_state.omega;
    ds_t0[s] = 0.0;
    ds_t0[l] = 0.0;

    for (int j = 0; j < N + 1; j++){
        tu = nU * j;
        td = nX * j;

        // Reference point on the center lane:
        s_ref = s_t0[s];
        sr_t0[x] = vehicle_state.centerlane.spline_x(s_ref);
        sr_t0[y] =vehicle_state.centerlane.spline_y(s_ref);
        sr_t0[psi] = vehicle_state.centerlane.compute_heading(s_ref);
        sr_t0[v] = vehicle_state.v + j * (vehicle_state.v_target - vehicle_state.v) / (N + 1); 

        // Input control:
        u_t0[d] = U_[tu + d];
        u_t0[F] = U_[tu + F];
        
        // Dynamic step: 
        bicycle_dynamics_step(ds_t1, ds_t0, s_t0, sr_t0, u_t0);

        // Integration to compute the new state: 
        s_t0[x] += dt * ds_t1[x];
        s_t0[y] += dt * ds_t1[y];
        s_t0[v] += dt * ds_t1[v];
        s_t0[psi] += dt * ds_t1[psi];
        s_t0[s] += dt * ds_t1[s];
        s_t0[l] += dt * ds_t1[l];

        ds_t0[x] = ds_t1[x];
        ds_t0[y] = ds_t1[y];
        ds_t0[v] = ds_t1[v];
        ds_t0[psi] = ds_t1[psi];
        ds_t0[s] = ds_t1[s];
        ds_t0[l] = ds_t1[l];

        if (s_t0[v] < 0.0){s_t0[v] = 0.0;}

        // Save the state in the trajectory
        X_[td + x] = s_t0[x];
        X_[td + y] = s_t0[y];
        X_[td + v] = s_t0[v];
        X_[td + psi] = s_t0[psi];
        X_[td + s] = s_t0[s];
        X_[td + l] = s_t0[l];

        t+= dt;
    }
    
}

/** Dyanamic step */
void TrustRegionMPCSolver::bicycle_dynamics_step(double* d_state, const double* d_state_old, const double* state, const double* ref_state, const double* control)
{
    double throttle_val = 0.0;
    double brake_val = 0.0;
    double u = control[F];            // signed command
    double slope = 5.0;               // slope of the switch; bigger = sharper

    double sigmoid  = 0.5 * (1.0 + std::tanh(slope*u));   // ∈ (0,1)

    throttle_val =  sigmoid *  u;          // positive half
    brake_val    = (1.0 - sigmoid) * -u;   // negative half, always ≥ 0
    
    // -------- Vehicle longitudinal dynamics --------
    double engine_force = throttle_val * max_engine_force;
    double brake_force = brake_val * max_brake_force;
    double drag_force = 0.5 * air_density * drag_coeff * frontal_area * state[v] * state[v];
    //double rolling_force = rolling_resistance_coeff * mass * g;
    double rolling_force = 0.0;

    double total_force = engine_force - brake_force - drag_force - rolling_force;
    
    // -------- Vehicle lateral/slip dynamics --------
    double beta = control[d] * Lr / L;
    double yaw_rate = d_state_old[psi];

    if (state[v] > 1.0)
    {
        double alpha_f = control[d] - std::atan2((yaw_rate * Lf + state[v] * std::sin(beta)), (state[v] * std::cos(beta) + 1e-6));
        double alpha_r = - std::atan2((yaw_rate * Lr - state[v] * std::sin(beta)), (state[v] * std::cos(beta) + 1e-6));

        double Fyf = - D_f * std::sin(C_f * std::atan(B_f * alpha_f - E_f * (B_f * alpha_f - std::atan(B_f * alpha_f))));
        double Fyr = - D_r * std::sin(C_r * std::atan(B_r * alpha_r - E_r * (B_r * alpha_r - std::atan(B_r * alpha_r))));

        double yaw_acc = (Lf * Fyf - Lr * Fyr) / Iz;
        yaw_rate = yaw_rate + yaw_acc * dt;

        double lat_acceleration = (Fyf + Fyr) / mass;
        double v_lat = state[v] * std::sin(beta) + lat_acceleration * dt;

        beta = std::atan2(v_lat, state[v] + 1e-6);
    }else{
        yaw_rate = state[v] * std::tan(control[d]) * cos((Lr / L) * control[d])/ L;
    }

    /* Derivatives computation:*/
    d_state[x] = state[v] * cos(state[psi] + beta);
    d_state[y] = state[v] * sin(state[psi] + beta);
    d_state[v] = total_force / mass;
    d_state[psi] = yaw_rate;
    d_state[l] =  weight_target_speed * (state[v] - ref_state[v]) * (state[v] - ref_state[v])
            + weight_center_lane * ((ref_state[x] - state[x]) * (ref_state[x] - state[x]) + (ref_state[y] - state[y]) * (ref_state[y] - state[y]))
            + weight_heading * ((std::cos(ref_state[psi]) - std::cos(state[psi]))*(std::cos(ref_state[psi]) - std::cos(state[psi]))
            +        (std::sin(ref_state[psi]) - std::sin(state[psi]))*(std::sin(ref_state[psi]) - std::sin(state[psi])));
    d_state[s] = state[v];
}

/** SR1 Hessian matrix update*/
void TrustRegionMPCSolver::hessian_SR1_update(Eigen::MatrixXd & H_, const Eigen::MatrixXd & s_, const Eigen::MatrixXd & y_, double r_)
{
    if (abs((s_.transpose() * (y_ - H_ * s_))(0,0)) > 
                r_ * (s_.transpose() * s_)(0,0) * ((y_ - H_* s_).transpose() * (y_ - H_ * s_))(0,0))
    {
        H_ = H_ + ((y_ - H_ * s_) * (y_ - H_ * s_).transpose()) / ((y_ - H_ * s_).transpose() * s_)(0,0);
    }
}

/** function to increase rho = rho * gamma at each iteration */
void TrustRegionMPCSolver::increasing_schedule()
{
    rho = gamma * rho;
}

/** function to save the lagrangian multipliers in the general variable */
void TrustRegionMPCSolver::save_lagrangian_multipliers(double* lagrangian_multipliers_)
{
    for (int i = 0; i < nC; i++){
        lagrangian_multipliers(i,0) = lagrangian_multipliers_[i];
    }
}

/* computation of lambda (without update)*/
void TrustRegionMPCSolver::compute_lagrangian_multipliers(double* lagrangian_multipliers_, const double* constraints_)
{
    double l;
    for (int i = 0; i < nC; i++){
        l = lagrangian_multipliers(i,0) + rho * constraints_[i];
        lagrangian_multipliers_[i] = std::max(l, 0.0);
    }
}

/** computation of the inequality constraints (target: constraints < 0) */
void TrustRegionMPCSolver::compute_constraints(double* constraints, const double* X_, const double* U_)
{
    int ind = 0;
    int indCu;
    int indCl;
    int indU;
    int indCto;
    int indClk;
    int indCuk;
    int indf;
    int n1;
    int n2;
    double dist2t[N + 1];
    double rad2[N + 1];
    double latdist2t[N + 1];
    double r_lane_ = r_lane;

    // constraints for the inputs 
    indU = 0;
    indCu = nU * (N + 1);
    indCl = indCu + nU * (N + 1);
    for (int k = 0; k < N + 1; k++){
        constraints[nU * k + d] = 1e3 * (U_[indU + nU * k + d] - uu(nU * k + d,0));
        constraints[nU * k + F] = 1e3 * (U_[indU + nU * k + F] - uu(nU * k + F,0));
        constraints[indCu + nU * k + d] = 1e3 * (ul(nU * k + d,0) - U_[indU + nU * k + d]);
        constraints[indCu + nU * k + F] = 1e3 * (ul(nU * k + F,0) - U_[indU + nU * k + F]);
    }

    indCto = indCl;

    // constraints to remain in the lane
    compute_squared_lateral_distance_vector(latdist2t, X_);
    for (int k = 0; k < N + 1; k++){
        constraints[indCto + k] = (latdist2t[k] - r_lane_ * r_lane_);
    }
    indf = indCto + (N + 1);
}

/** computes a vector of the squared lateral distance between the i-th trajectory and the allowed center lines at each time step*/
void TrustRegionMPCSolver::compute_squared_lateral_distance_vector(double* squared_distances_, const double* X_)
{
    double s_;
    double x_;
    double y_;
    double x_c;
    double y_c;
    double x_r;
    double y_r;
    double x_l;
    double y_l;
    double psi_c;
    double psi_l;
    double psi_r;
    double dist_c;
    double dist_l = 1e3;
    double dist_r = 1e3;
    double dist_long_c;
    double dist_long_l;
    double dist_long_r;
    double dist2_c[N + 1];
    double dist2_l[N + 1];
    double dist2_r[N + 1];
    double dist2_rl_min;
    for (int j = 0; j < N + 1; j++){
        s_ = X_[nX * j + s];
        x_ = X_[nX * j + x];
        y_ = X_[nX * j + y];
        dist2_c[j] = 1e3;
        dist2_l[j] = 1e3;
        dist2_r[j] = 1e3;
        if (s_ < vehicle_state.centerlane.s_max){
            x_c = vehicle_state.centerlane.spline_x(s_);
            y_c = vehicle_state.centerlane.spline_y(s_);
            psi_c = vehicle_state.centerlane.compute_heading(s_);
            dist_c = ((x_ - x_c) * (x_ - x_c) + (y_ - y_c) * (y_ - y_c));
            dist_long_c = ((x_ - x_c) * std::cos(psi_c) + (y_ - y_c) * std::sin(psi_c)) * ((x_ - x_c) * std::cos(psi_c) + (y_ - y_c) * std::sin(psi_c));
            dist2_c[j] = dist_c - dist_long_c;
        }
        if (vehicle_state.leftlane.present == true && s_ < vehicle_state.leftlane.s_max && vehicle_state.leftlane.s_max > 10.0){
            x_l = vehicle_state.leftlane.spline_x(s_);
            y_l = vehicle_state.leftlane.spline_y(s_);
            psi_l = vehicle_state.leftlane.compute_heading(s_);
            dist_l = ((x_ - x_l) * (x_ - x_l) + (y_ - y_l) * (y_ - y_l));
            dist_long_l = ((x_ - x_l) * std::cos(psi_l) + (y_ - y_l) * std::sin(psi_l)) * ((x_ - x_l) * std::cos(psi_l) + (y_ - y_l) * std::sin(psi_l));
            dist2_l[j] = dist_l - dist_long_l;
        }
        if (vehicle_state.rightlane.present == true && s_ < vehicle_state.rightlane.s_max && vehicle_state.rightlane.s_max > 10.0){
            x_r = vehicle_state.rightlane.spline_x(s_);
            y_r = vehicle_state.rightlane.spline_y(s_);
            psi_r = vehicle_state.rightlane.compute_heading(s_);
            dist_r = ((x_ - x_r) * (x_ - x_r) + (y_ - y_r) * (y_ - y_r));
            dist_long_r = ((x_ - x_r) * std::cos(psi_r) + (y_ - y_r) * std::sin(psi_r)) * ((x_ - x_r) * std::cos(psi_r) + (y_ - y_r) * std::sin(psi_r));
            dist2_r[j] = dist_r - dist_long_r;
        }
        dist2_rl_min = std::min(dist2_l[j], dist2_r[j]);
        squared_distances_[j] = std::min(dist2_rl_min, dist2_c[j]);
    }
}

/** compute the cost function */
double TrustRegionMPCSolver::compute_cost_function(const double* X_, const double* U_)
{
    double final_lagrangian = X_[nX * N + l];
    double cost = 0.5 * final_lagrangian * qf * final_lagrangian;
    return cost;
}

/** computes of the augmented lagrangian vector  L = <L_1, ..., L_M> L_i = cost_i + lagrangian_multipliers * constraints */
double TrustRegionMPCSolver::compute_lagrangian(const double* X_, const double* U_)
{
    double cost;
    double constraints[nC];
    
    cost = compute_cost_function( X_, U_);
    compute_constraints(constraints, X_, U_);
    double lagrangian = cost;
    double constraint;
    for (int k = 0; k < nC; k++){
        constraint = std::max(0.0, constraints[k]);
        lagrangian += 0.5 * rho * constraint * constraint + lagrangian_multipliers(k,0) * constraints[k];
    }
    
    return lagrangian;
}

/** computation of the gradient of lagrangian_i with respect to U_i for each i with parallelization on cpu*/
void TrustRegionMPCSolver::compute_gradient(double* gradient, const double* U_)
{
    const int num_threads = std::thread::hardware_concurrency();
    std::mutex mutex;
    std::vector<std::thread> threads(num_threads);

    // Definition of the work for each thread:
    auto computeGradient = [&](int start, int end) {
        double dU[nu];
        double dX[nx];
        double X_[nx];
        double lagrangian;
        double d_lagrangian;
        double cost;
        double constraints[nC];
        int index;
        for (int i = 0; i < nu; i++){
            dU[i] = U_[i];
        }
        integrate(X_, U_);
        lagrangian = compute_lagrangian(X_, U_);
        for (int i = start; i < end; i++) {
            dU[i] = U_[i] + eps;
            integrate(dX, dU);
            compute_constraints(constraints, dX, dU);
            cost = compute_cost_function( dX, dU);
            d_lagrangian = compute_lagrangian( dX, dU);
            {
                std::lock_guard<std::mutex> lock(mutex);
                gradient[i] = (d_lagrangian - lagrangian) / eps;
            }
            dU[i] = U_[i];
        }
    };

    // Parallelize:
    int work_per_thread = nu / num_threads;
    int start_index = 0;
    int end_index = 0;
    for (int i = 0; i < num_threads; ++i) {
        start_index = i * work_per_thread;
        end_index = (i == num_threads - 1) ? nu : start_index + work_per_thread;
        threads[i] = std::thread(computeGradient, start_index, end_index);
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

/** it solves the quadratic problem (GT * s + 0.5 * sT * H * s) with solution included in the trust region ||s|| < Delta */
void TrustRegionMPCSolver::quadratic_problem_solver(Eigen::MatrixXd & s_, const Eigen::MatrixXd & G_, const Eigen::MatrixXd & H_, double Delta)
{
    Eigen::MatrixXd ps(nG,1);
    double tau;
    double tau_c;
    double normG;
    double GTHG;
    GTHG = (G_.transpose() * H_ * G_)(0,0);
    normG = sqrt((G_.transpose() * G_)(0,0));
    ps = - Delta * (G_ /(normG));
    if ( GTHG <= 0.0){
        tau = 1.0;
    }else{
        tau_c = (normG * normG * normG)/(Delta * GTHG);
        tau = std::min(tau_c, 1.0);
    }
    s_ = tau * ps;
}

/** prints if some constraints are violated */
void TrustRegionMPCSolver::constraints_diagnostic(const double* constraints, bool print = false)
{
    bool flag0 = false;
    bool flag1 = false;
    bool flag2 = false;
    bool flag3 = false;
    flag0 = false;
    flag1 = false;
    flag2 = false;
    flag3 = false;
    for (int j = 0; j < nC; j++){
        if (constraints[j] > 0){
            if (j < (2 * nU * (N + 1))) {
                std::cerr<<"The vehicle violates input constraints: "<<constraints[j]<<"\n";
                flag0 = true;
            }
            if (j > (2 * nU * (N + 1))){
                std::cerr<<"The vehicle violates lane constraints: "<<constraints[j]<<"\n";
                flag2 = true;
            }
        }
    }
    if (print == true){
        std::cerr<<"input constraint: \n";
        for (int j = 0; j < 2 * nU * (N + 1); j++){
            std::cerr<<constraints[j]<<"\t";
        }
        std::cerr<<"\nlane constraint: \n";
        for (int j = (2 * nU * (N + 1)); j < (2 * nU * (N + 1) + (N + 1)); j++){
            std::cerr<<constraints[j]<<"\t";
        }
        std::cerr<<"\n";
    }
    
}

void TrustRegionMPCSolver::print_trajectories(const double* X, const double* U)
{
    const int col_width = 12;  // Adjust as needed

    // Print header
    std::cerr << std::left
              << std::setw(col_width) << "Time"
              << std::setw(col_width) << "X"
              << std::setw(col_width) << "Y"
              << std::setw(col_width) << "V"
              << std::setw(col_width) << "PSI"
              << std::setw(col_width) << "S"
              << std::setw(col_width) << "L"
              << std::setw(col_width) << "F"
              << std::setw(col_width) << "d"
              << "\n";

    // Print separator
    std::cerr << std::string(col_width * 9, '-') << "\n";

    // Print trajectory with time
    for (int j = 0; j < N + 1; j++) {
        double t = j * dt;
        std::cerr << std::fixed << std::left
                  << std::setw(col_width) << t
                  << std::setw(col_width) << X[nX * j + x]
                  << std::setw(col_width) << X[nX * j + y]
                  << std::setw(col_width) << X[nX * j + v]
                  << std::setw(col_width) << X[nX * j + psi]
                  << std::setw(col_width) << X[nX * j + s]
                  << std::setw(col_width) << X[nX * j + l]
                  << std::setw(col_width) << U[nU * j + F]
                  << std::setw(col_width) << U[nU * j + d]
                  << "\n";
    }

    std::cerr << "\n";
}

/** sets the prediction to the traffic structure*/
void TrustRegionMPCSolver::set_prediction(const double* X_, const double* U_) 
{
    Trajectory trajectory;
    Control control;
    double time = 0.0;

    for (int j = 0; j < N; j++) {
        TrajectoryPoint point;
        Input input;
        input.T = U_[nU * j + F];
        input.delta = U_[nU * j + d];
        point.x = X_[nX * j + x];
        point.y = X_[nX * j + y];
        point.psi = X_[nX * j + psi];
        point.v = X_[nX * j + v];
        point.a = X_[nX * (j + 1) + v] - X_[nX * j + v] / dt;
        point.omega = X_[nX * (j + 1) + psi] - X_[nX * j + psi] / dt;
        point.t_start = time;
        point.t_end = time + dt;
        trajectory.push_back(point);
        control.push_back(input);
        time += dt;
    }

    // Now update the actual member variable
    vehicle_state.predicted_trajectory = trajectory;
    vehicle_state.predicted_control = control;
}

/** computes the norm of the gradient */
double TrustRegionMPCSolver::gradient_norm(const double* gradient)
{
    double norm = 0.0;
    for (int j = 0; j < nU; j++){
        norm += gradient[j] * gradient[j];
    }
    return norm;
}

/** Trust-Region solver of the dynamic game*/
void TrustRegionMPCSolver::trust_region_solver(double* U_)
{
    bool convergence = false;

    // Parameters:
    double eta = 1e-4;
    double r_ = 1e-8;
    double threshold_gradient_norm = 1e-2;
    int iter = 1;
    int iter_lim = 20;

    // Variables definition:
    double gradient[nu];
    double dU[nu]; 
    double dU_[nu]; 
    double dX[nx]; 
    double dX_[nx];
    double d_gradient[nu];
    double d_lagrangian;
    double lagrangian;
    double constraints[nC];
    double lagrangian_multipliers[nC];

    double actual_reduction;
    double predicted_reduction;
    double delta;
    Eigen::MatrixXd H_;
    Eigen::MatrixXd g_;
    Eigen::MatrixXd p_;
    Eigen::MatrixXd s_;
    Eigen::MatrixXd y_;

    // Variables initialization:
    integrate(dX, U_);
    for (int i = 0; i < nu; i++){
        dU[i] = U_[i];
        dU_[i] = U_[i];
    }
    for (int i = 0; i < nx; i++){
        dX_[i] = dX[i];
    }
    H_.resize(nu, nu);
    g_.resize(nu, 1);
    p_.resize(nu, 1);
    s_.resize(nu, 1);
    y_.resize(nu, 1);
    delta = 1.0;
    H_ = Eigen::MatrixXd::Identity(nu, nu);
    
    compute_gradient(gradient, dU_);
    
    // Check for convergence:
    if (gradient_norm(gradient) < threshold_gradient_norm){
        convergence = true;
    }

    // Iteration loop:
    while (convergence == false && iter < iter_lim ){

        // Compute the grandient and the lagrangian
        integrate(dX_, dU_);
        compute_gradient(gradient, dU_);
        lagrangian = compute_lagrangian(dX_, dU_);

        // Solves the quadratic subproblem and compute the possible step dU:
        for (int j = 0; j < N + 1; j++){
            g_(j * nU + d,0) = gradient[j * nU + d];
            g_(j * nU + F,0) = gradient[j * nU + F];
        }
        quadratic_problem_solver(s_, g_, H_, delta);
        for (int j = 0; j < N + 1; j++){
            dU[j * nU + d] = dU_[j * nU + d] + s_(j * nU + d,0);
            dU[j * nU + F] = dU_[j * nU + F] + s_(j * nU + F,0);
        }
        

        // Compute the new grandient and the new lagrangian with the possible step dU:
        integrate(dX, dU);
        compute_gradient(d_gradient, dU);
        d_lagrangian = compute_lagrangian(dX, dU);

        // Check if to accept the step:
        // Compute the actual reduction and of the predicted reduction:
        actual_reduction = lagrangian - d_lagrangian;
        predicted_reduction = - (g_.transpose() *  s_ + 0.5 * s_.transpose() * H_ * s_)(0,0);

        // In case of very low or negative actual reduction, reject the step:
        if ( actual_reduction / predicted_reduction < eta){ 
            for (int j = 0; j < N + 1; j++){
                dU[j * nU + d] = dU_[j * nU + d];
                dU[j * nU + F] = dU_[j * nU + F];
            }
        }

        // In case of great reduction, and solution close to the trust region, increase the trust region:
        if ( actual_reduction / predicted_reduction > 0.75){ 
            if (std::sqrt((s_.transpose() * s_)(0,0)) > 0.8 * delta){
                delta = 2.0 * delta;
            }
        }

        // In case of low actual reduction, decrease the step:
        if ( actual_reduction / predicted_reduction < 0.1){
            delta = 0.5 * delta;
        }

        // Compute the difference of the gradients, then the Hessian matrix update:
        for (int j = 0; j < N + 1; j++){
            y_(j * nU + d,0) = d_gradient[j * nU + d] - gradient[j * nU + d];
            y_(j * nU + F,0) = d_gradient[j * nU + F] - gradient[j * nU + F];
        }
        hessian_SR1_update(H_, s_, y_, r_);

        // Save the solution for the next iteration:
        for (int j = 0; j < N + 1; j++){
            dU_[j * nU + d] = dU[j * nU + d];
            dU_[j * nU + F] = dU[j * nU + F];
        }
        
        // Check for convergence:
        if (gradient_norm(gradient) < threshold_gradient_norm){
            convergence = true;
        }

        // Compute the new state: 
        integrate(dX_, dU_);

        // Compute the constraints with the new solution:
        compute_constraints(constraints, dX_, dU_);

        // Compute and save in the general variable the lagrangian multipliers with the new solution:
        compute_lagrangian_multipliers(lagrangian_multipliers, constraints);
        save_lagrangian_multipliers(lagrangian_multipliers);

        // Increase the weight of the constraints in the lagrangian multipliers:
        increasing_schedule();
        iter++;
    }

    std::cerr<<"number of iterations: "<<iter<<"\n";

    //Correct the final solution:
    correctionU(dU_);

    // Save the solution:
    for(int k = 0; k < nu; k++){
        U_[k] = dU_[k];
    }
}

void TrustRegionMPCSolver::correctionU(double* U_)
{
    for (int j = 0; j < N + 1; j++){
        if (j == N){
            U_[nU * j + d] = U_[nU * (j - 1) + d];
            U_[nU * j + F] = U_[nU * (j - 1) + F];
        }
        if (U_[nU * j + d] > d_up){
            U_[nU * j + d] = d_up;
        }
        if (U_[nU * j + d] < d_low){
            U_[nU * j + d] = d_low;
        }
    }
}