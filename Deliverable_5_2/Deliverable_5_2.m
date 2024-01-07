addpath(fullfile('..', 'src'));

close all
clear all
clc


%% Setup our system and compute the equilibrim state and corresponding control
Ts = 1/20; % Sample time
rocket = Rocket(Ts);
[xs, us] = rocket.trim();

%% Linearize the system
sys = rocket.linearize(xs, us);
[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us);


%% Setup our controller
H = 2; % Horizon length in seconds
mpc_x = MpcControl_x(sys_x, Ts, H);
mpc_y = MpcControl_y(sys_y, Ts, H);
mpc_z = MpcControl_z(sys_z, Ts, H);
mpc_roll = MpcControl_roll(sys_roll, Ts, H);

% Merge four sub−system controllers into one full−system controller
mpc = rocket.merge_lin_controllers(xs, us, mpc_x, mpc_y, mpc_z, mpc_roll);

% Setup the initial paramertes and the reference 
x0 = [zeros(1, 9), 1 0 3]';
ref = [1.2, 0, 3, 0]';
Tf = 8;

%% Perform simulations

% Manipulate mass for simulation
rocket.mass = 2.13; 
rocket.mass_rate = -0.27;

% Simulate
[T, X, U, Ref] = rocket.simulate(x0, Tf, @mpc.get_u, ref);
ph = rocket.plotvis(T, X, U, Ref);
set(ph.fig, 'Name', 'Simulation after mass manipulation');

% Simulate with disturbance estimation
[T, X, U, Ref, Z_hat] = rocket.simulate_est_z(x0, Tf, @mpc.get_u, ref, mpc_z, sys_z);
ph = rocket.plotvis(T, X, U, Ref);
set(ph.fig, 'Name', 'Simulation with disturbance estimation');


%% PLot the distubance estimation
fig = figure;
plot(Z_hat(13,:))
title("Disturbance Estimation");
set(fig, 'Name', 'Disturbance Estimation');

%% Bonus question
Tf = 20;
[T, X, U, Ref] = rocket.simulate(x0, Tf, @mpc.get_u, ref);
ph = rocket.plotvis(T, X, U, Ref);
set(ph.fig, 'Name', 'Simulation after mass manipulation for 20 sec (BONUS)');