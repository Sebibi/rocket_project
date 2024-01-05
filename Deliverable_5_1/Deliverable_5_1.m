addpath(fullfile('..', 'src'));

close all
clear all
clc

%% TODO: This file should produce all the plots for the deliverable

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

%% Code du prof todo 5.1
% Merge four sub−system controllers into one full−system controller
mpc = rocket.merge_lin_controllers(xs, us, mpc_x, mpc_y, mpc_z, mpc_roll);

% z0 = [0; 0];
% 
% uz = mpc_z.get_u(z0)

%%
x0 = [zeros(1, 9), 1 0 3]';
ref = [1.2, 0, 3, 0]';
Tf = 8;


%% Manipulate mass for simulation
rocket.mass = 2.13;

%% 
%[T, X, U, Ref] = rocket.simulate(x0, Tf, @mpc.get_u, ref);

[T, X, U, Ref, Z_hat] = rocket.simulate_est_z(x0, Tf, @mpc.get_u, ref, mpc_z, sys_z)






