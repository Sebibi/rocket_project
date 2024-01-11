addpath(fullfile('..', 'src'));

close all
clear all
clc

Ts = 1/20;
rocket = Rocket(Ts);

nmpc = NmpcControl(rocket, H);

%% Evaluate once and plot optimal open−loop trajectory,
% pad last input to get consistent size with time and state
x0 = zeros(12, 1);
ref4 = [2 2 2 deg2rad(40)]';
[u, T_opt, X_opt, U_opt] = nmpc.get_u(x0, ref4);
U_opt(:,end+1) = nan;
ph = rocket.plotvis(T_opt, X_opt, U_opt, ref4);

% MPC reference with default maximum roll = 15 deg
ref = @(t_, x_) ref_TVC(t_);

%% Simulate close loop system and track TVC
Tf = 30;
[T, X, U, Ref] = rocket.simulate(x0, Tf, @nmpc.get_u, ref);
% Visualize
rocket.anim_rate = 2; % Increase this to make the animation faster
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = 'NMPC in nonlinear simulation roll_max=50'; % Set a figure title


% MPC reference with specified maximum roll = 50 deg
roll_max = deg2rad(50);
ref = @(t_, x_) ref_TVC(t_, roll_max);

%% Simulate close loop system and track TVC
Tf = 30;
[T, X, U, Ref] = rocket.simulate(x0, Tf, @nmpc.get_u, ref);
% Visualize
rocket.anim_rate = 2; % Increase this to make the animation faster
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = 'NMPC in nonlinear simulation roll_max=50'; % Set a figure title


%% Simulate the linear controller on TVC roll_max = 50
% Linearize the system
[xs, us] = rocket.trim();
sys = rocket.linearize(xs, us);
[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us);

% Setup linear controller
H = 2; % Horizon length in seconds
mpc_x = MpcControl_x(sys_x, Ts, H);
mpc_y = MpcControl_y(sys_y, Ts, H);
mpc_z = MpcControl_z(sys_z, Ts, H);
mpc_roll = MpcControl_roll(sys_roll, Ts, H);
mpc = rocket.merge_lin_controllers(xs, us, mpc_x, mpc_y, mpc_z, mpc_roll);

% Simulate and visualize closed-loop system 
Tf = 30;
[T, X, U, Ref] = rocket.simulate(x0, Tf, @mpc.get_u, ref);
rocket.anim_rate = 2; % Increase this to make the animation faster
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = 'Merged linear MPC in nonlinear simulation roll_max=50'; 