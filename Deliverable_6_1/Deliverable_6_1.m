addpath(fullfile('..', 'src'));

close all
clear all
clc

Ts = 1/20;
rocket = Rocket(Ts);

H = 1.5; % Horizon length in seconds
nmpc = NmpcControl(rocket, H);

% MPC reference with default maximum roll = 15 deg
% ref = @(t_, x_) ref_TVC(t_);

% MPC reference with specified maximum roll = 50 deg
roll_max = deg2rad(50);
ref = @(t_, x_) ref_TVC(t_, roll_max);

%% Evaluate once and plot optimal open−loop trajectory,
% pad last input to get consistent size with time and state
x0 = zeros(12, 1);
ref4 = [2 2 2 deg2rad(40)]';
[u, T_opt, X_opt, U_opt] = nmpc.get_u(x0, ref4);
U_opt(:,end+1) = nan;
ph = rocket.plotvis(T_opt, X_opt, U_opt, ref4);

%% Simulate the T%% Simulate close loop system and track TVC
Tf = 30;
[T, X, U, Ref] = rocket.simulate(x0, Tf, @nmpc.get_u, ref);

% Visualize
rocket.anim_rate = 2; % Increase this to make the animation faster
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = 'Merged lin. MPC in nonlinear simulation'; % Set a figure title