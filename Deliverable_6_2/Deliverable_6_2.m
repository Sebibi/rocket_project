addpath(fullfile('..', 'src'));

close all
clear all
clc

Ts = 1/40; % Higher sampling rate for this part!
delay = 5;

rocket = Rocket(Ts);
H = 2; % Horizon length in seconds
nmpc = NmpcControl(rocket, H);
nmpc_d = NmpcControl(rocket, H, delay);

x0 = zeros(12, 1);
ref = [0.5, 0, 1, deg2rad(65)]';
Tf = 2.5;
rocket.mass = 1.75;
rocket.delay = delay; % 0 if not specified

%% Simulate open-loop system
[u, T_opt, X_opt, U_opt] = nmpc.get_u(x0, ref);
U_opt(:,end+1) = nan;
ph = rocket.plotvis(T_opt, X_opt, U_opt, ref);
set(ph.fig, 'Name', 'Open-Loop simulation');

[u, T_opt, X_opt, U_opt] = nmpc_d.get_u(x0, ref);
U_opt(:,end+1) = nan;
ph = rocket.plotvis(T_opt, X_opt, U_opt, ref);
set(ph.fig, 'Name', 'Open-Loop simulation with delay compensation');

%% Simulate closed-loop system 
[T, X, U, Ref] = rocket.simulate(x0, Tf, @nmpc.get_u, ref);
ph = rocket.plotvis(T, X, U, Ref);
set(ph.fig, 'Name', 'Closed-Loop simulation');

[T, X, U, Ref] = rocket.simulate(x0, Tf, @nmpc_d.get_u, ref);
ph = rocket.plotvis(T, X, U, Ref);
set(ph.fig, 'Name', 'Closed-Loop simulation with delay compensation');

