clc;
clear all;
close all;


%% todo1

Ts = 1/20;
rocket = Rocket(Ts);

d1 = 1;
d2 = 0;
Pavg = 65;
Pdiff = 0;

u = [d1, d2, Pavg, Pdiff]'; % (Assign appropriately)
[b_F, b_M] = rocket.getForceAndMomentFromThrust(u)

w = [0, 0, 0];
phi = [0 0 0];
v = [0 0 1];
p = [0 0 0];

x = [w, phi, v, p]'; % (Assign appropriately)
x_dot = rocket.f(x, u)

%% todo2
% info = rendererinfo(gca);

rocket = Rocket(Ts);
Tf = 5.0; % Simulation end time
%x0 = [deg2rad([2 -2 0, -2 2 0]), 0 0 0, 0 0 0]'; % (w, phi, v, p) Initial state
x0 = x
u = [deg2rad([d1 d2]), Pavg, Pdiff ]'; % (d1 d2 Pavg Pdiff) Constant input
[T, X, U] = rocket.simulate(x0, Tf, u); % Simulate unknown, nonlinear model
rocket.anim_rate = 1.0; % Visualize at 1.0x real−time
rocket.vis(T, X, U);

%% todo3
rocket = Rocket(Ts);
[xs, us] = rocket.trim() % Compute steady−state for which 0 = f(xs,us)
sys = rocket.linearize(xs, us) % Linearize the nonlinear model about trim point
[T, X, U] = rocket.simulate(xs, Tf, us); % Simulate unknown, nonlinear model
rocket.anim_rate = 1.0; % Visualize at 1.0x real−time
rocket.vis(T, X, U);

%% todo4
[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us)


