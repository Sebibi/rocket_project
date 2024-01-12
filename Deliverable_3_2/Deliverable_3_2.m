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

%% Initial state for the controllers 
x0_x = [0; 0; 0; 0]; % System x
x0_y = [0; 0; 0; 0]; %System y
x0_z = [0; 0]; % System z
x0_roll = [0; deg2rad(0)]; % System roll

%% Open loop 
plot_open_loop(x0_x, mpc_x, sys_x, xs, us, rocket, -4);
plot_open_loop(x0_y, mpc_y, sys_y, xs, us, rocket, -4);
plot_open_loop(x0_z, mpc_z, sys_z, xs, us, rocket, -4);
plot_open_loop(x0_roll, mpc_roll, sys_roll, xs, us, rocket, deg2rad(35));

%% Closed Loop
Tf = 10;
plot_closed_loop(Tf, x0_x, sys_x, mpc_x, xs, us, rocket, -4);
plot_closed_loop(Tf, x0_y, sys_y, mpc_y, xs, us, rocket, -4);
plot_closed_loop(Tf, x0_z, sys_z, mpc_z, xs, us, rocket, -4);
plot_closed_loop(Tf, x0_roll, sys_roll, mpc_roll, xs, us, rocket, deg2rad(35));

%% Function definition
function plot_open_loop(x0, mpc, sys, xs, us, rocket, pos_ref)
    [u, T_opt, X_opt, U_opt] = mpc.get_u(x0, pos_ref);
    U_opt(:,end+1) = NaN;
    X_opt = X_opt + sys.UserData.xs;
    U_opt = U_opt + sys.UserData.us;
    ph = rocket.plotvis_sub(T_opt, X_opt, U_opt, sys, xs, us); % Plot as usual
end

function plot_closed_loop(Tf, x0, sys, mpc, xs, us, rocket, pos_ref)
    [T, X_sub, U_sub] = rocket.simulate_f(sys, x0, Tf, @mpc.get_u, pos_ref);
    ph = rocket.plotvis_sub(T, X_sub, U_sub, sys, xs, us);
end