addpath('/home/sebastien/Documents/EPFL/MA1/MPC')
addpath('/home/sebastien/Documents/EPFL/MA1/MPC/tbxmanager')
import casadi.*

tbxmanager restorepath
x = sdpvar
optimize (abs ( x ) <= 1 , (x -3)*(x-3))
value(x)