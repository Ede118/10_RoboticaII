% Requiere Robotics Toolbox for MATLAB (Peter Corke)
% https://petercorke.com/toolboxes/robotics-toolbox/

clc; close all;


%% parametros del robot

% uso esta notacion de denavit-hartenberg (es la del informe)
% las unidades son en [mm]
DH = [
    0       0       0        pi/2   0;   % Joint 1
    0       0       150      0      0;   % Joint 2
    0       0       150      0      0];  % Joint 3
%  theta    d       a      alpha  sigma

name = 'VenganzaDeLaFATEFI';
qlim = deg2rad([ ...
   -180     180;       % qlim1
   -180     180;       % qlim2
   -135     135]);     % qlim3


offset = deg2rad([
    180; 
    90; 
    -90]);

base = transl(0,0,0);
d_tool = 0; 
tool = transl(0,0,d_tool);

%% armo el objeto del robot (SerialLink)

Robot = SerialLink(DH);
Robot.name = name;
Robot.qlim = qlim;
Robot.offset = offset;
Robot.base = base;
Robot.tool = tool;

% Mantener variable 'robot' (con minúscula) para compatibilidad
robot = Robot;

%% limites del plot (workspace)
x1lim = -350;
x2lim = 350;
y1lim = -350;
y2lim = 350;
z1lim = -100;
z2lim = 350;

WS = [x1lim x2lim y1lim y2lim z1lim z2lim];