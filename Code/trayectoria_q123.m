% TRAYECTORIA_Q123  Se mueven q1, q2 y q3 al mismo tiempo.
% - q1: la base gira 90 grados (quintico).
% - q2: levanta el brazo de x0 (horizontal) a z0 (vertical) (quintico).
% - q3: el codo se dobla y vuelve a 0 (ida y vuelta con sin^2), asi que el
%       brazo termina estirado y alineado con z0.
% Velocidad nula en los extremos para las tres articulaciones.
% Genera en el workspace t (Nx1), q, dq, ddq (Nx3) y calcula los torques.
clear; clc;

%% Parametros de la trayectoria
T      = 2;                 % duracion [s]
dt     = 0.01;              % paso de tiempo [s]
offset = [0 0 0];           % IGUAL que P.offset en dinamica_inversa.m
th0    = [0 0 0];           % angulos DH iniciales [rad]
thf    = [pi/2 pi/2 0];     % angulos DH finales   [rad]
amp3   = pi/4;              % maximo doblado del codo (q3) [rad]
calcular_torque = true;     % false = solo generar la trayectoria

%% Perfil quintico normalizado (q1 y q2)
t   = (0:dt:T).';
tn  = t/T;
s   = 10*tn.^3 - 15*tn.^4 + 6*tn.^5;
ds  = (30*tn.^2 - 60*tn.^3 + 30*tn.^4)/T;
dds = (60*tn - 180*tn.^2 + 120*tn.^3)/T^2;

th   = th0  + (thf - th0).*s;            % angulos DH
dth  = (thf - th0).*ds;
ddth = (thf - th0).*dds;

%% q3: ida y vuelta  th3 = amp3*(1 - cos(2*pi*t/T))/2
th(:,3)   = th0(3) + amp3*(1 - cos(2*pi*tn))/2;
dth(:,3)  = amp3*pi*sin(2*pi*tn)/T;
ddth(:,3) = amp3*2*pi^2*cos(2*pi*tn)/T^2;

%% Trayectoria articular (q = theta_DH - offset)
q   = th - offset;                       % [rad]
dq  = dth;                               % [rad/s]
ddq = ddth;                              % [rad/s^2]

%% Graficos
figure;
subplot(3,1,1); plot(t, q,   'LineWidth', 1.2); grid on; ylabel('q [rad]');
legend('q_1','q_2','q_3'); title('q_1, q_2 y q_3 simultaneos');
subplot(3,1,2); plot(t, dq,  'LineWidth', 1.2); grid on; ylabel('dq [rad/s]');
subplot(3,1,3); plot(t, ddq, 'LineWidth', 1.2); grid on; ylabel('ddq [rad/s^2]');
xlabel('t [s]');

%% Torques (dinamica inversa)
if calcular_torque
    tau = dinamica_inversa(t, q, dq, ddq);
end
