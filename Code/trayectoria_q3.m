% TRAYECTORIA_Q3  Solo se mueve q3.
% El eslabon 2 queda horizontal (a lo largo de x0) y el eslabon 3 se dobla
% 90 grados hacia arriba. Perfil quintico (velocidad y aceleracion nulas en
% los extremos). Genera en el workspace t (Nx1), q, dq, ddq (Nx3) y calcula
% los torques.
clear; clc;

%% Parametros de la trayectoria
T      = 2;                 % duracion [s]
dt     = 0.01;              % paso de tiempo [s]
offset = [0 0 0];           % IGUAL que P.offset en dinamica_inversa.m
th0    = [0 0 0];           % angulos DH iniciales [rad]
thf    = [0 0 pi/2];        % angulos DH finales   [rad]
calcular_torque = true;     % false = solo generar la trayectoria

%% Perfil quintico normalizado
t   = (0:dt:T).';
tn  = t/T;
s   = 10*tn.^3 - 15*tn.^4 + 6*tn.^5;
ds  = (30*tn.^2 - 60*tn.^3 + 30*tn.^4)/T;
dds = (60*tn - 180*tn.^2 + 120*tn.^3)/T^2;

%% Trayectoria articular (q = theta_DH - offset)
q   = th0 + (thf - th0).*s - offset;     % [rad]
dq  = (thf - th0).*ds;                   % [rad/s]
ddq = (thf - th0).*dds;                  % [rad/s^2]

%% Graficos
figure;
subplot(3,1,1); plot(t, q,   'LineWidth', 1.2); grid on; ylabel('q [rad]');
legend('q_1','q_2','q_3'); title('Solo q_3');
subplot(3,1,2); plot(t, dq,  'LineWidth', 1.2); grid on; ylabel('dq [rad/s]');
subplot(3,1,3); plot(t, ddq, 'LineWidth', 1.2); grid on; ylabel('ddq [rad/s^2]');
xlabel('t [s]');

%% Torques (dinamica inversa)
if calcular_torque
    tau = dinamica_inversa(t, q, dq, ddq);
end
