function tau = dinamica_inversa(t, q, dq, ddq)
% DINAMICA_INVERSA  Torques articulares de un robot serie RRR de 3 gdl.
%
% Formulacion: Le Ngoc Truc, Nguyen Van Quyen, Nguyen Phung Quang,
% "Dynamic model with a new formulation of Coriolis/centrifugal matrix for
% robot manipulators", J. Comput. Sci. Cybern. 36(1), 2020, pp. 89-104.
% Se usan las Secs. 2, 3 y el inicio de la Sec. 4 (hasta la Ec. 46):
%   - velocidades de eslabon:   S(w_i) = R0i' * dR0i/dt          (Ec. 10)
%   - jacobianos:               JT_i = dp0Ci/dq,  JR_i = dw_i/ddq (Ec. 24-25)
%   - matriz de inercia:        M = sum(m_i JT'JT + JR' I_i JR)   (Ec. 28)
%   - energia potencial:        P = -sum(m_i g0' p0Ci)            (Ec. 17)
%   - gravedad:                 g = (dP/dq)'                      (Ec. 45)
%   - vector de Coriolis:       C*dq segun la Ec. 46 (Kronecker)
%   - modelo:                   tau = M*ddq + C*dq + g            (Ec. 44)
%
% ENTRADAS (una fila por instante de tiempo):
%   t   : N x 1   vector de tiempo [s]
%   q   : N x 3   posiciones articulares [rad]
%   dq  : N x 3   velocidades articulares [rad/s]
%   ddq : N x 3   aceleraciones articulares [rad/s^2]
% SALIDA:
%   tau : N x 3   torque de cada articulacion [N*m]  (tambien se grafica)
%
% Requiere Symbolic Toolbox. Completar la seccion PARAMETROS.

% Modo verificacion: llamar  dinamica_inversa  sin argumentos. Construye el
% modelo y deja en el workspace M, Cdq, g (simbolicos) y fM, fCdq, fg.
if nargin == 0
    [fM, fCdq, fg, M, Cdq, g] = construir_modelo();
    assignin('base', 'fM', fM);     assignin('base', 'fCdq', fCdq);
    assignin('base', 'fg', fg);     assignin('base', 'M', M);
    assignin('base', 'Cdq', Cdq);   assignin('base', 'g', g);
    return
end
[fM, fCdq, fg] = construir_modelo();

N   = numel(t);
tau = zeros(N, 3);
for k = 1:N
    qk   = q(k,:).';
    dqk  = dq(k,:).';
    ddqk = ddq(k,:).';
    M    = fM(qk(1), qk(2), qk(3));
    Cdq  = fCdq(qk(1), qk(2), qk(3), dqk(1), dqk(2), dqk(3));
    g    = fg(qk(1), qk(2), qk(3));
    tau(k,:) = (M*ddqk + Cdq + g).';
end

figure;
plot(t, tau, 'LineWidth', 1.2); grid on;
xlabel('t [s]'); ylabel('\tau [N m]');
legend('\tau_1', '\tau_2', '\tau_3');
title('Torques articulares - dinamica inversa');
end


%% ========================================================================
%  PARAMETROS DEL ROBOT  (completar; NaN = pendiente)
%  ========================================================================
function P = parametros()
% --- Tabla DH (convencion clasica: Rz(theta) Tz(d) Tx(a) Rx(alpha)) ------
d1 = 0.03;            % [m] altura del eslabon 1 (d1 de la tabla DH)
L2 = 0.124;            % [m] a2 (longitud del eslabon 2)
L3 = 0.128;            % [m] a3 (longitud del eslabon 3)

P.d      = [d1, 0, 0];
P.a      = [0, L2, L3];
P.alpha  = [-sym(pi)/2, 0, 0];   % simbolico para que cos(pi/2) sea exactamente 0
P.offset = [0, 0, 0];           % theta_i = q_i + offset_i. Ej.: si tu q2 se
                                % mide desde la vertical, prueba [0, pi/2, 0]

% --- Gravedad en el sistema base (z0 hacia arriba) -----------------------
P.g0 = [0; 0; -9.807];

% --- Masas [kg] ----------------------------------------------------------
P.m = [0.069, 0.044, 0.097];

% --- Centro de masa de cada eslabon, expresado en SU PROPIO frame i ------
% (vector desde el origen del frame i hasta el CM, en ejes del frame i).
% En el paper: rC2 = [r1 - l1; 0; 0], es decir, [lc - a; 0; 0], donde lc es
% la distancia al CM medida desde el inicio del eslabon (origen del frame
% anterior) y a la longitud del eslabon.
P.rC = {[0; -0.013; 0], ...    % eslabon 1
        [-0.056; 0; 0], ...    % eslabon 2
        [-0.104; 0; 0]};       % eslabon 3

% --- Tensor de inercia respecto del CM, en ejes paralelos al frame i -----
% (Ec. 19-22 del paper: I_i es el tensor en el CM, con ejes del frame i).
% Unidades kg*m^2. Los productos de inercia van en cero si los ejes del
% frame coinciden con los ejes principales.
P.I = { [2956 0 0; 0 2711 0; 0 0 2521]*(1e-10), ...   % eslabon 1
        [6214 0 0; 0 20224 0; 0 0 14473]*(1e-10), ...   % eslabon 2
        [11963 0 0; 0 68870 0; 0 0 58373]*(1e-10) };     % eslabon 3

% --- Control de que no quede nada sin completar --------------------------
rc  = vertcat(P.rC{:});
II  = vertcat(P.I{:});
chk = [P.d(:); P.a(:); P.m(:); rc(:); II(:)];
if any(isnan(chk))
    error('dinamica_inversa:parametros', ...
          'Faltan parametros por completar en parametros().');
end
end


%% ========================================================================
%  CONSTRUCCION SIMBOLICA DEL MODELO
%  ========================================================================
function [fM, fCdq, fg, M, Cdq, g] = construir_modelo()
n = 3;
P = parametros();

syms q1 q2 q3 dq1 dq2 dq3 real
q  = [q1; q2; q3];
dq = [dq1; dq2; dq3];

theta = q + P.offset(:);

M  = sym(zeros(n));
Pot = sym(0);
T0 = sym(eye(4));

for i = 1:n
    % Transformacion homogenea 0 -> i  (Ec. 15)
    T0 = T0 * dh(theta(i), P.d(i), P.a(i), P.alpha(i));
    R  = simplify(T0(1:3,1:3));

    % Velocidad angular en el frame i: S(w_i) = R' * dR/dt   (Ec. 10)
    dR = sym(zeros(3));
    for k = 1:n
        dR = dR + diff(R, q(k)) * dq(k);
    end
    S = simplify(R.' * dR);
    w = [S(3,2); S(1,3); S(2,1)];            % vector axial de S(w)
    JR = jacobian(w, dq);                    % Ec. 25

    % Posicion del CM en la base: [p;1] = T0i*[rC;1]          (Ec. 14)
    pC = T0 * [sym(P.rC{i}); 1];
    pC = simplify(pC(1:3));
    JT = jacobian(pC, q);                    % Ec. 24

    % Matriz de inercia generalizada                          (Ec. 28)
    M = M + P.m(i) * (JT.' * JT) + JR.' * sym(P.I{i}) * JR;

    % Energia potencial                                       (Ec. 17)
    Pot = Pot - P.m(i) * (sym(P.g0).' * pC);
end
M = simplify(M);

% Vector de gravedad                                          (Ec. 45)
g = simplify(jacobian(Pot, q).');

% dM/dq segun la Definicion 2 (Ec. 31): matriz n x (n*n), donde el bloque
% de la columna j de M contiene las derivadas respecto de q1..qn.
dMdq = sym(zeros(n, n*n));
for j = 1:n
    for k = 1:n
        dMdq(:, (j-1)*n + k) = diff(M(:, j), q(k));
    end
end

% Vector de Coriolis/centrifugo                               (Ec. 46)
In  = eye(n);
Cdq = dMdq * kron(In, dq) * dq - 0.5 * (dMdq * kron(dq, In)).' * dq;
Cdq = simplify(Cdq);

% Funciones numericas rapidas
fM   = matlabFunction(M,   'Vars', {q1, q2, q3});
fCdq = matlabFunction(Cdq, 'Vars', {q1, q2, q3, dq1, dq2, dq3});
fg   = matlabFunction(g,   'Vars', {q1, q2, q3});
end


function T = dh(th, d, a, al)
% Transformacion DH clasica: Rz(th) * Tz(d) * Tx(a) * Rx(al)
T = [cos(th), -sin(th)*cos(al),  sin(th)*sin(al), a*cos(th);
     sin(th),  cos(th)*cos(al), -cos(th)*sin(al), a*sin(th);
     0,        sin(al),          cos(al),         d;
     0,        0,                0,               1];
end
