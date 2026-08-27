function LeyTemporal(guardar_cartesiano, guardar_articular, guardar_jacobiano, guardar_video)
% Requiere Robotics Toolbox for MATLAB (Peter Corke)
% https://petercorke.com/toolboxes/robotics-toolbox/

arguments (Input)
    guardar_cartesiano = false;
    guardar_articular = false;
    guardar_jacobiano = false;
    guardar_video = false;
end

%% ============================================================
% PARÁMETROS GEOMÉTRICOS Y TEMPORALES
% =============================================================

L = 0.10;      % Longitud recorrida por T1 y T4 [m]
R = 0.025;     % Radio / amplitud [m]
H = 2*R;       % Descenso total de T2 [m]

pasos = 400;   % Cantidad de puntos por tramo

% Punto inicial de toda la trayectoria
P0 = [-L, 4*R, 0];

% Vector de tiempo interpolado (ley temporal trapezoidal lspb)
[u, ud, udd] = tpoly(0, 1, pasos); 

%% ============================================================
% GENERACIÓN DE TRAMOS CON LEY TEMPORAL Y VECTORES DE AVANCE
% =============================================================

% --- Tramo 1 (T1) ---
% Avance L en X. Cuarto de onda sinusoidal en Y y Z.
theta1 = (pi/2) * u;
xT1 = P0(1) + (2*L/pi)*theta1;
yT1 = P0(2) + R*sin(theta1);
zT1 = P0(3) + R*sin(theta1);
% Vector de avance tangente
t_adv_T1 = [ (2*L/pi)*ones(size(theta1)), R*cos(theta1), R*cos(theta1) ];

P1 = [xT1(end), yT1(end), zT1(end)];

% --- Tramo 2 (T2) ---
% Circunferencia completa en XY. Centro XY: C2 = [P1_x, P1_y - R]
% Descenso H en Z.
theta2 = 2*pi * u;
xT2 = P1(1) + R*sin(theta2);
yT2 = P1(2) - R + R*cos(theta2);
zT2 = P1(3) - (H/2)*(1 - cos(theta2/2));
% Vector de avance tangente
t_adv_T2 = [ R*cos(theta2), -R*sin(theta2), -(H/4)*sin(theta2/2) ];

P2 = [xT2(end), yT2(end), zT2(end)];

% --- Tramo 3 (T3) ---
% Desplazamiento recto de -2R en Y.
s3 = u;
xT3 = P2(1)*ones(size(s3));
yT3 = P2(2) - 2*R*s3;
zT3 = P2(3)*ones(size(s3));
% Vector de avance tangente
t_adv_T3 = [ zeros(size(s3)), -2*R*ones(size(s3)), zeros(size(s3)) ];

P3 = [xT3(end), yT3(end), zT3(end)];

% --- Tramo 4 (T4) ---
% Avance L en X. Incremento R en Y y Z con concavidad opuesta a T1.
theta4 = (pi/2) * u;
xT4 = P3(1) + (2*L/pi)*theta4;
yT4 = P3(2) + R*(1 - cos(theta4));
zT4 = P3(3) + R*(1 - cos(theta4));
% Vector de avance tangente
t_adv_T4 = [ (2*L/pi)*ones(size(theta4)), R*sin(theta4), R*sin(theta4) ];

P4 = [xT4(end), yT4(end), zT4(end)];


%% ============================================================
% CONCATENACIÓN DE LOS TRAMOS
% =============================================================

Trayectoria_X = [xT1; xT2; xT3; xT4];
Trayectoria_Y = [yT1; yT2; yT3; yT4];
Trayectoria_Z = [zT1; zT2; zT3; zT4];
Trayectoria_Adv   = [t_adv_T1; t_adv_T2; t_adv_T3; t_adv_T4];

Total_Pasos = length(Trayectoria_X);


%% ============================================================
% BLENDING DE LAS ESQUINAS (SUAVIZADO)
% =============================================================

porcentaje_empalme = 0.02;
ancho_empalme = round(pasos * porcentaje_empalme);

% Índices de transiciones entre tramos
transiciones = [pasos, 2*pasos, 3*pasos];

% Máscara gaussiana para el blending (1 en la esquina, 0 lejos)
peso = zeros(Total_Pasos, 1);
for idx = transiciones
    rango = max(1, idx - ancho_empalme) : min(Total_Pasos, idx + ancho_empalme);
    w = exp(-0.5 * ((rango - idx) / (ancho_empalme / 3)).^2);
    peso(rango) = max(peso(rango), w');
end

% Suavizado de la posición en las esquinas
Suave_X = smoothdata(Trayectoria_X, 'gaussian', 2 * ancho_empalme);
Suave_Y = smoothdata(Trayectoria_Y, 'gaussian', 2 * ancho_empalme);
Suave_Z = smoothdata(Trayectoria_Z, 'gaussian', 2 * ancho_empalme);

Trayectoria_X = (1 - peso) .* Trayectoria_X + peso .* Suave_X;
Trayectoria_Y = (1 - peso) .* Trayectoria_Y + peso .* Suave_Y;
Trayectoria_Z = (1 - peso) .* Trayectoria_Z + peso .* Suave_Z;

% Suavizado del vector de avance
Suave_Adv = smoothdata(Trayectoria_Adv, 1, 'gaussian', 2 * ancho_empalme);
Trayectoria_Adv = (1 - peso) .* Trayectoria_Adv + peso .* Suave_Adv;

% Normalización de los vectores de avance
Trayectoria_Adv = Trayectoria_Adv ./ vecnorm(Trayectoria_Adv, 2, 2);


%% ============================================================
% POSICIÓN DEL VASO (SOLO POSICIÓN CARTESIANA)
% =============================================================

% Dado que el robot tiene 3 GDL, solo se requiere la posición X, Y, Z
CPosition = [Trayectoria_X, Trayectoria_Y, Trayectoria_Z];


%% ============================================================
% CINEMÁTICA INVERSA
% =============================================================

% Cargar el modelo del robot
RobotSerie;

% Semilla inicial (Front - Elbow Up)
q_semilla_inicial = [0, pi/4, -pi/4]; 

[Q_middle] = CinematicaInversa(Robot, CPosition, q_semilla_inicial);

% Homing para acercamiento y retirada
q_home = [0, 0, 0]; 
pasos_homing = 50; 
tiempo_homing = 5; % segundos

% Interpolación de trayectorias de homing
Q_approach = jtraj(q_home, Q_middle(1,:), pasos_homing);
Q_retreat  = jtraj(Q_middle(end,:), q_home, pasos_homing);

% Unión de todas las secuencias
Q = [Q_approach; Q_middle(2:end, :); Q_retreat(2:end, :)];
Q_animacion = Q; % Se define para evitar el error del script original

fprintf('Matriz Target_Poses (%dx6) generada con éxito.\n', size(Q, 1));


%% ============================================================
% ANÁLISIS EN EL ESPACIO CARTESIANO
% =============================================================

tiempo_por_tramo = 60; 
t = linspace(0, tiempo_por_tramo, pasos);
dt = t(2) - t(1); 

Tramos_X = {xT1, xT2, xT3, xT4};
Tramos_Y = {yT1, yT2, yT3, yT4};
Tramos_Z = {zT1, zT2, zT3, zT4};
Nombres = {'Tramo 1 (Subida)', 'Tramo 2 (Círculo)', 'Tramo 3 (Recta Y)', 'Tramo 4 (Salida)'};

carpeta_destino_C = 'Graficos_Cinematica_Cartesiana_LT';
if ~exist(carpeta_destino_C, 'dir')
    mkdir(carpeta_destino_C);
end

for i = 1:4
    X = Tramos_X{i};
    Y = Tramos_Y{i};
    Z = Tramos_Z{i};

    % Derivadas de posición para obtener velocidad
    Vx = gradient(X, dt);
    Vy = gradient(Y, dt);
    Vz = gradient(Z, dt);

    % Derivadas de velocidad para obtener aceleración
    Ax = gradient(Vx, dt);
    Ay = gradient(Vy, dt);
    Az = gradient(Vz, dt);

    % Gráfico de posición
    fig1 = figure('Color', 'w', 'Name', ['Posicion - ' Nombres{i}]);
    plot(t, [X Y Z], 'LineWidth', 1.5);
    title(['Posición Cartesiana - ' Nombres{i}]);
    ylabel('Posición [m]'); xlabel('Tiempo [s]');
    lgdX = legend('X', 'Y', 'Z', 'Location', 'eastoutside'); 
    lgdX.ItemHitFcn = @toggleSignal;
    grid on; grid minor;

    % Gráfico de velocidad
    fig2 = figure('Color', 'w', 'Name', ['Velocidad - ' Nombres{i}]);
    plot(t, [Vx Vy Vz], 'LineWidth', 1.5);
    title(['Velocidad Cartesiana - ' Nombres{i}]);
    ylabel('Velocidad [m/s]'); xlabel('Tiempo [s]');
    lgdV = legend('V_x', 'V_y', 'V_z', 'Location', 'eastoutside'); 
    lgdV.ItemHitFcn = @toggleSignal;
    grid on; grid minor;

    % Gráfico de aceleración
    fig3 = figure('Color', 'w', 'Name', ['Aceleración - ' Nombres{i}]);
    plot(t, [Ax Ay Az], 'LineWidth', 1.5);
    title(['Aceleración Cartesiana - ' Nombres{i}]);
    ylabel('Acel. [m/s^2]'); xlabel('Tiempo [s]');
    lgdA = legend('A_x', 'A_y', 'A_z', 'Location', 'eastoutside'); 
    lgdA.ItemHitFcn = @toggleSignal;
    grid on; grid minor;
    
    if guardar_cartesiano
        nombre_pos = fullfile(carpeta_destino_C, sprintf('Tramo_%d_Posicion.png', i));
        exportgraphics(fig1, nombre_pos, 'Resolution', 300);

        nombre_vel = fullfile(carpeta_destino_C, sprintf('Tramo_%d_Velocidad.png', i));
        exportgraphics(fig2, nombre_vel, 'Resolution', 300);

        nombre_acc = fullfile(carpeta_destino_C, sprintf('Tramo_%d_Aceleracion.png', i));
        exportgraphics(fig3, nombre_acc, 'Resolution', 300);
    end
end

if guardar_cartesiano
    disp('Los 12 gráficos cartesianos han sido guardados.');
end


%% ============================================================
% ANÁLISIS EN EL ESPACIO ARTICULAR (MOTORES)
% =============================================================

cant_tramos = 4;
tiempo_traslado = cant_tramos * tiempo_por_tramo; 

t_app = linspace(0, tiempo_homing, pasos_homing)';
t_mid = linspace(tiempo_homing, tiempo_homing + tiempo_traslado, size(Q_middle, 1))';
t_ret = linspace(tiempo_homing + tiempo_traslado, tiempo_homing + tiempo_traslado + tiempo_homing, pasos_homing)';

t_total = [t_app; t_mid(2:end); t_ret(2:end)];
Total_Pasos_Articulares = length(t_total);

V_art = zeros(Total_Pasos_Articulares, 3);
A_art = zeros(Total_Pasos_Articulares, 3);

for j = 1:3
    V_art(:, j) = gradient(Q(:, j), t_total);
    A_art(:, j) = gradient(V_art(:, j), t_total);
end

nombres_ejes = {'q_1 (Base)', 'q_2 (Hombro)', 'q_3 (Codo)'};

% Gráfico de posición articular
fig_q = figure('Color', 'w', 'Name', 'Posición Articular');
plot(t_total, rad2deg(Q), 'LineWidth', 1.5); 
title('Evolución de la Posición Articular');
ylabel('Posición [deg]'); xlabel('Tiempo [s]');
lgdQ = legend(nombres_ejes, 'Location', 'eastoutside'); 
lgdQ.ItemHitFcn = @toggleSignal;
grid on; grid minor;

% Gráfico de velocidad articular
fig_vq = figure('Color', 'w', 'Name', 'Velocidad Articular');
plot(t_total, rad2deg(V_art), 'LineWidth', 1.5);
title('Evolución de la Velocidad Articular');
ylabel('Velocidad [deg/s]'); xlabel('Tiempo [s]');
lgdQd = legend(nombres_ejes, 'Location', 'eastoutside'); 
lgdQd.ItemHitFcn = @toggleSignal;
grid on; grid minor;

% Gráfico de aceleración articular
fig_aq = figure('Color', 'w', 'Name', 'Aceleración Articular');
plot(t_total, rad2deg(A_art), 'LineWidth', 1.5);
title('Evolución de la Aceleración Articular');
ylabel('Aceleración [deg/s^2]'); xlabel('Tiempo [s]');
lgdQdd = legend(nombres_ejes, 'Location', 'eastoutside'); 
lgdQdd.ItemHitFcn = @toggleSignal;
grid on; grid minor;

% Gráfico del determinante del jacobiano
disp('Calculando Determinante del Jacobiano...');
det_J = zeros(Total_Pasos_Articulares, 1);
for j = 1:Total_Pasos_Articulares
    J = Robot.jacob0(Q(j,:));
    det_J(j) = det(J(1:3, :)); % Determinante de la parte de traslación (3x3)
end

fig_m = figure('Color', 'w', 'Name', 'Determinante del Jacobiano');
plot(t_total, det_J, 'LineWidth', 1.5, 'Color', '#D95319');
title('Determinante del Jacobiano Geométrico (det(J))');
ylabel('det(J)'); xlabel('Tiempo [s]');
grid on; grid minor;

carpeta_destino_Q = 'Graficos_Cinematica_Articular_LT';
if ~exist(carpeta_destino_Q, 'dir')
    mkdir(carpeta_destino_Q);
end

if guardar_articular
    nombre_q = fullfile(carpeta_destino_Q, 'Articular_1_Posicion.png');
    exportgraphics(fig_q, nombre_q, 'Resolution', 300);

    nombre_vq = fullfile(carpeta_destino_Q, 'Articular_2_Velocidad.png');
    exportgraphics(fig_vq, nombre_vq, 'Resolution', 300);
    
    nombre_aq = fullfile(carpeta_destino_Q, 'Articular_3_Aceleracion.png');
    exportgraphics(fig_aq, nombre_aq, 'Resolution', 300);

    disp('Los 3 gráficos articulares han sido guardados.');
end

if guardar_jacobiano
    nombre_m = fullfile(carpeta_destino_Q, 'Articular_4_Maniobrabilidad.png');
    exportgraphics(fig_m, nombre_m, 'Resolution', 300);
    disp('Gráfico de Maniobrabilidad guardado.');
end

disp('Presione [ENTER] para continuar con la simulación.');
pause()


%% ============================================================
% ANIMACIÓN DEL ROBOT
% =============================================================

x1lim = -0.5; x2lim = 0.5;
y1lim = -0.5; y2lim = 0.5;
z1lim = -0.2; z2lim = 0.8;
WS = [x1lim x2lim y1lim y2lim z1lim z2lim];

figure('Color', 'w', 'Name', 'Simulación de Traslado de Vaso - Ley Temporal', ...
    'WindowStyle', 'normal', 'Units', 'pixels', 'Position', [100 100 1920 1080]); 
grid on; 
hold on;

% Trazo la trayectoria en rojo y punto verde al inicio
plot3(Trayectoria_X, Trayectoria_Y, Trayectoria_Z, 'r-', 'LineWidth', 2);
plot3(Trayectoria_X(1), Trayectoria_Y(1), Trayectoria_Z(1), 'g.', 'MarkerSize', 20); 

axis equal;
view(135, 25);

disp('Animando trayectoria...');

if guardar_video
    Robot.plot(...
        Q_animacion, ...
        'workspace', WS, ...
        'notiles', ...
        'scale', 0.5, ...
        'jointdiam', 1.0, ...
        'jointlen', 0.8, ...
        'linkcolor', [.2 .2 .2], ...
        'jointcolor', [1 .4 0], ...
        'fps', 60, ...
        'movie', 'Simulacion_Traslado_LT.mp4'...
    );
else
    Robot.plot(...
        Q_animacion, ...
        'workspace', WS, ...
        'notiles', ...
        'scale', 0.5, ...
        'jointdiam', 1.0, ...
        'jointlen', 0.8, ...
        'linkcolor', [.2 .2 .2], ...
        'jointcolor', [1 .4 0], ...
        'fps', 60 ...
    );
end

disp('Fin de la simulación.');

end


%% ============================================================
% FUNCIONES AUXILIARES
% =============================================================

function toggleSignal(~, event)
    % Alterna la visibilidad de la señal al hacer clic en la leyenda
    if strcmp(event.Peer.Visible, 'on')
        event.Peer.Visible = 'off';
    else
        event.Peer.Visible = 'on';
    end
end

function Q = CinematicaInversa(Robot, CPosition, ~)
    N = size(CPosition, 1);
    Q = zeros(N, 3);
    for i = 1:N
        Q(i, :) = FncCinematicaInversa(Robot, CPosition(i, :));
    end
end
