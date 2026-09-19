%% ENERGIA CINETICA SIMBOLICA - ROBOT RRR 3 GDL
% Convenciones usadas:
%   - w0(:,i): velocidad angular absoluta del eslabon i expresada en {0}
%   - wb(:,i): velocidad angular absoluta del eslabon i expresada en un
%              sistema solidario al propio eslabon
%   - vc(:,i): velocidad lineal del centro de masa expresada en {0}
%   - I_ci:    tensor de inercia respecto del centro de masa, expresado
%              en el sistema solidario al eslabon i
%
% IMPORTANTE:
% Para usar directamente un tensor de inercia obtenido de CAD, dicho tensor
% debe estar referido al centro de masa y expresado en el mismo sistema
% solidario utilizado para wb(:,i). Si el frame CAD tiene otra orientacion,
% debe aplicarse una rotacion CONSTANTE entre ambos frames.

clear;
clc;

%% Variables simbolicas
syms q1 q2 q3 q_dot_1 q_dot_2 q_dot_3 real
syms m1 m2 m3 l_c_2 l_c_3 L2 real

% Tensor de inercia del eslabon 1
syms I_1_xx I_1_xy I_1_xz I_1_yy I_1_yz I_1_zz real
I_c1 = [I_1_xx I_1_xy I_1_xz;
        I_1_xy I_1_yy I_1_yz;
        I_1_xz I_1_yz I_1_zz];

% Tensor de inercia del eslabon 2
syms I_2_xx I_2_xy I_2_xz I_2_yy I_2_yz I_2_zz real
I_c2 = [I_2_xx I_2_xy I_2_xz;
        I_2_xy I_2_yy I_2_yz;
        I_2_xz I_2_yz I_2_zz];

% Tensor de inercia del eslabon 3
syms I_3_xx I_3_xy I_3_xz I_3_yy I_3_yz I_3_zz real
I_c3 = [I_3_xx I_3_xy I_3_xz;
        I_3_xy I_3_yy I_3_yz;
        I_3_xz I_3_yz I_3_zz];

% Matrices auxiliares
w0 = sym(zeros(3,3));   % velocidades angulares en {0}
wb = sym(zeros(3,3));   % velocidades angulares en frame solidario
vc = sym(zeros(3,3));   % velocidades de centros de masa en {0}

% K(1,i): traslacional
% K(2,i): rotacional
% K(3,i): total
K = sym(zeros(3,3));

%% K1
w0(:,1) = [0; 0; q_dot_1];
wb(:,1) = w0(:,1);
vc(:,1) = sym([0;0;0]);

K(1,1) = sym(1)/2 * m1 * (vc(:,1).' * vc(:,1));
K(2,1) = sym(1)/2 * (wb(:,1).' * I_c1 * wb(:,1));
K(3,1) = simplificar(K(1,1) + K(2,1));

fprintf('\n================ K1 ================\n');
disp('omega_1 en {0}:'); disp(w0(:,1));
disp('K1:'); disp(K(3,1));

%% K2
% Velocidad angular relativa de 2 respecto de 1, expresada en {0}
w21 = q_dot_2 * [sin(q1); -cos(q1); 0];
w0(:,2) = simplificar(w0(:,1) + w21);

% Sistema solidario al eslabon 2:
% x2: longitudinal al eslabon
% z2: eje de la articulacion 2
% y2 = z2 x x2 para cerrar terna dextrógira
x2 = [-cos(q1)*sin(q2);
      -sin(q1)*sin(q2);
       cos(q2)];

z2 = [ sin(q1);
      -cos(q1);
       0];

y2 = simplificar(cross(z2,x2));

R02 = simplificar([x2 y2 z2]);

% Velocidad angular expresada en el frame solidario al eslabon 2
wb(:,2) = simplificar(R02.' * w0(:,2));

% Vector O2 -> C2 expresado en {0}
rO2C2 = l_c_2 * x2;

% O2 es fijo para la geometria adoptada
vc(:,2) = simplificar(cross(w0(:,2), rO2C2));

K(1,2) = simplificar(sym(1)/2 * m2 * (vc(:,2).' * vc(:,2)));
K(2,2) = simplificar(sym(1)/2 * (wb(:,2).' * I_c2 * wb(:,2)));
K(3,2) = simplificar(K(1,2) + K(2,2));

% Ordenar K2 segun velocidades articulares
q_dot_2vec = [q_dot_1; q_dot_2];
[K(3,2), A2] = quadraticCollect(K(3,2), q_dot_2vec);

% Expresion teorica esperada de K2 (tensor general en frame solidario)
K2_ref = ...
    sym(1)/2 * ( ...
        m2*l_c_2^2*sin(q2)^2 ...
        + I_2_xx*cos(q2)^2 ...
        + I_2_yy*sin(q2)^2 ...
        - 2*I_2_xy*sin(q2)*cos(q2) ) * q_dot_1^2 ...
    + sym(1)/2 * (m2*l_c_2^2 + I_2_zz) * q_dot_2^2 ...
    + (I_2_xz*cos(q2) - I_2_yz*sin(q2)) * q_dot_1*q_dot_2;

check_K2 = verificarCero(K(3,2) - K2_ref);

fprintf('\n================ K2 ================\n');
disp('omega_2 en {0}:'); disp(w0(:,2));
disp('omega_2 en frame solidario {2}:'); disp(wb(:,2));
disp('v_C2:'); disp(vc(:,2));
disp('K2 traslacional:'); disp(K(1,2));
disp('K2 rotacional:'); disp(K(2,2));
disp('K2 total ordenado:'); disp(K(3,2));
disp('Chequeo K2 - K2_ref (debe ser 0):'); disp(check_K2);

%% K3
% Eje 3 paralelo y con igual sentido que eje 2
w32 = q_dot_3 * [sin(q1); -cos(q1); 0];
w0(:,3) = simplificar(w0(:,2) + w32);

q23 = q2 + q3;

% Sistema solidario al eslabon 3
x3 = [-cos(q1)*sin(q23);
      -sin(q1)*sin(q23);
       cos(q23)];

z3 = [ sin(q1);
      -cos(q1);
       0];

y3 = simplificar(cross(z3,x3));
R03 = simplificar([x3 y3 z3]);

% Velocidad angular expresada en el frame solidario al eslabon 3
wb(:,3) = simplificar(R03.' * w0(:,3));

% Geometria
rO2O3 = L2 * [-cos(q1)*sin(q2);
              -sin(q1)*sin(q2);
               cos(q2)];

rO3C3 = l_c_3 * x3;

% Velocidad de O3 y del centro de masa C3
v_O3  = simplificar(cross(w0(:,2), rO2O3));
v_rel3 = simplificar(cross(w0(:,3), rO3C3));
vc(:,3) = simplificar(v_O3 + v_rel3);

K(1,3) = simplificar(sym(1)/2 * m3 * (vc(:,3).' * vc(:,3)));
K(2,3) = simplificar(sym(1)/2 * (wb(:,3).' * I_c3 * wb(:,3)));
K(3,3) = simplificar(K(1,3) + K(2,3));

% Ordenar K3 segun velocidades articulares
q_dot = [q_dot_1; q_dot_2; q_dot_3];
[K(3,3), A3] = quadraticCollect(K(3,3), q_dot);

% Expresion teorica esperada de K3
radial_C3 = L2*sin(q2) + l_c_3*sin(q23);
H3 = I_3_xz*cos(q23) - I_3_yz*sin(q23);

K3_ref = ...
    sym(1)/2 * ( ...
        m3*radial_C3^2 ...
        + I_3_xx*cos(q23)^2 ...
        + I_3_yy*sin(q23)^2 ...
        - 2*I_3_xy*sin(q23)*cos(q23) ) * q_dot_1^2 ...
    + sym(1)/2 * ( ...
        m3*(L2^2 + l_c_3^2 + 2*L2*l_c_3*cos(q3)) ...
        + I_3_zz ) * q_dot_2^2 ...
    + sym(1)/2 * (m3*l_c_3^2 + I_3_zz) * q_dot_3^2 ...
    + ( ...
        m3*(l_c_3^2 + L2*l_c_3*cos(q3)) ...
        + I_3_zz ) * q_dot_2*q_dot_3 ...
    + H3 * q_dot_1*q_dot_2 ...
    + H3 * q_dot_1*q_dot_3;

check_K3 = verificarCero(K(3,3) - K3_ref);

fprintf('\n================ K3 ================\n');
disp('omega_3 en {0}:'); disp(w0(:,3));
disp('omega_3 en frame solidario {3}:'); disp(wb(:,3));
disp('v_O3:'); disp(v_O3);
disp('v_C3:'); disp(vc(:,3));
disp('K3 traslacional:'); disp(K(1,3));
disp('K3 rotacional:'); disp(K(2,3));
disp('K3 total ordenado:'); disp(K(3,3));
disp('Chequeo K3 - K3_ref (debe ser 0):'); disp(check_K3);

%% Energia cinetica total
K_total = simplificar(K(3,1) + K(3,2) + K(3,3));
[K_total, A] = quadraticCollect(K_total, q_dot);

% Matriz de inercia del robot: K = 1/2*qdot.'*M(q)*qdot
M = hessian(K_total, q_dot);
M = simplificarMatriz(M);

fprintf('\n================ K TOTAL ================\n');
disp('Energia cinetica total K:');
disp(K_total);
fprintf('\nLaTeX de K:\n%s\n', latex(K_total));

fprintf('\nMatriz de inercia M(q):\n');
disp(M);
fprintf('\nLaTeX de M(q):\n%s\n', latex(M));

% Coeficientes de la forma cuadratica K = qdot^T*A*qdot
coef = {
    'q_dot_1^2',       A(1,1);
    'q_dot_2^2',       A(2,2);
    'q_dot_3^2',       A(3,3);
    'q_dot_1 q_dot_2', 2*A(1,2);
    'q_dot_1 q_dot_3', 2*A(1,3);
    'q_dot_2 q_dot_3', 2*A(2,3)
    };

for k = 1:size(coef,1)
    fprintf('\nCoeficiente de %s:\n', coef{k,1});
    disp(coef{k,2});
    fprintf('LaTeX: %s\n', latex(coef{k,2}));
end

%% Caso particular: ejes principales de inercia
% Si el frame solidario coincide con los ejes principales, los productos
% de inercia se anulan.
productos_inercia = [I_2_xy I_2_xz I_2_yz I_3_xy I_3_xz I_3_yz];
ceros_productos   = zeros(1,6);

K2_principal = simplificar(subs(K(3,2), productos_inercia, ceros_productos));
K3_principal = simplificar(subs(K(3,3), productos_inercia, ceros_productos));
K_principal  = simplificar(subs(K_total, productos_inercia, ceros_productos));

fprintf('\n================ CASO EJES PRINCIPALES ================\n');
disp('K2 con productos de inercia nulos:'); disp(K2_principal);
disp('K3 con productos de inercia nulos:'); disp(K3_principal);
disp('K total con productos de inercia nulos:'); disp(K_principal);

%% FUNCIONES AUXILIARES
function expr = simplificar(expr)
% Simplificacion pensada para conservar expresiones trigonometricas compactas.
% Se evita expand() porque suele destruir terminos como sin(q2+q3).

    expr = rewrite(expr, 'sincos');
    expr = simplify(expr, 'Steps', 200);

    try
        expr = combine(expr, 'sincos');
    catch
        % Algunas versiones/objetos simbolicos pueden no admitir combine.
    end

    expr = simplify(expr, 'Steps', 200);
end

function A = simplificarMatriz(A)
% Simplifica elemento por elemento para evitar salidas innecesariamente largas.

    for i = 1:size(A,1)
        for j = 1:size(A,2)
            A(i,j) = simplificar(A(i,j));
        end
    end
end

function [expr_ordered, A] = quadraticCollect(expr, q_dot)
% Para una forma cuadratica:
%       expr = q_dot.' * A * q_dot
% obtiene A mediante la Hessiana y reconstruye expr ordenada.

    expr = simplify(expr, 'Steps', 200);
    A = hessian(expr, q_dot) / 2;
    A = simplificarMatriz(A);

    n = length(q_dot);
    expr_ordered = sym(0);

    for i = 1:n
        expr_ordered = expr_ordered + A(i,i)*q_dot(i)^2;

        for j = i+1:n
            expr_ordered = expr_ordered ...
                + 2*A(i,j)*q_dot(i)*q_dot(j);
        end
    end

    expr_ordered = simplificar(expr_ordered);
end

function z = verificarCero(expr)
% Rutina separada para verificar identidades. Aqui SI se permite expandir,
% porque no se usa para presentar el resultado sino para comprobar igualdad.

    z = simplify(expand(rewrite(expr,'sincos')), 'Steps', 500);
end
