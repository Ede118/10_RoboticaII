%[text] # LAGRANGIANO DEL ROBOT - ENERGIA CINÉTICA Y ENERGÍA POTENCIAL
%[text] Convenciones usadas:
%[text] - $\\omega\_0(:,i)$: velocidad angular absoluta del eslabon i expresada en {0}
%[text] - $\\omega\_b(:,i)$: velocidad angular absoluta del eslabon i expresada en un sistema solidario al propio eslabon
%[text] - $v\_c(:,i)$: velocidad lineal del centro de masa expresada en {0}
%[text] - $I\_{C\_i}}$: tensor de inercia respecto del centro de masa, expresado en el sistema solidario al eslabon i \
%[text] ## IMPORTANTE
%[text] Para usar directamente un tensor de inercia obtenido de CAD, dicho tensor debe estar referido al centro de masa y expresado en el mismo sistema solidario utilizado para $\\omega\_b(:,i)$. Si el frame CAD tiene otra orientacion, debe aplicarse una rotacion CONSTANTE entre ambos frames.
clear;
clc;
%[text] ## Definición de variables simbólicas
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
%[text] Matriz de Energia cinética
%[text] - K(1, i) = Energia Cinética Lineal del eslabón i
%[text] - K(2, i) = Energia Cinética Rotacional del eslabón i
%[text] - K(3, i) = Energia Cinética Total del eslabón i \
% K(1,i): traslacional
% K(2,i): rotacional
% K(3,i): total
K = sym(zeros(3,3));
%%
%[text] ## Energía cinética $K\_1$
w(:, 1)  = q_dot_1*[0; 0; 1];
vc(:, 1) = [0; 0; 0];

K(1,1) = sym(1)/2 * m1 * vc(:,1)' * vc(:,1);
K(2,1) = 0.5 * w(:,1)' * I_c1 * w(:,1);      
K(3,1) = K(1, 1) + K(2, 1);                  

disp('omega_1 en {0}:'); disp(w(:,1)); %[output:30006ce4] %[output:77bea529]
disp('v_c_1'); disp(vc(:,1)); %[output:7bf32629] %[output:7bf08550]
disp('K1 lineal'); disp(K(1,1)); %[output:4a9c2d46] %[output:61f85703]
disp('K1 rotacional'); disp(K(2,1)); %[output:1e7fea41] %[output:6dd8d755]
disp('K1 total'); disp(K(3,1)); %[output:435fef61] %[output:54c612c2]
%%
%[text] ## Energía cinética $K\_2$
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
% mismo que rO2C2 = l_c_2 * [-cos(q1)*sin(q2); -sin(q1)*sin(q2); cos(q2)];
rO2C2 = l_c_2 * x2;

% O2 es fijo para la geometria adoptada
vc(:,2) = simplificar(cross(w0(:,2), rO2C2));

K(1,2) = simplificar(sym(1)/2 * m2 * (vc(:,2).' * vc(:,2)));
K(2,2) = simplificar(sym(1)/2 * (wb(:,2).' * I_c2 * wb(:,2)));
K(3,2) = simplificar(K(1,2) + K(2,2));
                
% Simplificación y reordenamiento de términos
q_dot = [q_dot_1; q_dot_2];
[best, A2] = simplifyCompact(K(3,2), 200);
[best, ~] = quadraticCollect(best, q_dot);
K(3,2) = best;

% Expresion teorica esperada de K2 (tensor general en frame solidario)
% derivada de la bibliografia

K2_ref = ...
    sym(1)/2 * ( ...
    m2*l_c_2^2*sin(q2)^2 ...
    + I_2_xx*cos(q2)^2 ...
    + I_2_yy*sin(q2)^2 ...
    - 2*I_2_xy*sin(q2)*cos(q2) ) * q_dot_1^2 ...
    + sym(1)/2 * (m2*l_c_2^2 + I_2_zz) * q_dot_2^2 ...
    + (I_2_xz*cos(q2) - I_2_yz*sin(q2)) * q_dot_1*q_dot_2;

check_K2 = verificarCero(K(3,2) - K2_ref);

disp('omega_2 en {0}:'); disp(w0(:,2)); %[output:5afcd529] %[output:951a5433]
disp('omega_2 en frame solidario {2}:'); disp(wb(:,2)); %[output:21da127d] %[output:6a515821]
disp('v_c_2'); disp(vc(:,2)); %[output:7341bb72] %[output:0bab4640]
disp('K2 traslacional:'); disp(K(1,2)); %[output:67fa81d8] %[output:4bfdbac8]
disp('K2 rotacional:'); disp(K(2,2)); %[output:61624bd3] %[output:4a396191]
disp('K2 total ordenado:'); disp(K(3,2)); %[output:65251f6f] %[output:1d5af13d]
disp('Chequeo K2 - K2_ref (debe ser 0):'); disp(check_K2); %[output:4b0871b6] %[output:8339c208]
%%
%[text] ## Energía cinética $K\_3$
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

% Simplificación y reordenamiento de términos
q_dot = [q_dot_1; q_dot_2; q_dot_3];
[best, A3] = simplifyCompact(K(3,3), 200);
[best, ~] = quadraticCollect(best, q_dot);
K(3,3) = best;

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

disp('omega_3 en {0}:'); disp(w0(:,3)); %[output:2848e350] %[output:454c30c4]
disp('omega_3 en frame solidario {3}:'); disp(wb(:,3)); %[output:76dabcbd] %[output:73ef5ed9]
disp('v_{O3}:'); disp(v_O3); %[output:733a687c] %[output:5d4d816b]
disp('v_c_3:'); disp(vc(:,3)); %[output:359c9dbf] %[output:9252f871]
disp('K3 traslacional:'); disp(K(1,3)); %[output:69e274e4] %[output:62745e20]
disp('K3 rotacional:'); disp(K(2,3)); %[output:866b1500] %[output:37dafd57]
disp('K3 total ordenado:'); disp(K(3,3)); %[output:25aa1212] %[output:9b011b1f]
disp('Chequeo K3 - K3_ref (debe ser 0):'); disp(check_K3); %[output:8a3e050c] %[output:688de766]

%%
%[text] ## Energía Cinética Total $K$
K_total = K(3,1) + K(3,2) + K(3,3);

q_dot = [q_dot_1; q_dot_2; q_dot_3];
[best, ranking] = simplifyCompact(K_total, 200);
[K_total, A] = quadraticCollect(best, q_dot);


% Matriz de inercia del robot: K = 1/2*qdot.'*M(q)*qdot
M = hessian(K_total, q_dot);
M = simplificarMatriz(M);

disp('Energia cinetica total K:'); %[output:911ec97e]
disp(K_total); %[output:9a910fd8]
fprintf('\nLaTeX de K:\n%s\n', latex(K_total)); %[output:44d3d894]

fprintf('\nMatriz de inercia M(q):\n'); %[output:9fa4849b]
disp(M); %[output:5fd21027]
fprintf('\nLaTeX de M(q):\n%s\n', latex(M)); %[output:608a565d]

% Coeficientes de la forma cuadratica K = qdot^T*A*qdot
coef = {
    'q_1^2',       A(1,1);
    'q_2^2',       A(2,2);
    'q_3^2',       A(3,3);
    'q_2 q_3',   2*A(2,3)
    };

for k = 1:size(coef,1) %[output:group:0733b04f]
    fprintf('\nCoeficiente de %s:\n', coef{k,1}); %[output:5929d911] %[output:6804be4b] %[output:34b23e0f] %[output:854b957e]
    disp(coef{k,2}); %[output:8f618880] %[output:0ab08848] %[output:378cbabd] %[output:5a79efb7]
    fprintf('LaTeX: %s\n', latex(coef{k,2})); %[output:673bbe1a] %[output:3414c340] %[output:753e059e] %[output:2f890b48]

end %[output:group:0733b04f]
%[text] ## FUNCIONES AUXILIARES
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
%%
function A = simplificarMatriz(A)
% Simplifica elemento por elemento para evitar salidas innecesariamente largas.

for i = 1:size(A,1)
    for j = 1:size(A,2)
        A(i,j) = simplificar(A(i,j));
    end
end
end
%%
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
%%
function z = verificarCero(expr)
% Rutina separada para verificar identidades. Aqui SI se permite expandir,
% porque no se usa para presentar el resultado sino para comprobar igualdad.

z = simplify(expand(rewrite(expr,'sincos')), 'Steps', 500);
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:30006ce4]
%   data: {"dataType":"text","outputData":{"text":"omega_1 en {0}:\n","truncated":false}}
%---
%[output:77bea529]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\n0\\\\\n0\\\\\n{\\dot{q} }_1 \n\\end{array}\\right)"}}
%---
%[output:7bf32629]
%   data: {"dataType":"text","outputData":{"text":"v_c_1\n","truncated":false}}
%---
%[output:7bf08550]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\n0\\\\\n0\\\\\n0\n\\end{array}\\right)"}}
%---
%[output:4a9c2d46]
%   data: {"dataType":"text","outputData":{"text":"K1 lineal\n","truncated":false}}
%---
%[output:61f85703]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"0"}}
%---
%[output:1e7fea41]
%   data: {"dataType":"text","outputData":{"text":"K1 rotacional\n","truncated":false}}
%---
%[output:6dd8d755]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{I_{1,\\textrm{zz}} \\,{{\\dot{q} }_1 }^2 }{2}"}}
%---
%[output:435fef61]
%   data: {"dataType":"text","outputData":{"text":"K1 total\n","truncated":false}}
%---
%[output:54c612c2]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{I_{1,\\textrm{zz}} \\,{{\\dot{q} }_1 }^2 }{2}"}}
%---
%[output:5afcd529]
%   data: {"dataType":"text","outputData":{"text":"omega_2 en {0}:\n","truncated":false}}
%---
%[output:951a5433]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\n{\\dot{q} }_2 \\,\\sin \\left(q_1 \\right)\\\\\n-{\\dot{q} }_2 \\,\\cos \\left(q_1 \\right)\\\\\n0\n\\end{array}\\right)"}}
%---
%[output:21da127d]
%   data: {"dataType":"text","outputData":{"text":"omega_2 en frame solidario {2}:\n","truncated":false}}
%---
%[output:6a515821]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\n0\\\\\n0\\\\\n{\\dot{q} }_2 \n\\end{array}\\right)"}}
%---
%[output:7341bb72]
%   data: {"dataType":"text","outputData":{"text":"v_c_2\n","truncated":false}}
%---
%[output:0bab4640]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\n-l_{c,2} \\,{\\dot{q} }_2 \\,\\cos \\left(q_1 \\right)\\,\\cos \\left(q_2 \\right)\\\\\n-l_{c,2} \\,{\\dot{q} }_2 \\,\\cos \\left(q_2 \\right)\\,\\sin \\left(q_1 \\right)\\\\\n-l_{c,2} \\,{\\dot{q} }_2 \\,\\sin \\left(q_2 \\right)\n\\end{array}\\right)"}}
%---
%[output:67fa81d8]
%   data: {"dataType":"text","outputData":{"text":"K2 traslacional:\n","truncated":false}}
%---
%[output:4bfdbac8]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{{l_{c,2} }^2 \\,m_2 \\,{{\\dot{q} }_2 }^2 }{2}"}}
%---
%[output:61624bd3]
%   data: {"dataType":"text","outputData":{"text":"K2 rotacional:\n","truncated":false}}
%---
%[output:4a396191]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{I_{2,\\textrm{zz}} \\,{{\\dot{q} }_2 }^2 }{2}"}}
%---
%[output:65251f6f]
%   data: {"dataType":"text","outputData":{"text":"K2 total ordenado:\n","truncated":false}}
%---
%[output:1d5af13d]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{{{\\dot{q} }_2 }^2 \\,{\\left(m_2 \\,{l_{c,2} }^2 +I_{2,\\textrm{zz}} \\right)}}{2}"}}
%---
%[output:4b0871b6]
%   data: {"dataType":"text","outputData":{"text":"Chequeo K2 - K2_ref (debe ser 0):\n","truncated":false}}
%---
%[output:8339c208]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"-\\frac{m_2 \\,{l_{c,2} }^2 \\,{{\\dot{q} }_1 }^2 \\,{\\sin \\left(q_2 \\right)}^2 }{2}-\\frac{I_{2,\\textrm{xx}} \\,{{\\dot{q} }_1 }^2 \\,{\\cos \\left(q_2 \\right)}^2 }{2}+I_{2,\\textrm{xy}} \\,{{\\dot{q} }_1 }^2 \\,\\cos \\left(q_2 \\right)\\,\\sin \\left(q_2 \\right)-\\frac{I_{2,\\textrm{yy}} \\,{{\\dot{q} }_1 }^2 \\,{\\sin \\left(q_2 \\right)}^2 }{2}-I_{2,\\textrm{xz}} \\,{\\dot{q} }_2 \\,{\\dot{q} }_1 \\,\\cos \\left(q_2 \\right)+I_{2,\\textrm{yz}} \\,{\\dot{q} }_2 \\,{\\dot{q} }_1 \\,\\sin \\left(q_2 \\right)"}}
%---
%[output:2848e350]
%   data: {"dataType":"text","outputData":{"text":"omega_3 en {0}:\n","truncated":false}}
%---
%[output:454c30c4]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\n\\sin \\left(q_1 \\right)\\,{\\left({\\dot{q} }_2 +{\\dot{q} }_3 \\right)}\\\\\n-\\cos \\left(q_1 \\right)\\,{\\left({\\dot{q} }_2 +{\\dot{q} }_3 \\right)}\\\\\n0\n\\end{array}\\right)"}}
%---
%[output:76dabcbd]
%   data: {"dataType":"text","outputData":{"text":"omega_3 en frame solidario {3}:\n","truncated":false}}
%---
%[output:73ef5ed9]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\n0\\\\\n0\\\\\n{\\dot{q} }_2 +{\\dot{q} }_3 \n\\end{array}\\right)"}}
%---
%[output:733a687c]
%   data: {"dataType":"text","outputData":{"text":"v_{O3}:\n","truncated":false}}
%---
%[output:5d4d816b]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\n-L_2 \\,{\\dot{q} }_2 \\,\\cos \\left(q_1 \\right)\\,\\cos \\left(q_2 \\right)\\\\\n-L_2 \\,{\\dot{q} }_2 \\,\\cos \\left(q_2 \\right)\\,\\sin \\left(q_1 \\right)\\\\\n-L_2 \\,{\\dot{q} }_2 \\,\\sin \\left(q_2 \\right)\n\\end{array}\\right)"}}
%---
%[output:359c9dbf]
%   data: {"dataType":"text","outputData":{"text":"v_c_3:\n","truncated":false}}
%---
%[output:9252f871]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\n-\\cos \\left(q_1 \\right)\\,{\\left(l_{c,3} \\,{\\dot{q} }_2 \\,\\cos \\left(q_2 +q_3 \\right)+l_{c,3} \\,{\\dot{q} }_3 \\,\\cos \\left(q_2 +q_3 \\right)+L_2 \\,{\\dot{q} }_2 \\,\\cos \\left(q_2 \\right)\\right)}\\\\\n-\\sin \\left(q_1 \\right)\\,{\\left(l_{c,3} \\,{\\dot{q} }_2 \\,\\cos \\left(q_2 +q_3 \\right)+l_{c,3} \\,{\\dot{q} }_3 \\,\\cos \\left(q_2 +q_3 \\right)+L_2 \\,{\\dot{q} }_2 \\,\\cos \\left(q_2 \\right)\\right)}\\\\\n-L_2 \\,{\\dot{q} }_2 \\,\\sin \\left(q_2 \\right)-l_{c,3} \\,\\sin \\left(q_2 +q_3 \\right)\\,{\\left({\\dot{q} }_2 +{\\dot{q} }_3 \\right)}\n\\end{array}\\right)"}}
%---
%[output:69e274e4]
%   data: {"dataType":"text","outputData":{"text":"K3 traslacional:\n","truncated":false}}
%---
%[output:62745e20]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{m_3 \\,{\\left({L_2 }^2 \\,{{\\dot{q} }_2 }^2 +2\\,\\cos \\left(q_3 \\right)\\,L_2 \\,l_{c,3} \\,{{\\dot{q} }_2 }^2 +2\\,\\cos \\left(q_3 \\right)\\,L_2 \\,l_{c,3} \\,{\\dot{q} }_2 \\,{\\dot{q} }_3 +{l_{c,3} }^2 \\,{{\\dot{q} }_2 }^2 +2\\,{l_{c,3} }^2 \\,{\\dot{q} }_2 \\,{\\dot{q} }_3 +{l_{c,3} }^2 \\,{{\\dot{q} }_3 }^2 \\right)}}{2}"}}
%---
%[output:866b1500]
%   data: {"dataType":"text","outputData":{"text":"K3 rotacional:\n","truncated":false}}
%---
%[output:37dafd57]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{I_{3,\\textrm{zz}} \\,{{\\left({\\dot{q} }_2 +{\\dot{q} }_3 \\right)}}^2 }{2}"}}
%---
%[output:25aa1212]
%   data: {"dataType":"text","outputData":{"text":"K3 total ordenado:\n","truncated":false}}
%---
%[output:9b011b1f]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"{\\left(\\frac{m_3 \\,{L_2 }^2 }{2}+m_3 \\,\\cos \\left(q_3 \\right)\\,L_2 \\,l_{c,3} +\\frac{m_3 \\,{l_{c,3} }^2 }{2}+\\frac{I_{3,\\textrm{zz}} }{2}\\right)}\\,{{\\dot{q} }_2 }^2 +{\\left(m_3 \\,{l_{c,3} }^2 +L_2 \\,m_3 \\,\\cos \\left(q_3 \\right)\\,l_{c,3} +I_{3,\\textrm{zz}} \\right)}\\,{\\dot{q} }_2 \\,{\\dot{q} }_3 +{\\left(\\frac{m_3 \\,{l_{c,3} }^2 }{2}+\\frac{I_{3,\\textrm{zz}} }{2}\\right)}\\,{{\\dot{q} }_3 }^2"}}
%---
%[output:8a3e050c]
%   data: {"dataType":"text","outputData":{"text":"Chequeo K3 - K3_ref (debe ser 0):\n","truncated":false}}
%---
%[output:688de766]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\begin{array}{l}\n\\frac{I_{3,\\textrm{xx}} \\,{{\\dot{q} }_1 }^2 \\,\\sigma_1 }{2}-I_{3,\\textrm{xz}} \\,{\\dot{q} }_1 \\,{\\dot{q} }_2 -I_{3,\\textrm{xz}} \\,{\\dot{q} }_1 \\,{\\dot{q} }_3 -\\frac{I_{3,\\textrm{xx}} \\,{{\\dot{q} }_1 }^2 }{2}-\\frac{I_{3,\\textrm{yy}} \\,{{\\dot{q} }_1 }^2 \\,\\sigma_1 }{2}+\\frac{I_{3,\\textrm{xy}} \\,{{\\dot{q} }_1 }^2 \\,\\sin \\left(2\\,q_2 +2\\,q_3 \\right)}{2}+I_{3,\\textrm{yz}} \\,{\\dot{q} }_1 \\,{\\dot{q} }_2 \\,\\sin \\left(q_2 +q_3 \\right)+I_{3,\\textrm{yz}} \\,{\\dot{q} }_1 \\,{\\dot{q} }_3 \\,\\sin \\left(q_2 +q_3 \\right)+2\\,I_{3,\\textrm{xz}} \\,{\\dot{q} }_1 \\,{\\dot{q} }_2 \\,\\sigma_2 +2\\,I_{3,\\textrm{xz}} \\,{\\dot{q} }_1 \\,{\\dot{q} }_3 \\,\\sigma_2 -\\frac{{l_{c,3} }^2 \\,m_3 \\,{{\\dot{q} }_1 }^2 \\,\\sigma_1 }{2}-\\frac{{L_2 }^2 \\,m_3 \\,{{\\dot{q} }_1 }^2 \\,{\\sin \\left(q_2 \\right)}^2 }{2}-L_2 \\,l_{c,3} \\,m_3 \\,{{\\dot{q} }_1 }^2 \\,{\\sin \\left(q_2 +\\frac{q_3 }{2}\\right)}^2 +L_2 \\,l_{c,3} \\,m_3 \\,{{\\dot{q} }_1 }^2 \\,{\\sin \\left(\\frac{q_3 }{2}\\right)}^2 \\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 ={\\sin \\left(q_2 +q_3 \\right)}^2 \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 ={\\sin \\left(\\frac{q_2 }{2}+\\frac{q_3 }{2}\\right)}^2 \n\\end{array}"}}
%---
%[output:911ec97e]
%   data: {"dataType":"text","outputData":{"text":"Energia cinetica total K:\n","truncated":false}}
%---
%[output:9a910fd8]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\begin{array}{l}\n\\frac{I_{1,\\textrm{zz}} \\,{{\\dot{q} }_1 }^2 }{2}+{\\left(\\frac{m_3 \\,{L_2 }^2 }{2}+m_3 \\,\\cos \\left(q_3 \\right)\\,L_2 \\,l_{c,3} +\\frac{m_2 \\,{l_{c,2} }^2 }{2}+\\sigma_1 +\\frac{I_{2,\\textrm{zz}} }{2}+\\frac{I_{3,\\textrm{zz}} }{2}\\right)}\\,{{\\dot{q} }_2 }^2 +{\\left(m_3 \\,{l_{c,3} }^2 +L_2 \\,m_3 \\,\\cos \\left(q_3 \\right)\\,l_{c,3} +I_{3,\\textrm{zz}} \\right)}\\,{\\dot{q} }_2 \\,{\\dot{q} }_3 +{\\left(\\sigma_1 +\\frac{I_{3,\\textrm{zz}} }{2}\\right)}\\,{{\\dot{q} }_3 }^2 \\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\frac{m_3 \\,{l_{c,3} }^2 }{2}\n\\end{array}"}}
%---
%[output:44d3d894]
%   data: {"dataType":"text","outputData":{"text":"\nLaTeX de K:\n\\frac{I_{1,\\mathrm{zz}}\\,{\\dot{q}_{1}}^2}{2}+\\left(\\frac{m_{3}\\,{L_{2}}^2}{2}+m_{3}\\,\\cos\\left(q_{3}\\right)\\,L_{2}\\,l_{c,3}+\\frac{m_{2}\\,{l_{c,2}}^2}{2}+\\frac{m_{3}\\,{l_{c,3}}^2}{2}+\\frac{I_{2,\\mathrm{zz}}}{2}+\\frac{I_{3,\\mathrm{zz}}}{2}\\right)\\,{\\dot{q}_{2}}^2+\\left(m_{3}\\,{l_{c,3}}^2+L_{2}\\,m_{3}\\,\\cos\\left(q_{3}\\right)\\,l_{c,3}+I_{3,\\mathrm{zz}}\\right)\\,\\dot{q}_{2}\\,\\dot{q}_{3}+\\left(\\frac{m_{3}\\,{l_{c,3}}^2}{2}+\\frac{I_{3,\\mathrm{zz}}}{2}\\right)\\,{\\dot{q}_{3}}^2\n","truncated":false}}
%---
%[output:9fa4849b]
%   data: {"dataType":"text","outputData":{"text":"\nMatriz de inercia M(q):\n","truncated":false}}
%---
%[output:5fd21027]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{ccc}\nI_{1,\\textrm{zz}}  & 0 & 0\\\\\n0 & m_3 \\,{L_2 }^2 +2\\,m_3 \\,\\cos \\left(q_3 \\right)\\,L_2 \\,l_{c,3} +m_2 \\,{l_{c,2} }^2 +m_3 \\,{l_{c,3} }^2 +I_{2,\\textrm{zz}} +I_{3,\\textrm{zz}}  & m_3 \\,{l_{c,3} }^2 +L_2 \\,m_3 \\,\\cos \\left(q_3 \\right)\\,l_{c,3} +I_{3,\\textrm{zz}} \\\\\n0 & m_3 \\,{l_{c,3} }^2 +L_2 \\,m_3 \\,\\cos \\left(q_3 \\right)\\,l_{c,3} +I_{3,\\textrm{zz}}  & m_3 \\,{l_{c,3} }^2 +I_{3,\\textrm{zz}} \n\\end{array}\\right)"}}
%---
%[output:608a565d]
%   data: {"dataType":"text","outputData":{"text":"\nLaTeX de M(q):\n\\left(\\begin{array}{ccc} I_{1,\\mathrm{zz}} & 0 & 0\\\\ 0 & m_{3}\\,{L_{2}}^2+2\\,m_{3}\\,\\cos\\left(q_{3}\\right)\\,L_{2}\\,l_{c,3}+m_{2}\\,{l_{c,2}}^2+m_{3}\\,{l_{c,3}}^2+I_{2,\\mathrm{zz}}+I_{3,\\mathrm{zz}} & m_{3}\\,{l_{c,3}}^2+L_{2}\\,m_{3}\\,\\cos\\left(q_{3}\\right)\\,l_{c,3}+I_{3,\\mathrm{zz}}\\\\ 0 & m_{3}\\,{l_{c,3}}^2+L_{2}\\,m_{3}\\,\\cos\\left(q_{3}\\right)\\,l_{c,3}+I_{3,\\mathrm{zz}} & m_{3}\\,{l_{c,3}}^2+I_{3,\\mathrm{zz}} \\end{array}\\right)\n","truncated":false}}
%---
%[output:5929d911]
%   data: {"dataType":"text","outputData":{"text":"\nCoeficiente de q_1^2:\n","truncated":false}}
%---
%[output:8f618880]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{I_{1,\\textrm{zz}} }{2}"}}
%---
%[output:673bbe1a]
%   data: {"dataType":"text","outputData":{"text":"LaTeX: \\frac{I_{1,\\mathrm{zz}}}{2}\n","truncated":false}}
%---
%[output:6804be4b]
%   data: {"dataType":"text","outputData":{"text":"\nCoeficiente de q_2^2:\n","truncated":false}}
%---
%[output:0ab08848]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{m_3 \\,{L_2 }^2 }{2}+m_3 \\,\\cos \\left(q_3 \\right)\\,L_2 \\,l_{c,3} +\\frac{m_2 \\,{l_{c,2} }^2 }{2}+\\frac{m_3 \\,{l_{c,3} }^2 }{2}+\\frac{I_{2,\\textrm{zz}} }{2}+\\frac{I_{3,\\textrm{zz}} }{2}"}}
%---
%[output:3414c340]
%   data: {"dataType":"text","outputData":{"text":"LaTeX: \\frac{m_{3}\\,{L_{2}}^2}{2}+m_{3}\\,\\cos\\left(q_{3}\\right)\\,L_{2}\\,l_{c,3}+\\frac{m_{2}\\,{l_{c,2}}^2}{2}+\\frac{m_{3}\\,{l_{c,3}}^2}{2}+\\frac{I_{2,\\mathrm{zz}}}{2}+\\frac{I_{3,\\mathrm{zz}}}{2}\n","truncated":false}}
%---
%[output:34b23e0f]
%   data: {"dataType":"text","outputData":{"text":"\nCoeficiente de q_3^2:\n","truncated":false}}
%---
%[output:378cbabd]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{m_3 \\,{l_{c,3} }^2 }{2}+\\frac{I_{3,\\textrm{zz}} }{2}"}}
%---
%[output:753e059e]
%   data: {"dataType":"text","outputData":{"text":"LaTeX: \\frac{m_{3}\\,{l_{c,3}}^2}{2}+\\frac{I_{3,\\mathrm{zz}}}{2}\n","truncated":false}}
%---
%[output:854b957e]
%   data: {"dataType":"text","outputData":{"text":"\nCoeficiente de q_2 q_3:\n","truncated":false}}
%---
%[output:5a79efb7]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"m_3 \\,{l_{c,3} }^2 +L_2 \\,m_3 \\,\\cos \\left(q_3 \\right)\\,l_{c,3} +I_{3,\\textrm{zz}}"}}
%---
%[output:2f890b48]
%   data: {"dataType":"text","outputData":{"text":"LaTeX: m_{3}\\,{l_{c,3}}^2+L_{2}\\,m_{3}\\,\\cos\\left(q_{3}\\right)\\,l_{c,3}+I_{3,\\mathrm{zz}}\n","truncated":false}}
%---
