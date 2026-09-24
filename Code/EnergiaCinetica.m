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

disp('omega_1 en {0}:'); disp(w(:,1)); %[output:2db6dc2c] %[output:8dea2fb7]
disp('v_c_1'); disp(vc(:,1)); %[output:51725f39] %[output:3c7f7e7f]
disp('K1 lineal'); disp(K(1,1)); %[output:74f03d23] %[output:306771cb]
disp('K1 rotacional'); disp(K(2,1)); %[output:759144e8] %[output:0782b9f1]
disp('K1 total'); disp(K(3,1)); %[output:76afcc6d] %[output:6ed73920]
%%
%[text] ## Energía cinética $K\_2$
% Velocidad angular relativa de 2 respecto de 1, expresada en {0}
w21 = q_dot_2 * [sin(q1); -cos(q1); 0];
w0(:,2) = simplificar(w0(:,1) + w21);

% Sistema solidario al eslabon 2:
% x2: longitudinal al eslabon
% z2: eje de la articulacion 2
% y2 = z2 x x2 para cerrar terna dextrógira
x2 = [...
    -cos(q1)*sin(q2);
    -sin(q1)*sin(q2);
    cos(q2)];

z2 = [...
    sin(q1);
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

                
% % Simplificación y reordenamiento de términos
% q_dot = [q_dot_1; q_dot_2];
% [best, A2] = simplifyCompact(K(3,2), 200);
% [best, ~] = quadraticCollect(best, q_dot);
% K(3,2) = best;

% Ordenar K2 segun velocidades articulares
q_dot_2vec = [q_dot_1; q_dot_2];
[K(3,2), A2] = quadraticCollect(K(3,2), q_dot_2vec);

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

%[text] $\\omega\_2$ en el sistema $\\{S\_0\\}$
disp(w0(:,2)); %[output:6e2146ee]
%[text] $\\omega\_2$ en el sistema solidario $\\{S\_2\\}$
disp(wb(:,2)); %[output:8f234988]
%[text] $v\_{C\_2}$
disp(vc(:,2)); %[output:8dbafc72]
%[text] $K\_2$ traslacional
disp(K(1,2)); %[output:8beec4c1]
%[text] $K\_2$ rotacional
disp(K(2,2)); %[output:451d162d]
%[text] $K\_2$ total ordenado
disp(K(3,2)); disp(latex(K(3,2))); %[output:4aa3cad5] %[output:79822bdf]
%[text] Chequeo $K\_2 - K\_{2,\\text{ref}}$
disp(check_K2); %[output:2de62887]
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

%[text] $\\omega\_3$ en el sistema $\\{S\_0\\}$
disp(w0(:,3)); %[output:13b5f5aa]
%[text] $\\omega\_3$ en el sistema solidario $\\{S\_3\\}$
disp(wb(:,3)); %[output:9b2e8203]
%[text] $v\_{O\_2}$
disp(v_O3); %[output:40886320]
%[text] $v\_{C\_2}$
disp(vc(:,3)); %[output:65d806fa]
%[text] $K\_3$ traslacional
disp(K(1,3)); %[output:15a139b8]
%[text] $K\_3$ rotacional
disp(K(2,3)); %[output:073a6037]
%[text] $K\_3$ total ordenado
disp(K(3,3)); disp(latex(K(3,3))); %[output:7ae8189e] %[output:74708ae9]
%[text] Chequeo $K\_3 - K\_{3,\\text{ref}}$
disp(check_K3); %[output:0519236d]

%%
%[text] ## Energía Cinética Total $K$
q_dot = [q_dot_1; q_dot_2; q_dot_3];
K_total = simplificar(K(3,1) + K(3,2) + K(3,3));
[K_total, A] = quadraticCollect(K_total, q_dot);

% Matriz de inercia del robot: K = 1/2*qdot.'*M(q)*qdot
M = hessian(K_total, q_dot);
M = simplificarMatriz(M);
%[text] Energia cinética total $K\_\\text{total}$
disp(K_total); disp(latex(K_total)); %[output:0627861c] %[output:28cc0848]
%[text] Matriz de inercia $\\mathbf{M}(\\underline{q})$
disp(M); disp(latex(M)); %[output:5fb14173] %[output:1b90dbce]

% Coeficientes de la forma cuadratica K = qdot^T*A*qdot
coef = {
    'q_dot_1^2',       A(1,1);
    'q_dot_2^2',       A(2,2);
    'q_dot_3^2',       A(3,3);
    'q_dot_1 q_dot_2', 2*A(1,2);
    'q_dot_1 q_dot_3', 2*A(1,3);
    'q_dot_2 q_dot_3', 2*A(2,3)
    };

for k = 1:size(coef,1) %[output:group:70742dd7]
    fprintf('\nCoeficiente de %s:\n', coef{k,1}); %[output:5da860fa] %[output:0928220f] %[output:6f412a29] %[output:0598359b] %[output:152394b3] %[output:9a81db9a]
    disp(coef{k,2}); %[output:5403eb8f] %[output:71673c66] %[output:41f70f41] %[output:9856851f] %[output:80cd910e] %[output:94681b93]
    fprintf('LaTeX: %s\n', latex(coef{k,2})); %[output:7b1e4103] %[output:474aaa01] %[output:0dd59a79] %[output:45ff9be2] %[output:3c8f4d56] %[output:3f456646]
end %[output:group:70742dd7]

%[text] ### Caso particular: ejes principales de inercia
% Si el frame solidario coincide con los ejes principales, los productos
% de inercia se anulan.
productos_inercia = [I_2_xy I_2_xz I_2_yz I_3_xy I_3_xz I_3_yz];
ceros_productos   = zeros(1,6);

K2_principal = simplificar(subs(K(3,2), productos_inercia, ceros_productos));
K3_principal = simplificar(subs(K(3,3), productos_inercia, ceros_productos));
K_principal  = simplificar(subs(K_total, productos_inercia, ceros_productos));
%[text] Caso de Ejes Principales
%[text] - $K\_2$ con productos de inercia nulos \
disp(K2_principal); %[output:059e82df]
%[text] - $K\_3$ con productos de inercia nulos \
disp(K3_principal); %[output:968358ea]
%[text] - $K\_{\\text{total}}$ con productos de inercia nulos \
disp(K_principal); %[output:249b1dc6]
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
%[output:2db6dc2c]
%   data: {"dataType":"text","outputData":{"text":"omega_1 en {0}:\n","truncated":false}}
%---
%[output:8dea2fb7]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\n0\\\\\n0\\\\\n{\\dot{q} }_1 \n\\end{array}\\right)"}}
%---
%[output:51725f39]
%   data: {"dataType":"text","outputData":{"text":"v_c_1\n","truncated":false}}
%---
%[output:3c7f7e7f]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\n0\\\\\n0\\\\\n0\n\\end{array}\\right)"}}
%---
%[output:74f03d23]
%   data: {"dataType":"text","outputData":{"text":"K1 lineal\n","truncated":false}}
%---
%[output:306771cb]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"0"}}
%---
%[output:759144e8]
%   data: {"dataType":"text","outputData":{"text":"K1 rotacional\n","truncated":false}}
%---
%[output:0782b9f1]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{I_{1,\\textrm{zz}} \\,{{\\dot{q} }_1 }^2 }{2}"}}
%---
%[output:76afcc6d]
%   data: {"dataType":"text","outputData":{"text":"K1 total\n","truncated":false}}
%---
%[output:6ed73920]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{I_{1,\\textrm{zz}} \\,{{\\dot{q} }_1 }^2 }{2}"}}
%---
%[output:6e2146ee]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\n{\\dot{q} }_2 \\,\\sin \\left(q_1 \\right)\\\\\n-{\\dot{q} }_2 \\,\\cos \\left(q_1 \\right)\\\\\n{\\dot{q} }_1 \n\\end{array}\\right)"}}
%---
%[output:8f234988]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\n{\\dot{q} }_1 \\,\\cos \\left(q_2 \\right)\\\\\n-{\\dot{q} }_1 \\,\\sin \\left(q_2 \\right)\\\\\n{\\dot{q} }_2 \n\\end{array}\\right)"}}
%---
%[output:8dbafc72]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\nl_{c,2} \\,{\\dot{q} }_1 \\,\\sin \\left(q_1 \\right)\\,\\sin \\left(q_2 \\right)-l_{c,2} \\,{\\dot{q} }_2 \\,\\cos \\left(q_1 \\right)\\,\\cos \\left(q_2 \\right)\\\\\n-l_{c,2} \\,{\\dot{q} }_1 \\,\\cos \\left(q_1 \\right)\\,\\sin \\left(q_2 \\right)-l_{c,2} \\,{\\dot{q} }_2 \\,\\cos \\left(q_2 \\right)\\,\\sin \\left(q_1 \\right)\\\\\n-l_{c,2} \\,{\\dot{q} }_2 \\,\\sin \\left(q_2 \\right)\n\\end{array}\\right)"}}
%---
%[output:8beec4c1]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{{l_{c,2} }^2 \\,m_2 \\,{\\left({{\\dot{q} }_1 }^2 \\,{\\sin \\left(q_2 \\right)}^2 +{{\\dot{q} }_2 }^2 \\right)}}{2}"}}
%---
%[output:451d162d]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{I_{2,\\textrm{zz}} \\,{{\\dot{q} }_2 }^2 }{2}-\\frac{I_{2,\\textrm{xx}} \\,{{\\dot{q} }_1 }^2 \\,{\\left({\\sin \\left(q_2 \\right)}^2 -1\\right)}}{2}-\\frac{I_{2,\\textrm{xy}} \\,{{\\dot{q} }_1 }^2 \\,\\sin \\left(2\\,q_2 \\right)}{2}+\\frac{I_{2,\\textrm{yy}} \\,{{\\dot{q} }_1 }^2 \\,{\\sin \\left(q_2 \\right)}^2 }{2}-I_{2,\\textrm{yz}} \\,{\\dot{q} }_1 \\,{\\dot{q} }_2 \\,\\sin \\left(q_2 \\right)-I_{2,\\textrm{xz}} \\,{\\dot{q} }_1 \\,{\\dot{q} }_2 \\,{\\left(2\\,{\\sin \\left(\\frac{q_2 }{2}\\right)}^2 -1\\right)}"}}
%---
%[output:4aa3cad5]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\begin{array}{l}\n{\\left(\\frac{I_{2,\\textrm{yy}} }{4}+\\frac{{l_{c,2} }^2 \\,m_2 }{4}-\\frac{I_{2,\\textrm{yy}} \\,\\sigma_1 }{4}+\\frac{I_{2,\\textrm{xx}} \\,{\\cos \\left(q_2 \\right)}^2 }{2}-I_{2,\\textrm{xy}} \\,\\cos \\left(q_2 \\right)\\,\\sin \\left(q_2 \\right)-\\frac{{l_{c,2} }^2 \\,m_2 \\,\\sigma_1 }{4}\\right)}\\,{{\\dot{q} }_1 }^2 +{\\left(I_{2,\\textrm{xz}} \\,\\cos \\left(q_2 \\right)-I_{2,\\textrm{yz}} \\,\\sin \\left(q_2 \\right)\\right)}\\,{\\dot{q} }_1 \\,{\\dot{q} }_2 +{\\left(\\frac{m_2 \\,{l_{c,2} }^2 }{2}+\\frac{I_{2,\\textrm{zz}} }{2}\\right)}\\,{{\\dot{q} }_2 }^2 \\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =2\\,{\\cos \\left(q_2 \\right)}^2 -1\n\\end{array}"}}
%---
%[output:79822bdf]
%   data: {"dataType":"text","outputData":{"text":"\\left(\\frac{I_{2,\\mathrm{yy}}}{4}+\\frac{{l_{c,2}}^2\\,m_{2}}{4}-\\frac{I_{2,\\mathrm{yy}}\\,\\left(2\\,{\\cos\\left(q_{2}\\right)}^2-1\\right)}{4}+\\frac{I_{2,\\mathrm{xx}}\\,{\\cos\\left(q_{2}\\right)}^2}{2}-I_{2,\\mathrm{xy}}\\,\\cos\\left(q_{2}\\right)\\,\\sin\\left(q_{2}\\right)-\\frac{{l_{c,2}}^2\\,m_{2}\\,\\left(2\\,{\\cos\\left(q_{2}\\right)}^2-1\\right)}{4}\\right)\\,{\\dot{q}_{1}}^2+\\left(I_{2,\\mathrm{xz}}\\,\\cos\\left(q_{2}\\right)-I_{2,\\mathrm{yz}}\\,\\sin\\left(q_{2}\\right)\\right)\\,\\dot{q}_{1}\\,\\dot{q}_{2}+\\left(\\frac{m_{2}\\,{l_{c,2}}^2}{2}+\\frac{I_{2,\\mathrm{zz}}}{2}\\right)\\,{\\dot{q}_{2}}^2\n","truncated":false}}
%---
%[output:2de62887]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"0"}}
%---
%[output:13b5f5aa]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\n\\sin \\left(q_1 \\right)\\,{\\left({\\dot{q} }_2 +{\\dot{q} }_3 \\right)}\\\\\n-\\cos \\left(q_1 \\right)\\,{\\left({\\dot{q} }_2 +{\\dot{q} }_3 \\right)}\\\\\n{\\dot{q} }_1 \n\\end{array}\\right)"}}
%---
%[output:9b2e8203]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\n{\\dot{q} }_1 \\,\\cos \\left(q_2 +q_3 \\right)\\\\\n-{\\dot{q} }_1 \\,\\sin \\left(q_2 +q_3 \\right)\\\\\n{\\dot{q} }_2 +{\\dot{q} }_3 \n\\end{array}\\right)"}}
%---
%[output:40886320]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\nL_2 \\,{\\dot{q} }_1 \\,\\sin \\left(q_1 \\right)\\,\\sin \\left(q_2 \\right)-L_2 \\,{\\dot{q} }_2 \\,\\cos \\left(q_1 \\right)\\,\\cos \\left(q_2 \\right)\\\\\n-L_2 \\,{\\dot{q} }_1 \\,\\cos \\left(q_1 \\right)\\,\\sin \\left(q_2 \\right)-L_2 \\,{\\dot{q} }_2 \\,\\cos \\left(q_2 \\right)\\,\\sin \\left(q_1 \\right)\\\\\n-L_2 \\,{\\dot{q} }_2 \\,\\sin \\left(q_2 \\right)\n\\end{array}\\right)"}}
%---
%[output:65d806fa]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\nL_2 \\,{\\dot{q} }_1 \\,\\sin \\left(q_1 \\right)\\,\\sin \\left(q_2 \\right)-L_2 \\,{\\dot{q} }_2 \\,\\cos \\left(q_1 \\right)\\,\\cos \\left(q_2 \\right)-l_{c,3} \\,{\\dot{q} }_2 \\,\\cos \\left(q_1 \\right)\\,\\cos \\left(q_2 \\right)\\,\\cos \\left(q_3 \\right)-l_{c,3} \\,{\\dot{q} }_3 \\,\\cos \\left(q_1 \\right)\\,\\cos \\left(q_2 \\right)\\,\\cos \\left(q_3 \\right)+l_{c,3} \\,{\\dot{q} }_1 \\,\\cos \\left(q_2 \\right)\\,\\sin \\left(q_1 \\right)\\,\\sin \\left(q_3 \\right)+l_{c,3} \\,{\\dot{q} }_1 \\,\\cos \\left(q_3 \\right)\\,\\sin \\left(q_1 \\right)\\,\\sin \\left(q_2 \\right)+l_{c,3} \\,{\\dot{q} }_2 \\,\\cos \\left(q_1 \\right)\\,\\sin \\left(q_2 \\right)\\,\\sin \\left(q_3 \\right)+l_{c,3} \\,{\\dot{q} }_3 \\,\\cos \\left(q_1 \\right)\\,\\sin \\left(q_2 \\right)\\,\\sin \\left(q_3 \\right)\\\\\nl_{c,3} \\,{\\dot{q} }_2 \\,\\sin \\left(q_1 \\right)\\,\\sin \\left(q_2 \\right)\\,\\sin \\left(q_3 \\right)-L_2 \\,{\\dot{q} }_2 \\,\\cos \\left(q_2 \\right)\\,\\sin \\left(q_1 \\right)-l_{c,3} \\,{\\dot{q} }_1 \\,\\cos \\left(q_1 \\right)\\,\\cos \\left(q_2 \\right)\\,\\sin \\left(q_3 \\right)-l_{c,3} \\,{\\dot{q} }_1 \\,\\cos \\left(q_1 \\right)\\,\\cos \\left(q_3 \\right)\\,\\sin \\left(q_2 \\right)-l_{c,3} \\,{\\dot{q} }_2 \\,\\cos \\left(q_2 \\right)\\,\\cos \\left(q_3 \\right)\\,\\sin \\left(q_1 \\right)-l_{c,3} \\,{\\dot{q} }_3 \\,\\cos \\left(q_2 \\right)\\,\\cos \\left(q_3 \\right)\\,\\sin \\left(q_1 \\right)-L_2 \\,{\\dot{q} }_1 \\,\\cos \\left(q_1 \\right)\\,\\sin \\left(q_2 \\right)+l_{c,3} \\,{\\dot{q} }_3 \\,\\sin \\left(q_1 \\right)\\,\\sin \\left(q_2 \\right)\\,\\sin \\left(q_3 \\right)\\\\\n-L_2 \\,{\\dot{q} }_2 \\,\\sin \\left(q_2 \\right)-l_{c,3} \\,\\sin \\left(q_2 +q_3 \\right)\\,{\\left({\\dot{q} }_2 +{\\dot{q} }_3 \\right)}\n\\end{array}\\right)"}}
%---
%[output:15a139b8]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{m_3 \\,{\\left({L_2 }^2 \\,{{\\dot{q} }_1 }^2 +2\\,{L_2 }^2 \\,{{\\dot{q} }_2 }^2 +{l_{c,3} }^2 \\,{{\\dot{q} }_1 }^2 +2\\,{l_{c,3} }^2 \\,{{\\dot{q} }_2 }^2 +2\\,{l_{c,3} }^2 \\,{{\\dot{q} }_3 }^2 +4\\,{l_{c,3} }^2 \\,{\\dot{q} }_2 \\,{\\dot{q} }_3 -{L_2 }^2 \\,{{\\dot{q} }_1 }^2 \\,\\cos \\left(2\\,q_2 \\right)-{l_{c,3} }^2 \\,{{\\dot{q} }_1 }^2 \\,\\cos \\left(2\\,q_2 +2\\,q_3 \\right)+2\\,L_2 \\,l_{c,3} \\,{{\\dot{q} }_1 }^2 \\,\\cos \\left(q_3 \\right)+4\\,L_2 \\,l_{c,3} \\,{{\\dot{q} }_2 }^2 \\,\\cos \\left(q_3 \\right)-2\\,L_2 \\,l_{c,3} \\,{{\\dot{q} }_1 }^2 \\,\\cos \\left(2\\,q_2 +q_3 \\right)+4\\,L_2 \\,l_{c,3} \\,{\\dot{q} }_2 \\,{\\dot{q} }_3 \\,\\cos \\left(q_3 \\right)\\right)}}{4}"}}
%---
%[output:073a6037]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{{\\left({\\dot{q} }_2 +{\\dot{q} }_3 \\right)}\\,{\\left(I_{3,\\textrm{zz}} \\,{\\left({\\dot{q} }_2 +{\\dot{q} }_3 \\right)}+I_{3,\\textrm{xz}} \\,{\\dot{q} }_1 \\,\\cos \\left(q_2 +q_3 \\right)-I_{3,\\textrm{yz}} \\,{\\dot{q} }_1 \\,\\sin \\left(q_2 +q_3 \\right)\\right)}}{2}+\\frac{{\\dot{q} }_1 \\,\\cos \\left(q_2 +q_3 \\right)\\,{\\left(I_{3,\\textrm{xz}} \\,{\\left({\\dot{q} }_2 +{\\dot{q} }_3 \\right)}+I_{3,\\textrm{xx}} \\,{\\dot{q} }_1 \\,\\cos \\left(q_2 +q_3 \\right)-I_{3,\\textrm{xy}} \\,{\\dot{q} }_1 \\,\\sin \\left(q_2 +q_3 \\right)\\right)}}{2}-\\frac{{\\dot{q} }_1 \\,\\sin \\left(q_2 +q_3 \\right)\\,{\\left(I_{3,\\textrm{yz}} \\,{\\left({\\dot{q} }_2 +{\\dot{q} }_3 \\right)}+I_{3,\\textrm{xy}} \\,{\\dot{q} }_1 \\,\\cos \\left(q_2 +q_3 \\right)-I_{3,\\textrm{yy}} \\,{\\dot{q} }_1 \\,\\sin \\left(q_2 +q_3 \\right)\\right)}}{2}"}}
%---
%[output:7ae8189e]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\begin{array}{l}\n{\\left(\\frac{I_{3,\\textrm{xx}} }{4}+\\frac{I_{3,\\textrm{yy}} }{4}+\\frac{I_{3,\\textrm{xx}} \\,\\cos \\left(\\sigma_2 \\right)}{4}-\\frac{I_{3,\\textrm{yy}} \\,\\cos \\left(\\sigma_2 \\right)}{4}-\\frac{I_{3,\\textrm{xy}} \\,\\sin \\left(\\sigma_2 \\right)}{2}+\\frac{{L_2 }^2 \\,m_3 }{4}+\\frac{{l_{c,3} }^2 \\,m_3 }{4}-\\frac{{L_2 }^2 \\,m_3 \\,\\cos \\left(2\\,q_2 \\right)}{4}-\\frac{{l_{c,3} }^2 \\,m_3 \\,\\cos \\left(\\sigma_2 \\right)}{4}+\\frac{L_2 \\,l_{c,3} \\,m_3 \\,\\cos \\left(q_3 \\right)}{2}-\\frac{L_2 \\,l_{c,3} \\,m_3 \\,\\cos \\left(2\\,q_2 +q_3 \\right)}{2}\\right)}\\,{{\\dot{q} }_1 }^2 +\\sigma_1 \\,{\\dot{q} }_1 \\,{\\dot{q} }_2 +\\sigma_1 \\,{\\dot{q} }_1 \\,{\\dot{q} }_3 +{\\left(\\frac{I_{3,\\textrm{zz}} }{2}+\\frac{m_3 \\,{\\left(4\\,{L_2 }^2 +8\\,\\cos \\left(q_3 \\right)\\,L_2 \\,l_{c,3} +4\\,{l_{c,3} }^2 \\right)}}{8}\\right)}\\,{{\\dot{q} }_2 }^2 +{\\left(I_{3,\\textrm{zz}} +\\frac{m_3 \\,{\\left(4\\,{l_{c,3} }^2 +4\\,L_2 \\,\\cos \\left(q_3 \\right)\\,l_{c,3} \\right)}}{4}\\right)}\\,{\\dot{q} }_2 \\,{\\dot{q} }_3 +{\\left(\\frac{m_3 \\,{l_{c,3} }^2 }{2}+\\frac{I_{3,\\textrm{zz}} }{2}\\right)}\\,{{\\dot{q} }_3 }^2 \\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =I_{3,\\textrm{xz}} \\,\\cos \\left(q_2 +q_3 \\right)-I_{3,\\textrm{yz}} \\,\\sin \\left(q_2 +q_3 \\right)\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 =2\\,q_2 +2\\,q_3 \n\\end{array}"}}
%---
%[output:74708ae9]
%   data: {"dataType":"text","outputData":{"text":"\\left(\\frac{I_{3,\\mathrm{xx}}}{4}+\\frac{I_{3,\\mathrm{yy}}}{4}+\\frac{I_{3,\\mathrm{xx}}\\,\\cos\\left(2\\,q_{2}+2\\,q_{3}\\right)}{4}-\\frac{I_{3,\\mathrm{yy}}\\,\\cos\\left(2\\,q_{2}+2\\,q_{3}\\right)}{4}-\\frac{I_{3,\\mathrm{xy}}\\,\\sin\\left(2\\,q_{2}+2\\,q_{3}\\right)}{2}+\\frac{{L_{2}}^2\\,m_{3}}{4}+\\frac{{l_{c,3}}^2\\,m_{3}}{4}-\\frac{{L_{2}}^2\\,m_{3}\\,\\cos\\left(2\\,q_{2}\\right)}{4}-\\frac{{l_{c,3}}^2\\,m_{3}\\,\\cos\\left(2\\,q_{2}+2\\,q_{3}\\right)}{4}+\\frac{L_{2}\\,l_{c,3}\\,m_{3}\\,\\cos\\left(q_{3}\\right)}{2}-\\frac{L_{2}\\,l_{c,3}\\,m_{3}\\,\\cos\\left(2\\,q_{2}+q_{3}\\right)}{2}\\right)\\,{\\dot{q}_{1}}^2+\\left(I_{3,\\mathrm{xz}}\\,\\cos\\left(q_{2}+q_{3}\\right)-I_{3,\\mathrm{yz}}\\,\\sin\\left(q_{2}+q_{3}\\right)\\right)\\,\\dot{q}_{1}\\,\\dot{q}_{2}+\\left(I_{3,\\mathrm{xz}}\\,\\cos\\left(q_{2}+q_{3}\\right)-I_{3,\\mathrm{yz}}\\,\\sin\\left(q_{2}+q_{3}\\right)\\right)\\,\\dot{q}_{1}\\,\\dot{q}_{3}+\\left(\\frac{I_{3,\\mathrm{zz}}}{2}+\\frac{m_{3}\\,\\left(4\\,{L_{2}}^2+8\\,\\cos\\left(q_{3}\\right)\\,L_{2}\\,l_{c,3}+4\\,{l_{c,3}}^2\\right)}{8}\\right)\\,{\\dot{q}_{2}}^2+\\left(I_{3,\\mathrm{zz}}+\\frac{m_{3}\\,\\left(4\\,{l_{c,3}}^2+4\\,L_{2}\\,\\cos\\left(q_{3}\\right)\\,l_{c,3}\\right)}{4}\\right)\\,\\dot{q}_{2}\\,\\dot{q}_{3}+\\left(\\frac{m_{3}\\,{l_{c,3}}^2}{2}+\\frac{I_{3,\\mathrm{zz}}}{2}\\right)\\,{\\dot{q}_{3}}^2\n","truncated":false}}
%---
%[output:0519236d]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"0"}}
%---
%[output:0627861c]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\begin{array}{l}\n\\frac{I_{3,\\textrm{xx}} \\,{{\\dot{q} }_1 }^2 }{4}+\\frac{I_{2,\\textrm{yy}} \\,{{\\dot{q} }_1 }^2 }{4}+\\frac{I_{1,\\textrm{zz}} \\,{{\\dot{q} }_1 }^2 }{2}+\\frac{I_{3,\\textrm{yy}} \\,{{\\dot{q} }_1 }^2 }{4}+\\frac{I_{2,\\textrm{zz}} \\,{{\\dot{q} }_2 }^2 }{2}+\\frac{I_{3,\\textrm{zz}} \\,{{\\dot{q} }_2 }^2 }{2}+\\frac{I_{3,\\textrm{zz}} \\,{{\\dot{q} }_3 }^2 }{2}+\\frac{{L_2 }^2 \\,m_3 \\,{{\\dot{q} }_1 }^2 }{4}+\\frac{{L_2 }^2 \\,m_3 \\,{{\\dot{q} }_2 }^2 }{2}+\\frac{{l_{c,2} }^2 \\,m_2 \\,{{\\dot{q} }_1 }^2 }{4}+\\frac{{l_{c,2} }^2 \\,m_2 \\,{{\\dot{q} }_2 }^2 }{2}+\\frac{{l_{c,3} }^2 \\,m_3 \\,{{\\dot{q} }_1 }^2 }{4}+\\frac{{l_{c,3} }^2 \\,m_3 \\,{{\\dot{q} }_2 }^2 }{2}+\\frac{{l_{c,3} }^2 \\,m_3 \\,{{\\dot{q} }_3 }^2 }{2}+I_{3,\\textrm{zz}} \\,{\\dot{q} }_2 \\,{\\dot{q} }_3 -\\frac{I_{2,\\textrm{yy}} \\,{{\\dot{q} }_1 }^2 \\,\\cos \\left(2\\,q_2 \\right)}{4}+\\frac{I_{2,\\textrm{xx}} \\,{{\\dot{q} }_1 }^2 \\,{\\cos \\left(q_2 \\right)}^2 }{2}+\\frac{I_{3,\\textrm{xx}} \\,{{\\dot{q} }_1 }^2 \\,\\cos \\left(\\sigma_1 \\right)}{4}-\\frac{I_{3,\\textrm{yy}} \\,{{\\dot{q} }_1 }^2 \\,\\cos \\left(\\sigma_1 \\right)}{4}-\\frac{I_{3,\\textrm{xy}} \\,{{\\dot{q} }_1 }^2 \\,\\sin \\left(\\sigma_1 \\right)}{2}+I_{3,\\textrm{xz}} \\,{\\dot{q} }_1 \\,{\\dot{q} }_2 \\,\\cos \\left(q_2 +q_3 \\right)+I_{3,\\textrm{xz}} \\,{\\dot{q} }_1 \\,{\\dot{q} }_3 \\,\\cos \\left(q_2 +q_3 \\right)-I_{2,\\textrm{xy}} \\,{{\\dot{q} }_1 }^2 \\,\\cos \\left(q_2 \\right)\\,\\sin \\left(q_2 \\right)-I_{3,\\textrm{yz}} \\,{\\dot{q} }_1 \\,{\\dot{q} }_2 \\,\\sin \\left(q_2 +q_3 \\right)-I_{3,\\textrm{yz}} \\,{\\dot{q} }_1 \\,{\\dot{q} }_3 \\,\\sin \\left(q_2 +q_3 \\right)+I_{2,\\textrm{xz}} \\,{\\dot{q} }_1 \\,{\\dot{q} }_2 \\,\\cos \\left(q_2 \\right)-I_{2,\\textrm{yz}} \\,{\\dot{q} }_1 \\,{\\dot{q} }_2 \\,\\sin \\left(q_2 \\right)+{l_{c,3} }^2 \\,m_3 \\,{\\dot{q} }_2 \\,{\\dot{q} }_3 -\\frac{{L_2 }^2 \\,m_3 \\,{{\\dot{q} }_1 }^2 \\,\\cos \\left(2\\,q_2 \\right)}{4}-\\frac{{l_{c,2} }^2 \\,m_2 \\,{{\\dot{q} }_1 }^2 \\,\\cos \\left(2\\,q_2 \\right)}{4}-\\frac{{l_{c,3} }^2 \\,m_3 \\,{{\\dot{q} }_1 }^2 \\,\\cos \\left(\\sigma_1 \\right)}{4}+\\frac{L_2 \\,l_{c,3} \\,m_3 \\,{{\\dot{q} }_1 }^2 \\,\\cos \\left(q_3 \\right)}{2}+L_2 \\,l_{c,3} \\,m_3 \\,{{\\dot{q} }_2 }^2 \\,\\cos \\left(q_3 \\right)-\\frac{L_2 \\,l_{c,3} \\,m_3 \\,{{\\dot{q} }_1 }^2 \\,\\cos \\left(2\\,q_2 +q_3 \\right)}{2}+L_2 \\,l_{c,3} \\,m_3 \\,{\\dot{q} }_2 \\,{\\dot{q} }_3 \\,\\cos \\left(q_3 \\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =2\\,q_2 +2\\,q_3 \n\\end{array}"}}
%---
%[output:28cc0848]
%   data: {"dataType":"text","outputData":{"text":"\\frac{I_{3,\\mathrm{xx}}\\,{\\dot{q}_{1}}^2}{4}+\\frac{I_{2,\\mathrm{yy}}\\,{\\dot{q}_{1}}^2}{4}+\\frac{I_{1,\\mathrm{zz}}\\,{\\dot{q}_{1}}^2}{2}+\\frac{I_{3,\\mathrm{yy}}\\,{\\dot{q}_{1}}^2}{4}+\\frac{I_{2,\\mathrm{zz}}\\,{\\dot{q}_{2}}^2}{2}+\\frac{I_{3,\\mathrm{zz}}\\,{\\dot{q}_{2}}^2}{2}+\\frac{I_{3,\\mathrm{zz}}\\,{\\dot{q}_{3}}^2}{2}+\\frac{{L_{2}}^2\\,m_{3}\\,{\\dot{q}_{1}}^2}{4}+\\frac{{L_{2}}^2\\,m_{3}\\,{\\dot{q}_{2}}^2}{2}+\\frac{{l_{c,2}}^2\\,m_{2}\\,{\\dot{q}_{1}}^2}{4}+\\frac{{l_{c,2}}^2\\,m_{2}\\,{\\dot{q}_{2}}^2}{2}+\\frac{{l_{c,3}}^2\\,m_{3}\\,{\\dot{q}_{1}}^2}{4}+\\frac{{l_{c,3}}^2\\,m_{3}\\,{\\dot{q}_{2}}^2}{2}+\\frac{{l_{c,3}}^2\\,m_{3}\\,{\\dot{q}_{3}}^2}{2}+I_{3,\\mathrm{zz}}\\,\\dot{q}_{2}\\,\\dot{q}_{3}-\\frac{I_{2,\\mathrm{yy}}\\,{\\dot{q}_{1}}^2\\,\\cos\\left(2\\,q_{2}\\right)}{4}+\\frac{I_{2,\\mathrm{xx}}\\,{\\dot{q}_{1}}^2\\,{\\cos\\left(q_{2}\\right)}^2}{2}+\\frac{I_{3,\\mathrm{xx}}\\,{\\dot{q}_{1}}^2\\,\\cos\\left(2\\,q_{2}+2\\,q_{3}\\right)}{4}-\\frac{I_{3,\\mathrm{yy}}\\,{\\dot{q}_{1}}^2\\,\\cos\\left(2\\,q_{2}+2\\,q_{3}\\right)}{4}-\\frac{I_{3,\\mathrm{xy}}\\,{\\dot{q}_{1}}^2\\,\\sin\\left(2\\,q_{2}+2\\,q_{3}\\right)}{2}+I_{3,\\mathrm{xz}}\\,\\dot{q}_{1}\\,\\dot{q}_{2}\\,\\cos\\left(q_{2}+q_{3}\\right)+I_{3,\\mathrm{xz}}\\,\\dot{q}_{1}\\,\\dot{q}_{3}\\,\\cos\\left(q_{2}+q_{3}\\right)-I_{2,\\mathrm{xy}}\\,{\\dot{q}_{1}}^2\\,\\cos\\left(q_{2}\\right)\\,\\sin\\left(q_{2}\\right)-I_{3,\\mathrm{yz}}\\,\\dot{q}_{1}\\,\\dot{q}_{2}\\,\\sin\\left(q_{2}+q_{3}\\right)-I_{3,\\mathrm{yz}}\\,\\dot{q}_{1}\\,\\dot{q}_{3}\\,\\sin\\left(q_{2}+q_{3}\\right)+I_{2,\\mathrm{xz}}\\,\\dot{q}_{1}\\,\\dot{q}_{2}\\,\\cos\\left(q_{2}\\right)-I_{2,\\mathrm{yz}}\\,\\dot{q}_{1}\\,\\dot{q}_{2}\\,\\sin\\left(q_{2}\\right)+{l_{c,3}}^2\\,m_{3}\\,\\dot{q}_{2}\\,\\dot{q}_{3}-\\frac{{L_{2}}^2\\,m_{3}\\,{\\dot{q}_{1}}^2\\,\\cos\\left(2\\,q_{2}\\right)}{4}-\\frac{{l_{c,2}}^2\\,m_{2}\\,{\\dot{q}_{1}}^2\\,\\cos\\left(2\\,q_{2}\\right)}{4}-\\frac{{l_{c,3}}^2\\,m_{3}\\,{\\dot{q}_{1}}^2\\,\\cos\\left(2\\,q_{2}+2\\,q_{3}\\right)}{4}+\\frac{L_{2}\\,l_{c,3}\\,m_{3}\\,{\\dot{q}_{1}}^2\\,\\cos\\left(q_{3}\\right)}{2}+L_{2}\\,l_{c,3}\\,m_{3}\\,{\\dot{q}_{2}}^2\\,\\cos\\left(q_{3}\\right)-\\frac{L_{2}\\,l_{c,3}\\,m_{3}\\,{\\dot{q}_{1}}^2\\,\\cos\\left(2\\,q_{2}+q_{3}\\right)}{2}+L_{2}\\,l_{c,3}\\,m_{3}\\,\\dot{q}_{2}\\,\\dot{q}_{3}\\,\\cos\\left(q_{3}\\right)\n","truncated":false}}
%---
%[output:5fb14173]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\begin{array}{l}\n\\left(\\begin{array}{ccc}\n\\frac{I_{3,\\textrm{xx}} }{2}+I_{1,\\textrm{zz}} +\\frac{I_{3,\\textrm{yy}} }{2}+\\frac{I_{3,\\textrm{xx}} \\,\\cos \\left(\\sigma_2 \\right)}{2}-\\frac{I_{3,\\textrm{yy}} \\,\\cos \\left(\\sigma_2 \\right)}{2}-I_{3,\\textrm{xy}} \\,\\sin \\left(\\sigma_2 \\right)+{L_2 }^2 \\,m_3 +{l_{c,2} }^2 \\,m_2 +\\frac{{l_{c,3} }^2 \\,m_3 }{2}-I_{2,\\textrm{yy}} \\,{\\left(\\frac{\\cos \\left(2\\,q_2 \\right)}{2}-\\frac{1}{2}\\right)}+I_{2,\\textrm{xx}} \\,{\\cos \\left(q_2 \\right)}^2 -I_{2,\\textrm{xy}} \\,\\sin \\left(2\\,q_2 \\right)-{L_2 }^2 \\,m_3 \\,{\\cos \\left(q_2 \\right)}^2 -{l_{c,2} }^2 \\,m_2 \\,{\\cos \\left(q_2 \\right)}^2 -\\frac{{l_{c,3} }^2 \\,m_3 \\,\\cos \\left(\\sigma_2 \\right)}{2}+L_2 \\,l_{c,3} \\,m_3 \\,\\cos \\left(q_3 \\right)-L_2 \\,l_{c,3} \\,m_3 \\,\\cos \\left(2\\,q_2 +q_3 \\right) & \\sigma_3  & \\sigma_5 -\\sigma_4 \\\\\n\\sigma_3  & m_3 \\,{L_2 }^2 +2\\,m_3 \\,\\cos \\left(q_3 \\right)\\,L_2 \\,l_{c,3} +m_2 \\,{l_{c,2} }^2 +m_3 \\,{l_{c,3} }^2 +I_{2,\\textrm{zz}} +I_{3,\\textrm{zz}}  & \\sigma_1 \\\\\n\\sigma_5 -\\sigma_4  & \\sigma_1  & m_3 \\,{l_{c,3} }^2 +I_{3,\\textrm{zz}} \n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =m_3 \\,{l_{c,3} }^2 +L_2 \\,m_3 \\,\\cos \\left(q_3 \\right)\\,l_{c,3} +I_{3,\\textrm{zz}} \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 =2\\,q_2 +2\\,q_3 \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_3 =\\sigma_5 -\\sigma_4 +I_{2,\\textrm{xz}} \\,\\cos \\left(q_2 \\right)-I_{2,\\textrm{yz}} \\,\\sin \\left(q_2 \\right)\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_4 =I_{3,\\textrm{yz}} \\,\\sin \\left(q_2 +q_3 \\right)\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_5 =I_{3,\\textrm{xz}} \\,\\cos \\left(q_2 +q_3 \\right)\n\\end{array}"}}
%---
%[output:1b90dbce]
%   data: {"dataType":"text","outputData":{"text":"\\left(\\begin{array}{ccc} \\frac{I_{3,\\mathrm{xx}}}{2}+I_{1,\\mathrm{zz}}+\\frac{I_{3,\\mathrm{yy}}}{2}+\\frac{I_{3,\\mathrm{xx}}\\,\\cos\\left(2\\,q_{2}+2\\,q_{3}\\right)}{2}-\\frac{I_{3,\\mathrm{yy}}\\,\\cos\\left(2\\,q_{2}+2\\,q_{3}\\right)}{2}-I_{3,\\mathrm{xy}}\\,\\sin\\left(2\\,q_{2}+2\\,q_{3}\\right)+{L_{2}}^2\\,m_{3}+{l_{c,2}}^2\\,m_{2}+\\frac{{l_{c,3}}^2\\,m_{3}}{2}-I_{2,\\mathrm{yy}}\\,\\left(\\frac{\\cos\\left(2\\,q_{2}\\right)}{2}-\\frac{1}{2}\\right)+I_{2,\\mathrm{xx}}\\,{\\cos\\left(q_{2}\\right)}^2-I_{2,\\mathrm{xy}}\\,\\sin\\left(2\\,q_{2}\\right)-{L_{2}}^2\\,m_{3}\\,{\\cos\\left(q_{2}\\right)}^2-{l_{c,2}}^2\\,m_{2}\\,{\\cos\\left(q_{2}\\right)}^2-\\frac{{l_{c,3}}^2\\,m_{3}\\,\\cos\\left(2\\,q_{2}+2\\,q_{3}\\right)}{2}+L_{2}\\,l_{c,3}\\,m_{3}\\,\\cos\\left(q_{3}\\right)-L_{2}\\,l_{c,3}\\,m_{3}\\,\\cos\\left(2\\,q_{2}+q_{3}\\right) & I_{3,\\mathrm{xz}}\\,\\cos\\left(q_{2}+q_{3}\\right)-I_{3,\\mathrm{yz}}\\,\\sin\\left(q_{2}+q_{3}\\right)+I_{2,\\mathrm{xz}}\\,\\cos\\left(q_{2}\\right)-I_{2,\\mathrm{yz}}\\,\\sin\\left(q_{2}\\right) & I_{3,\\mathrm{xz}}\\,\\cos\\left(q_{2}+q_{3}\\right)-I_{3,\\mathrm{yz}}\\,\\sin\\left(q_{2}+q_{3}\\right)\\\\ I_{3,\\mathrm{xz}}\\,\\cos\\left(q_{2}+q_{3}\\right)-I_{3,\\mathrm{yz}}\\,\\sin\\left(q_{2}+q_{3}\\right)+I_{2,\\mathrm{xz}}\\,\\cos\\left(q_{2}\\right)-I_{2,\\mathrm{yz}}\\,\\sin\\left(q_{2}\\right) & m_{3}\\,{L_{2}}^2+2\\,m_{3}\\,\\cos\\left(q_{3}\\right)\\,L_{2}\\,l_{c,3}+m_{2}\\,{l_{c,2}}^2+m_{3}\\,{l_{c,3}}^2+I_{2,\\mathrm{zz}}+I_{3,\\mathrm{zz}} & m_{3}\\,{l_{c,3}}^2+L_{2}\\,m_{3}\\,\\cos\\left(q_{3}\\right)\\,l_{c,3}+I_{3,\\mathrm{zz}}\\\\ I_{3,\\mathrm{xz}}\\,\\cos\\left(q_{2}+q_{3}\\right)-I_{3,\\mathrm{yz}}\\,\\sin\\left(q_{2}+q_{3}\\right) & m_{3}\\,{l_{c,3}}^2+L_{2}\\,m_{3}\\,\\cos\\left(q_{3}\\right)\\,l_{c,3}+I_{3,\\mathrm{zz}} & m_{3}\\,{l_{c,3}}^2+I_{3,\\mathrm{zz}} \\end{array}\\right)\n","truncated":false}}
%---
%[output:5da860fa]
%   data: {"dataType":"text","outputData":{"text":"\nCoeficiente de q_dot_1^2:\n","truncated":false}}
%---
%[output:5403eb8f]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\begin{array}{l}\n\\frac{I_{3,\\textrm{xx}} }{4}+\\frac{I_{1,\\textrm{zz}} }{2}+\\frac{I_{3,\\textrm{yy}} }{4}+\\frac{I_{3,\\textrm{xx}} \\,\\cos \\left(\\sigma_1 \\right)}{4}-\\frac{I_{3,\\textrm{yy}} \\,\\cos \\left(\\sigma_1 \\right)}{4}-\\frac{I_{3,\\textrm{xy}} \\,\\sin \\left(\\sigma_1 \\right)}{2}+\\frac{{L_2 }^2 \\,m_3 }{2}+\\frac{{l_{c,2} }^2 \\,m_2 }{4}+\\frac{{l_{c,3} }^2 \\,m_3 }{4}-\\frac{I_{2,\\textrm{yy}} \\,{\\left(\\frac{\\cos \\left(2\\,q_2 \\right)}{2}-\\frac{1}{2}\\right)}}{2}+\\frac{I_{2,\\textrm{xx}} \\,{\\cos \\left(q_2 \\right)}^2 }{2}-\\frac{I_{2,\\textrm{xy}} \\,\\sin \\left(2\\,q_2 \\right)}{2}-\\frac{{L_2 }^2 \\,m_3 \\,{\\cos \\left(q_2 \\right)}^2 }{2}-\\frac{{l_{c,2} }^2 \\,m_2 \\,\\cos \\left(2\\,q_2 \\right)}{4}-\\frac{{l_{c,3} }^2 \\,m_3 \\,\\cos \\left(\\sigma_1 \\right)}{4}+\\frac{L_2 \\,l_{c,3} \\,m_3 \\,\\cos \\left(q_3 \\right)}{2}-\\frac{L_2 \\,l_{c,3} \\,m_3 \\,\\cos \\left(2\\,q_2 +q_3 \\right)}{2}\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =2\\,q_2 +2\\,q_3 \n\\end{array}"}}
%---
%[output:7b1e4103]
%   data: {"dataType":"text","outputData":{"text":"LaTeX: \\frac{I_{3,\\mathrm{xx}}}{4}+\\frac{I_{1,\\mathrm{zz}}}{2}+\\frac{I_{3,\\mathrm{yy}}}{4}+\\frac{I_{3,\\mathrm{xx}}\\,\\cos\\left(2\\,q_{2}+2\\,q_{3}\\right)}{4}-\\frac{I_{3,\\mathrm{yy}}\\,\\cos\\left(2\\,q_{2}+2\\,q_{3}\\right)}{4}-\\frac{I_{3,\\mathrm{xy}}\\,\\sin\\left(2\\,q_{2}+2\\,q_{3}\\right)}{2}+\\frac{{L_{2}}^2\\,m_{3}}{2}+\\frac{{l_{c,2}}^2\\,m_{2}}{4}+\\frac{{l_{c,3}}^2\\,m_{3}}{4}-\\frac{I_{2,\\mathrm{yy}}\\,\\left(\\frac{\\cos\\left(2\\,q_{2}\\right)}{2}-\\frac{1}{2}\\right)}{2}+\\frac{I_{2,\\mathrm{xx}}\\,{\\cos\\left(q_{2}\\right)}^2}{2}-\\frac{I_{2,\\mathrm{xy}}\\,\\sin\\left(2\\,q_{2}\\right)}{2}-\\frac{{L_{2}}^2\\,m_{3}\\,{\\cos\\left(q_{2}\\right)}^2}{2}-\\frac{{l_{c,2}}^2\\,m_{2}\\,\\cos\\left(2\\,q_{2}\\right)}{4}-\\frac{{l_{c,3}}^2\\,m_{3}\\,\\cos\\left(2\\,q_{2}+2\\,q_{3}\\right)}{4}+\\frac{L_{2}\\,l_{c,3}\\,m_{3}\\,\\cos\\left(q_{3}\\right)}{2}-\\frac{L_{2}\\,l_{c,3}\\,m_{3}\\,\\cos\\left(2\\,q_{2}+q_{3}\\right)}{2}\n","truncated":false}}
%---
%[output:0928220f]
%   data: {"dataType":"text","outputData":{"text":"\nCoeficiente de q_dot_2^2:\n","truncated":false}}
%---
%[output:71673c66]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{m_3 \\,{L_2 }^2 }{2}+m_3 \\,\\cos \\left(q_3 \\right)\\,L_2 \\,l_{c,3} +\\frac{m_2 \\,{l_{c,2} }^2 }{2}+\\frac{m_3 \\,{l_{c,3} }^2 }{2}+\\frac{I_{2,\\textrm{zz}} }{2}+\\frac{I_{3,\\textrm{zz}} }{2}"}}
%---
%[output:474aaa01]
%   data: {"dataType":"text","outputData":{"text":"LaTeX: \\frac{m_{3}\\,{L_{2}}^2}{2}+m_{3}\\,\\cos\\left(q_{3}\\right)\\,L_{2}\\,l_{c,3}+\\frac{m_{2}\\,{l_{c,2}}^2}{2}+\\frac{m_{3}\\,{l_{c,3}}^2}{2}+\\frac{I_{2,\\mathrm{zz}}}{2}+\\frac{I_{3,\\mathrm{zz}}}{2}\n","truncated":false}}
%---
%[output:6f412a29]
%   data: {"dataType":"text","outputData":{"text":"\nCoeficiente de q_dot_3^2:\n","truncated":false}}
%---
%[output:41f70f41]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{m_3 \\,{l_{c,3} }^2 }{2}+\\frac{I_{3,\\textrm{zz}} }{2}"}}
%---
%[output:0dd59a79]
%   data: {"dataType":"text","outputData":{"text":"LaTeX: \\frac{m_{3}\\,{l_{c,3}}^2}{2}+\\frac{I_{3,\\mathrm{zz}}}{2}\n","truncated":false}}
%---
%[output:0598359b]
%   data: {"dataType":"text","outputData":{"text":"\nCoeficiente de q_dot_1 q_dot_2:\n","truncated":false}}
%---
%[output:9856851f]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"I_{3,\\textrm{xz}} \\,\\cos \\left(q_2 +q_3 \\right)-I_{3,\\textrm{yz}} \\,\\sin \\left(q_2 +q_3 \\right)+I_{2,\\textrm{xz}} \\,\\cos \\left(q_2 \\right)-I_{2,\\textrm{yz}} \\,\\sin \\left(q_2 \\right)"}}
%---
%[output:45ff9be2]
%   data: {"dataType":"text","outputData":{"text":"LaTeX: I_{3,\\mathrm{xz}}\\,\\cos\\left(q_{2}+q_{3}\\right)-I_{3,\\mathrm{yz}}\\,\\sin\\left(q_{2}+q_{3}\\right)+I_{2,\\mathrm{xz}}\\,\\cos\\left(q_{2}\\right)-I_{2,\\mathrm{yz}}\\,\\sin\\left(q_{2}\\right)\n","truncated":false}}
%---
%[output:152394b3]
%   data: {"dataType":"text","outputData":{"text":"\nCoeficiente de q_dot_1 q_dot_3:\n","truncated":false}}
%---
%[output:80cd910e]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"I_{3,\\textrm{xz}} \\,\\cos \\left(q_2 +q_3 \\right)-I_{3,\\textrm{yz}} \\,\\sin \\left(q_2 +q_3 \\right)"}}
%---
%[output:3c8f4d56]
%   data: {"dataType":"text","outputData":{"text":"LaTeX: I_{3,\\mathrm{xz}}\\,\\cos\\left(q_{2}+q_{3}\\right)-I_{3,\\mathrm{yz}}\\,\\sin\\left(q_{2}+q_{3}\\right)\n","truncated":false}}
%---
%[output:9a81db9a]
%   data: {"dataType":"text","outputData":{"text":"\nCoeficiente de q_dot_2 q_dot_3:\n","truncated":false}}
%---
%[output:94681b93]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"m_3 \\,{l_{c,3} }^2 +L_2 \\,m_3 \\,\\cos \\left(q_3 \\right)\\,l_{c,3} +I_{3,\\textrm{zz}}"}}
%---
%[output:3f456646]
%   data: {"dataType":"text","outputData":{"text":"LaTeX: m_{3}\\,{l_{c,3}}^2+L_{2}\\,m_{3}\\,\\cos\\left(q_{3}\\right)\\,l_{c,3}+I_{3,\\mathrm{zz}}\n","truncated":false}}
%---
%[output:059e82df]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"{\\left(\\frac{I_{2,\\textrm{yy}} }{4}+\\frac{{l_{c,2} }^2 \\,m_2 }{4}+\\frac{I_{2,\\textrm{xx}} \\,{\\left(\\frac{\\cos \\left(2\\,q_2 \\right)}{2}+\\frac{1}{2}\\right)}}{2}-\\frac{I_{2,\\textrm{yy}} \\,\\cos \\left(2\\,q_2 \\right)}{4}-\\frac{{l_{c,2} }^2 \\,m_2 \\,\\cos \\left(2\\,q_2 \\right)}{4}\\right)}\\,{{\\dot{q} }_1 }^2 +{\\left(\\frac{m_2 \\,{l_{c,2} }^2 }{2}+\\frac{I_{2,\\textrm{zz}} }{2}\\right)}\\,{{\\dot{q} }_2 }^2"}}
%---
%[output:968358ea]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\begin{array}{l}\n{\\left(\\frac{I_{3,\\textrm{xx}} }{4}+\\frac{I_{3,\\textrm{yy}} }{4}+\\frac{I_{3,\\textrm{xx}} \\,\\sigma_1 }{4}-\\frac{I_{3,\\textrm{yy}} \\,\\sigma_1 }{4}+\\frac{{L_2 }^2 \\,m_3 }{4}+\\frac{{l_{c,3} }^2 \\,m_3 }{4}-\\frac{{L_2 }^2 \\,m_3 \\,\\cos \\left(2\\,q_2 \\right)}{4}-\\frac{{l_{c,3} }^2 \\,m_3 \\,\\sigma_1 }{4}+\\frac{L_2 \\,l_{c,3} \\,m_3 \\,\\cos \\left(q_3 \\right)}{2}-\\frac{L_2 \\,l_{c,3} \\,m_3 \\,\\cos \\left(2\\,q_2 +q_3 \\right)}{2}\\right)}\\,{{\\dot{q} }_1 }^2 +{\\left(\\frac{I_{3,\\textrm{zz}} }{2}+\\frac{m_3 \\,{\\left(4\\,{L_2 }^2 +8\\,\\cos \\left(q_3 \\right)\\,L_2 \\,l_{c,3} +4\\,{l_{c,3} }^2 \\right)}}{8}\\right)}\\,{{\\dot{q} }_2 }^2 +{\\left(I_{3,\\textrm{zz}} +\\frac{m_3 \\,{\\left(4\\,{l_{c,3} }^2 +4\\,L_2 \\,\\cos \\left(q_3 \\right)\\,l_{c,3} \\right)}}{4}\\right)}\\,{\\dot{q} }_2 \\,{\\dot{q} }_3 +{\\left(\\frac{m_3 \\,{l_{c,3} }^2 }{2}+\\frac{I_{3,\\textrm{zz}} }{2}\\right)}\\,{{\\dot{q} }_3 }^2 \\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\cos \\left(2\\,q_2 +2\\,q_3 \\right)\n\\end{array}"}}
%---
%[output:249b1dc6]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\begin{array}{l}\n\\frac{I_{2,\\textrm{xx}} \\,{{\\dot{q} }_1 }^2 }{2}+\\frac{I_{3,\\textrm{xx}} \\,{{\\dot{q} }_1 }^2 }{2}+\\frac{I_{1,\\textrm{zz}} \\,{{\\dot{q} }_1 }^2 }{2}+\\frac{I_{2,\\textrm{zz}} \\,{{\\dot{q} }_2 }^2 }{2}+\\frac{I_{3,\\textrm{zz}} \\,{{\\dot{q} }_2 }^2 }{2}+\\frac{I_{3,\\textrm{zz}} \\,{{\\dot{q} }_3 }^2 }{2}+\\frac{{L_2 }^2 \\,m_3 \\,{{\\dot{q} }_2 }^2 }{2}+\\frac{{l_{c,2} }^2 \\,m_2 \\,{{\\dot{q} }_2 }^2 }{2}+\\frac{{l_{c,3} }^2 \\,m_3 \\,{{\\dot{q} }_2 }^2 }{2}+\\frac{{l_{c,3} }^2 \\,m_3 \\,{{\\dot{q} }_3 }^2 }{2}+I_{3,\\textrm{zz}} \\,{\\dot{q} }_2 \\,{\\dot{q} }_3 -\\frac{I_{3,\\textrm{xx}} \\,{{\\dot{q} }_1 }^2 \\,\\sigma_2 }{2}+\\frac{I_{3,\\textrm{yy}} \\,{{\\dot{q} }_1 }^2 \\,\\sigma_2 }{2}-\\frac{I_{2,\\textrm{xx}} \\,{{\\dot{q} }_1 }^2 \\,{\\sin \\left(q_2 \\right)}^2 }{2}+\\frac{I_{2,\\textrm{yy}} \\,{{\\dot{q} }_1 }^2 \\,{\\sin \\left(q_2 \\right)}^2 }{2}+L_2 \\,l_{c,3} \\,m_3 \\,{{\\dot{q} }_2 }^2 +{l_{c,3} }^2 \\,m_3 \\,{\\dot{q} }_2 \\,{\\dot{q} }_3 +\\frac{{l_{c,3} }^2 \\,m_3 \\,{{\\dot{q} }_1 }^2 \\,\\sigma_2 }{2}+\\frac{{L_2 }^2 \\,m_3 \\,{{\\dot{q} }_1 }^2 \\,{\\sin \\left(q_2 \\right)}^2 }{2}+\\frac{{l_{c,2} }^2 \\,m_2 \\,{{\\dot{q} }_1 }^2 \\,{\\sin \\left(q_2 \\right)}^2 }{2}+L_2 \\,l_{c,3} \\,m_3 \\,{\\dot{q} }_2 \\,{\\dot{q} }_3 +L_2 \\,l_{c,3} \\,m_3 \\,{{\\dot{q} }_1 }^2 \\,{\\sin \\left(q_2 +\\frac{q_3 }{2}\\right)}^2 -L_2 \\,l_{c,3} \\,m_3 \\,{{\\dot{q} }_1 }^2 \\,\\sigma_1 -2\\,L_2 \\,l_{c,3} \\,m_3 \\,{{\\dot{q} }_2 }^2 \\,\\sigma_1 -2\\,L_2 \\,l_{c,3} \\,m_3 \\,{\\dot{q} }_2 \\,{\\dot{q} }_3 \\,\\sigma_1 \\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 ={\\sin \\left(\\frac{q_3 }{2}\\right)}^2 \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 ={\\sin \\left(q_2 +q_3 \\right)}^2 \n\\end{array}"}}
%---
