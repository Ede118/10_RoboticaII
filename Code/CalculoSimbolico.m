%[text] # PRODUCTO PUNTO Y PRODUCTO VECTORIAL SIMBÓLICO EN $R^3$
clear;
clc;
%[text] Declaración de variables simbólicas
syms theta1 omega1 theta2 omega2 l_c_2 real
disp('Variables definidas.') %[output:8ce8099b]
%[text] Definición de vectores
w1 = [0; 0; omega1] %[output:1246dc62]
w2 = omega2*[-sin(theta1); cos(theta1); 0] %[output:5fd5104a]
v_c2 = l_c_2*[cos(theta1)*sin(theta2); sin(theta1)*sin(theta2); cos(theta2)] %[output:58c2f488]
%[text] Calculos:
v_1 = crossproduct(w1, v_c2)' %[output:5e2522a2]
v_2 = crossproduct(w2, v_c2)' %[output:6a7cb1c7]
v   = v_1 + v_2;
v_mod = v' * v;
v_mod = simpTrig(v_mod);


disp('Velocidad del centro de masa'); %[output:479f771a]
pretty(v); %[output:5217dff9]

disp('Modulo cuadrado de la velocidad'); %[output:328911bf]
pretty(v_mod); %[output:2753a5c5]


disp(latex(v)); %[output:2109c709]
disp(latex(v_mod)); %[output:39826c1a]
%[text] ## EJEMPLOS
syms t theta phi a b c real
u = [cos(theta), sin(theta), 0];
v = [-sin(theta), cos(theta), 1];

%% ============================================================
%  PRODUCTO PUNTO
% =============================================================

dot_uv = dotproduct(u, v);

%% ============================================================
%  PRODUCTO VECTORIAL
% =============================================================

cross_uv = crossproduct(u, v);

%% ============================================================
%  MOSTRAR RESULTADOS
% =============================================================

disp('Vector u:') %[output:5202c036]
pretty(u) %[output:50a21e0e]

disp('Vector v:') %[output:6909b2bc]
pretty(v) %[output:9b22d302]

disp('Producto punto u · v:') %[output:0ba472a2]
pretty(dot_uv) %[output:08b9fc5e]

disp('Producto vectorial u × v:') %[output:291efba6]
pretty(cross_uv) %[output:5959b3ca]

%[text] ## FUNCIONES AUXILIARES
function dot_uv = dotproduct(u, v)
% u · v
dot_raw = sum(u .* v);

% Simplificación
dot_uv = simpTrig(dot_raw);

end

function cross_uv = crossproduct(u, v)

cross_raw = [ ...
    u(2)*v(3) - u(3)*v(2), ...
    u(3)*v(1) - u(1)*v(3), ...
    u(1)*v(2) - u(2)*v(1)];

% Simplificación
cross_uv = simpTrig(cross_raw);
end

function expr = simpTrig(expr)

% Reescribe tan(), cot(), etc. en términos de seno y coseno.
expr = rewrite(expr, 'sincos');

% Expande productos y potencias cuando esto ayuda a exponer
% identidades trigonométricas.
expr = expand(expr);

% Simplificación simbólica más agresiva.
expr = simplify(expr, 'Steps', 100);

end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:8ce8099b]
%   data: {"dataType":"text","outputData":{"text":"Variables definidas.\n","truncated":false}}
%---
%[output:1246dc62]
%   data: {"dataType":"symbolic","outputData":{"name":"w1","value":"\\left(\\begin{array}{c}\n0\\\\\n0\\\\\n\\omega_1 \n\\end{array}\\right)"}}
%---
%[output:5fd5104a]
%   data: {"dataType":"symbolic","outputData":{"name":"w2","value":"\\left(\\begin{array}{c}\n-\\omega_2 \\,\\sin \\left(\\theta_1 \\right)\\\\\n\\omega_2 \\,\\cos \\left(\\theta_1 \\right)\\\\\n0\n\\end{array}\\right)"}}
%---
%[output:58c2f488]
%   data: {"dataType":"symbolic","outputData":{"name":"v_c2","value":"\\left(\\begin{array}{c}\nl_{c,2} \\,\\cos \\left(\\theta_1 \\right)\\,\\sin \\left(\\theta_2 \\right)\\\\\nl_{c,2} \\,\\sin \\left(\\theta_1 \\right)\\,\\sin \\left(\\theta_2 \\right)\\\\\nl_{c,2} \\,\\cos \\left(\\theta_2 \\right)\n\\end{array}\\right)"}}
%---
%[output:5e2522a2]
%   data: {"dataType":"symbolic","outputData":{"name":"v_1","value":"\\left(\\begin{array}{c}\n-l_{c,2} \\,\\omega_1 \\,\\sin \\left(\\theta_1 \\right)\\,\\sin \\left(\\theta_2 \\right)\\\\\nl_{c,2} \\,\\omega_1 \\,\\cos \\left(\\theta_1 \\right)\\,\\sin \\left(\\theta_2 \\right)\\\\\n0\n\\end{array}\\right)"}}
%---
%[output:6a7cb1c7]
%   data: {"dataType":"symbolic","outputData":{"name":"v_2","value":"\\left(\\begin{array}{c}\nl_{c,2} \\,\\omega_2 \\,\\cos \\left(\\theta_1 \\right)\\,\\cos \\left(\\theta_2 \\right)\\\\\nl_{c,2} \\,\\omega_2 \\,\\cos \\left(\\theta_2 \\right)\\,\\sin \\left(\\theta_1 \\right)\\\\\n-l_{c,2} \\,\\omega_2 \\,\\sin \\left(\\theta_2 \\right)\n\\end{array}\\right)"}}
%---
%[output:479f771a]
%   data: {"dataType":"text","outputData":{"text":"Velocidad del centro de masa\n","truncated":false}}
%---
%[output:5217dff9]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\nl_{c,2} \\,\\omega_2 \\,\\cos \\left(\\theta_1 \\right)\\,\\cos \\left(\\theta_2 \\right)-l_{c,2} \\,\\omega_1 \\,\\sin \\left(\\theta_1 \\right)\\,\\sin \\left(\\theta_2 \\right)\\\\\nl_{c,2} \\,\\omega_1 \\,\\cos \\left(\\theta_1 \\right)\\,\\sin \\left(\\theta_2 \\right)+l_{c,2} \\,\\omega_2 \\,\\cos \\left(\\theta_2 \\right)\\,\\sin \\left(\\theta_1 \\right)\\\\\n-l_{c,2} \\,\\omega_2 \\,\\sin \\left(\\theta_2 \\right)\n\\end{array}\\right)"}}
%---
%[output:328911bf]
%   data: {"dataType":"text","outputData":{"text":"Modulo cuadrado de la velocidad\n","truncated":false}}
%---
%[output:2753a5c5]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"{l_{c,2} }^2 \\,{\\left({\\omega_1 }^2 \\,{\\sin \\left(\\theta_2 \\right)}^2 +{\\omega_2 }^2 \\right)}"}}
%---
%[output:2109c709]
%   data: {"dataType":"text","outputData":{"text":"\\left(\\begin{array}{c} l_{c,2}\\,\\omega _{2}\\,\\cos\\left(\\theta _{1}\\right)\\,\\cos\\left(\\theta _{2}\\right)-l_{c,2}\\,\\omega _{1}\\,\\sin\\left(\\theta _{1}\\right)\\,\\sin\\left(\\theta _{2}\\right)\\\\ l_{c,2}\\,\\omega _{1}\\,\\cos\\left(\\theta _{1}\\right)\\,\\sin\\left(\\theta _{2}\\right)+l_{c,2}\\,\\omega _{2}\\,\\cos\\left(\\theta _{2}\\right)\\,\\sin\\left(\\theta _{1}\\right)\\\\ -l_{c,2}\\,\\omega _{2}\\,\\sin\\left(\\theta _{2}\\right) \\end{array}\\right)\n","truncated":false}}
%---
%[output:39826c1a]
%   data: {"dataType":"text","outputData":{"text":"{l_{c,2}}^2\\,\\left({\\omega _{1}}^2\\,{\\sin\\left(\\theta _{2}\\right)}^2+{\\omega _{2}}^2\\right)\n","truncated":false}}
%---
%[output:5202c036]
%   data: {"dataType":"text","outputData":{"text":"Vector u:\n","truncated":false}}
%---
%[output:50a21e0e]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{ccc}\n\\cos \\left(\\theta \\right) & \\sin \\left(\\theta \\right) & 0\n\\end{array}\\right)"}}
%---
%[output:6909b2bc]
%   data: {"dataType":"text","outputData":{"text":"Vector v:\n","truncated":false}}
%---
%[output:9b22d302]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{ccc}\n-\\sin \\left(\\theta \\right) & \\cos \\left(\\theta \\right) & 1\n\\end{array}\\right)"}}
%---
%[output:0ba472a2]
%   data: {"dataType":"text","outputData":{"text":"Producto punto u · v:\n","truncated":false}}
%---
%[output:08b9fc5e]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"0"}}
%---
%[output:291efba6]
%   data: {"dataType":"text","outputData":{"text":"Producto vectorial u × v:\n","truncated":false}}
%---
%[output:5959b3ca]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{ccc}\n\\sin \\left(\\theta \\right) & -\\cos \\left(\\theta \\right) & 1\n\\end{array}\\right)"}}
%---
