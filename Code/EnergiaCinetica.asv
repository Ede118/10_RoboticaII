%[text] # PRODUCTO PUNTO Y PRODUCTO VECTORIAL SIMBÓLICO EN $R^3$
clear;
clc;
%[text] ## Variables Compartidas
%[text] Matriz de velocidades angulares
%[text] $\\underline{w\_i} = w(:, i)$
w  = sym(zeros(3,3));

%[text] Matriz de velocidades lineales
%[text] $\\underline{v\_i} = v(:, i)$
vc = sym(zeros(3,3));
%[text] Matriz de Energia cinética
%[text] K(1, i) = Energia Cinética Lineal del eslabón i
%[text] K(2, i) = Energia Cinética Rotacional del eslabón i
%[text] K(3, i) = Energia Cinética Total del eslabón i
K  = sym(zeros(3,3));
%%
%[text] ## Energía cinética $K\_1$
syms q1 q_dot_1 I_1_xx I_1_yy I_1_zz m1 real

I_c1 = [I_1_xx, 0, 0; 0, I_1_yy, 0; 0, 0, I_1_zz];

disp('Variables definidas.') %[output:5888ea94]
disp(q1) %[output:33714adc]
disp(q_dot_1) %[output:221c45b4]
disp(I_c1) %[output:3c345edd]


w(:, 1)  = q_dot_1*[0; 0; 1];
vc(:, 1) = [0; 0; 0];

K(1,1) = sym(1)/2 * m1 * vc(:,1)' * vc(:,1);
K(2,1) = 0.5 * w(:,1)' * I_c1 * w(:,1);      
K(3,1) = K(1, 1) + K(2, 1);                  

disp('omega1'); disp(w(:,1)); %[output:572dba19] %[output:029ad126]
disp('v_c_1'); disp(vc(:,1)); %[output:7437d110] %[output:804edd63]
disp('K1 lineal'); disp(K(1,1)); %[output:79eee012] %[output:358bb781]
disp('K1 rotacional'); disp(K(2,1)); %[output:239a3d6f] %[output:983bdcb5]
disp('K1 total'); disp(K(3,1)); %[output:856b2e38] %[output:9633d7a1]
%%
%[text] ## Energía cinética $K\_2$
syms q2 q_dot_2 l_c_2 I_2_xx I_2_yy I_2_zz m2 real

I_c2 = [I_2_xx, 0, 0; 0, I_2_yy, 0; 0, 0, I_2_zz];

disp('Variables definidas.') %[output:5391f4dc]

disp(q2) %[output:66b1fe90]
disp(q_dot_2) %[output:0b284933]
disp(l_c_2) %[output:9d9d9639]
disp(I_c2) %[output:554cb999]

% Velocidad angular relativa 2/1
w21 = q_dot_2*[sin(q1); -cos(q1); 0] %[output:543a7623]

w(:,2) = w(:, 1) + w21;

w2_mod = w(:,2)' * w(:,2);
w2_mod = simplificar(w2_mod);

rO2C2 = l_c_2 * [-cos(q1)*sin(q2); -sin(q1)*sin(q2); cos(q2)];

vc(:,2) = simplify(cross(w(:,2), rO2C2)');


K(1,2) = sym(1)/2 * m2 * vc(:,2)' * vc(:,2);
K(2,2) = 0.5 * w(:,2)' * I_c2 * w(:,2);      
K(3,2) = K(1, 2) + K(2, 2);                  

q_dot = [q_dot_1; q_dot_2];
[best, ranking] = simplifyCompact(K(3,2), 200);
[best, ~] = quadraticCollect(best, q_dot);
K(3,2) = best;

disp('omega2'); disp(w(:,2)); %[output:7094eda1] %[output:6cfb581e]
disp('v_c_2'); disp(vc(:,2)); %[output:6a3fc521] %[output:9fec0789]
disp('K2 lineal'); disp(K(1,2)); %[output:2cb7cb6d] %[output:1e0ab48f]
disp('K2 rotacional'); disp(K(2,2)); %[output:4174e359] %[output:5fbde242]
disp('K2 total'); disp(K(3,2)); %[output:3f850793] %[output:5a61d626]
%%
%[text] ## Energía cinética $K\_3$
syms q3 q_dot_3 l_c_3 L2 I_3_xx I_3_yy I_3_zz m3 real

I_c3 = [I_3_xx, 0, 0; 0, I_3_yy, 0; 0, 0, I_3_zz];

disp('Variables definidas.') %[output:57b26a47]

w32 = q_dot_3*[sin(q1); -cos(q1); 0];
w(:,3) = w(:,2) + w32;
w3_mod = w(:,3)' * w(:,3);
w3_mod = simplificar(w3_mod);

rO2O3 = L2 * [-cos(q1)*sin(q2); -sin(q1)*sin(q2); cos(q2)];
rO3C3 = l_c_3 * [-cos(q1)*sin(q2+q3); -sin(q1)*sin(q2+q3); cos(q2+q3)];

v_O3 = simplify(cross(w(:,2), rO2O3)');
v32 = simplify(cross(w(:,3), rO3C3)');

vc(:,3) = v_O3 + v32;
v_c_3_mod = vc(:,3)' * vc(:,3);

K(1,3) = sym(1)/2 * m2 * vc(:,3)' * vc(:,3);
K(2,3) = 0.5 * w(:,3)' * I_c2 * w(:,3);      
K(3,3) = K(1, 3) + K(2, 3);  

q_dot = [q_dot_1; q_dot_2; q_dot_3];
[best, ranking] = simplifyCompact(K(3,3), 200);
[best, ~] = quadraticCollect(best, q_dot);

K(3,3) = best;

disp('omega3'); disp(w(:,3)); %[output:4279093a] %[output:807f5752]
disp('v_c_3'); disp(vc(:,3)); %[output:330d77ce] %[output:00d22f83]
disp('K3 lineal'); disp(K(1,3)); %[output:252ee814] %[output:35534f56]
disp('K3 rotacional'); disp(K(2,3)); %[output:954fb10c] %[output:9500456d]
disp('K3 total'); disp(K(3,3)); %[output:6a9de17b] %[output:3d5f4830]

%%
%[text] ## Energía Cinética Total $K$
K_total = K(3,1) + K(3,2) + K(3,3);

q_dot = [q_dot_1; q_dot_2; q_dot_3];
[best, ranking] = simplifyCompact(K_total, 200);
[best, A] = quadraticCollect(best, q_dot);
K_total = best;

disp('Energia Cinética K'); disp(K_total); disp(latex(K_total)); %[output:7ef652c3] %[output:142c571f] %[output:040d6e47]

coef = {
    'q_1^2',       A(1,1);
    'q_2^2',       A(2,2);
    'q_3^2',       A(3,3);
    'q_2 q_3',   2*A(2,3)
    };

for k = 1:size(coef,1) %[output:group:1ffc9c19]
    fprintf('\nCoeficiente de %s:\n', coef{k,1}); %[output:7f39a6fc] %[output:2d6a7ce0] %[output:712422de] %[output:951e5899]
    disp(coef{k,2}); %[output:077dc5c1] %[output:2e952745] %[output:3513f52d] %[output:9378bc7c]
    fprintf('LaTeX: %s\n', latex(coef{k,2})); %[output:704206c4] %[output:36b56336] %[output:6b8a4ea2] %[output:1fb8c5dd]

end %[output:group:1ffc9c19]
%[text] ## FUNCIONES AUXILIARES
function dot_uv = dotproduct(u, v)
% u · v
dot_raw = sum(u .* v);

% Simplificación
dot_uv = simplificar(dot_raw);

end

function cross_uv = crossproduct(u, v)

cross_raw = [ ...
    u(2)*v(3) - u(3)*v(2), ...
    u(3)*v(1) - u(1)*v(3), ...
    u(1)*v(2) - u(2)*v(1)];

% Simplificación
cross_uv = simplificar(cross_raw);
end

function expr = simplificar(expr)

% Reescribe tan(), cot(), etc. en términos de seno y coseno.
expr = rewrite(expr, 'sincos');

% Expande productos y potencias cuando esto ayuda a exponer
% identidades trigonométricas.
expr = expand(expr);

% Simplificación simbólica más agresiva.
expr = simplify(expr, 'Steps', 100);

end

function [expr_ordered, A] = quadraticCollect(expr, q_dot)

% Simplificación inicial
expr = simplify(expr, 'Steps', 200);

% Matriz de la forma cuadrática
% expr = q_dot.' * A * q_dot
A = hessian(expr, q_dot) / 2;

n = length(q_dot);

% Simplificar cada coeficiente de la matriz
for i = 1:n
    for j = 1:n

        aux = simplify(A(i,j), 'Steps', 100);

        % Intenta reagrupar identidades trigonométricas
        aux = combine(aux, 'sincos');

        % Simplificación final
        aux = simplify(aux, 'Steps', 200);

        A(i,j) = aux;
    end
end

% Reconstruir la expresión ordenada según
% q_dot_i^2 y q_dot_i*q_dot_j
expr_ordered = sym(0);

for i = 1:n

    % Términos diagonales: Aii*q_dot_i^2
    expr_ordered = expr_ordered ...
        + A(i,i)*q_dot(i)^2;

    % Términos cruzados: 2*Aij*q_dot_i*q_dot_j
    for j = i+1:n
        expr_ordered = expr_ordered ...
            + 2*A(i,j)*q_dot(i)*q_dot(j);
    end
end

expr_ordered = simplify(expr_ordered, 'Steps', 200);

end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
%[output:5888ea94]
%   data: {"dataType":"text","outputData":{"text":"Variables definidas.\n","truncated":false}}
%---
%[output:33714adc]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"q_1"}}
%---
%[output:221c45b4]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"{\\dot{q} }_1"}}
%---
%[output:3c345edd]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{ccc}\nI_{1,\\textrm{xx}}  & 0 & 0\\\\\n0 & I_{1,\\textrm{yy}}  & 0\\\\\n0 & 0 & I_{1,\\textrm{zz}} \n\\end{array}\\right)"}}
%---
%[output:572dba19]
%   data: {"dataType":"text","outputData":{"text":"omega1\n","truncated":false}}
%---
%[output:029ad126]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\n0\\\\\n0\\\\\n{\\dot{q} }_1 \n\\end{array}\\right)"}}
%---
%[output:7437d110]
%   data: {"dataType":"text","outputData":{"text":"v_c_1\n","truncated":false}}
%---
%[output:804edd63]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\n0\\\\\n0\\\\\n0\n\\end{array}\\right)"}}
%---
%[output:79eee012]
%   data: {"dataType":"text","outputData":{"text":"K1 lineal\n","truncated":false}}
%---
%[output:358bb781]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"0"}}
%---
%[output:239a3d6f]
%   data: {"dataType":"text","outputData":{"text":"K1 rotacional\n","truncated":false}}
%---
%[output:983bdcb5]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{I_{1,\\textrm{zz}} \\,{{\\dot{q} }_1 }^2 }{2}"}}
%---
%[output:856b2e38]
%   data: {"dataType":"text","outputData":{"text":"K1 total\n","truncated":false}}
%---
%[output:9633d7a1]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{I_{1,\\textrm{zz}} \\,{{\\dot{q} }_1 }^2 }{2}"}}
%---
%[output:5391f4dc]
%   data: {"dataType":"text","outputData":{"text":"Variables definidas.\n","truncated":false}}
%---
%[output:66b1fe90]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"q_2"}}
%---
%[output:0b284933]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"{\\dot{q} }_2"}}
%---
%[output:9d9d9639]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"l_{c,2}"}}
%---
%[output:554cb999]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{ccc}\nI_{2,\\textrm{xx}}  & 0 & 0\\\\\n0 & I_{2,\\textrm{yy}}  & 0\\\\\n0 & 0 & I_{2,\\textrm{zz}} \n\\end{array}\\right)"}}
%---
%[output:543a7623]
%   data: {"dataType":"symbolic","outputData":{"name":"w21","value":"\\left(\\begin{array}{c}\n{\\dot{q} }_2 \\,\\sin \\left(q_1 \\right)\\\\\n-{\\dot{q} }_2 \\,\\cos \\left(q_1 \\right)\\\\\n0\n\\end{array}\\right)"}}
%---
%[output:7094eda1]
%   data: {"dataType":"text","outputData":{"text":"omega2\n","truncated":false}}
%---
%[output:6cfb581e]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\n{\\dot{q} }_2 \\,\\sin \\left(q_1 \\right)\\\\\n-{\\dot{q} }_2 \\,\\cos \\left(q_1 \\right)\\\\\n{\\dot{q} }_1 \n\\end{array}\\right)"}}
%---
%[output:6a3fc521]
%   data: {"dataType":"text","outputData":{"text":"v_c_2\n","truncated":false}}
%---
%[output:9fec0789]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\nl_{c,2} \\,{\\dot{q} }_1 \\,\\sin \\left(q_1 \\right)\\,\\sin \\left(q_2 \\right)-l_{c,2} \\,{\\dot{q} }_2 \\,\\cos \\left(q_1 \\right)\\,\\cos \\left(q_2 \\right)\\\\\n-l_{c,2} \\,{\\dot{q} }_1 \\,\\cos \\left(q_1 \\right)\\,\\sin \\left(q_2 \\right)-l_{c,2} \\,{\\dot{q} }_2 \\,\\cos \\left(q_2 \\right)\\,\\sin \\left(q_1 \\right)\\\\\n-l_{c,2} \\,{\\dot{q} }_2 \\,\\sin \\left(q_2 \\right)\n\\end{array}\\right)"}}
%---
%[output:2cb7cb6d]
%   data: {"dataType":"text","outputData":{"text":"K2 lineal\n","truncated":false}}
%---
%[output:1e0ab48f]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{m_2 \\,{{\\left(l_{c,2} \\,{\\dot{q} }_1 \\,\\cos \\left(q_1 \\right)\\,\\sin \\left(q_2 \\right)+l_{c,2} \\,{\\dot{q} }_2 \\,\\cos \\left(q_2 \\right)\\,\\sin \\left(q_1 \\right)\\right)}}^2 }{2}+\\frac{m_2 \\,{{\\left(l_{c,2} \\,{\\dot{q} }_1 \\,\\sin \\left(q_1 \\right)\\,\\sin \\left(q_2 \\right)-l_{c,2} \\,{\\dot{q} }_2 \\,\\cos \\left(q_1 \\right)\\,\\cos \\left(q_2 \\right)\\right)}}^2 }{2}+\\frac{{l_{c,2} }^2 \\,m_2 \\,{{\\dot{q} }_2 }^2 \\,{\\sin \\left(q_2 \\right)}^2 }{2}"}}
%---
%[output:4174e359]
%   data: {"dataType":"text","outputData":{"text":"K2 rotacional\n","truncated":false}}
%---
%[output:5fbde242]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{I_{2,\\textrm{zz}} \\,{{\\dot{q} }_1 }^2 }{2}+\\frac{I_{2,\\textrm{yy}} \\,{{\\dot{q} }_2 }^2 \\,{\\cos \\left(q_1 \\right)}^2 }{2}+\\frac{I_{2,\\textrm{xx}} \\,{{\\dot{q} }_2 }^2 \\,{\\sin \\left(q_1 \\right)}^2 }{2}"}}
%---
%[output:3f850793]
%   data: {"dataType":"text","outputData":{"text":"K2 total\n","truncated":false}}
%---
%[output:5a61d626]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"{\\left(\\frac{m_2 \\,{l_{c,2} }^2 \\,{\\sin \\left(q_2 \\right)}^2 }{2}+\\frac{I_{2,\\textrm{zz}} }{2}\\right)}\\,{{\\dot{q} }_1 }^2 +{\\left(\\frac{I_{2,\\textrm{xx}} }{2}+\\frac{{l_{c,2} }^2 \\,m_2 }{2}-\\frac{I_{2,\\textrm{xx}} \\,{\\cos \\left(q_1 \\right)}^2 }{2}+\\frac{I_{2,\\textrm{yy}} \\,{\\cos \\left(q_1 \\right)}^2 }{2}\\right)}\\,{{\\dot{q} }_2 }^2"}}
%---
%[output:57b26a47]
%   data: {"dataType":"text","outputData":{"text":"Variables definidas.\n","truncated":false}}
%---
%[output:4279093a]
%   data: {"dataType":"text","outputData":{"text":"omega3\n","truncated":false}}
%---
%[output:807f5752]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\n{\\dot{q} }_2 \\,\\sin \\left(q_1 \\right)+{\\dot{q} }_3 \\,\\sin \\left(q_1 \\right)\\\\\n-{\\dot{q} }_2 \\,\\cos \\left(q_1 \\right)-{\\dot{q} }_3 \\,\\cos \\left(q_1 \\right)\\\\\n{\\dot{q} }_1 \n\\end{array}\\right)"}}
%---
%[output:330d77ce]
%   data: {"dataType":"text","outputData":{"text":"v_c_3\n","truncated":false}}
%---
%[output:00d22f83]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{c}\nl_{c,3} \\,{\\dot{q} }_1 \\,\\sin \\left(q_2 +q_3 \\right)\\,\\sin \\left(q_1 \\right)-L_2 \\,{\\dot{q} }_2 \\,\\cos \\left(q_1 \\right)\\,\\cos \\left(q_2 \\right)-l_{c,3} \\,\\cos \\left(q_2 +q_3 \\right)\\,{\\left({\\dot{q} }_2 \\,\\cos \\left(q_1 \\right)+{\\dot{q} }_3 \\,\\cos \\left(q_1 \\right)\\right)}+L_2 \\,{\\dot{q} }_1 \\,\\sin \\left(q_1 \\right)\\,\\sin \\left(q_2 \\right)\\\\\n-l_{c,3} \\,\\cos \\left(q_2 +q_3 \\right)\\,{\\left({\\dot{q} }_2 \\,\\sin \\left(q_1 \\right)+{\\dot{q} }_3 \\,\\sin \\left(q_1 \\right)\\right)}-l_{c,3} \\,{\\dot{q} }_1 \\,\\sin \\left(q_2 +q_3 \\right)\\,\\cos \\left(q_1 \\right)-L_2 \\,{\\dot{q} }_1 \\,\\cos \\left(q_1 \\right)\\,\\sin \\left(q_2 \\right)-L_2 \\,{\\dot{q} }_2 \\,\\cos \\left(q_2 \\right)\\,\\sin \\left(q_1 \\right)\\\\\n-L_2 \\,{\\dot{q} }_2 \\,\\sin \\left(q_2 \\right)-l_{c,3} \\,\\sin \\left(q_2 +q_3 \\right)\\,{\\left({\\dot{q} }_2 +{\\dot{q} }_3 \\right)}\n\\end{array}\\right)"}}
%---
%[output:252ee814]
%   data: {"dataType":"text","outputData":{"text":"K3 lineal\n","truncated":false}}
%---
%[output:35534f56]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{m_2 \\,{{\\left(l_{c,3} \\,\\cos \\left(q_2 +q_3 \\right)\\,{\\left({\\dot{q} }_2 \\,\\cos \\left(q_1 \\right)+{\\dot{q} }_3 \\,\\cos \\left(q_1 \\right)\\right)}+L_2 \\,{\\dot{q} }_2 \\,\\cos \\left(q_1 \\right)\\,\\cos \\left(q_2 \\right)-l_{c,3} \\,{\\dot{q} }_1 \\,\\sin \\left(q_2 +q_3 \\right)\\,\\sin \\left(q_1 \\right)-L_2 \\,{\\dot{q} }_1 \\,\\sin \\left(q_1 \\right)\\,\\sin \\left(q_2 \\right)\\right)}}^2 }{2}+\\frac{m_2 \\,{{\\left(l_{c,3} \\,\\cos \\left(q_2 +q_3 \\right)\\,{\\left({\\dot{q} }_2 \\,\\sin \\left(q_1 \\right)+{\\dot{q} }_3 \\,\\sin \\left(q_1 \\right)\\right)}+l_{c,3} \\,{\\dot{q} }_1 \\,\\sin \\left(q_2 +q_3 \\right)\\,\\cos \\left(q_1 \\right)+L_2 \\,{\\dot{q} }_1 \\,\\cos \\left(q_1 \\right)\\,\\sin \\left(q_2 \\right)+L_2 \\,{\\dot{q} }_2 \\,\\cos \\left(q_2 \\right)\\,\\sin \\left(q_1 \\right)\\right)}}^2 }{2}+\\frac{m_2 \\,{{\\left(L_2 \\,{\\dot{q} }_2 \\,\\sin \\left(q_2 \\right)+l_{c,3} \\,\\sin \\left(q_2 +q_3 \\right)\\,{\\left({\\dot{q} }_2 +{\\dot{q} }_3 \\right)}\\right)}}^2 }{2}"}}
%---
%[output:954fb10c]
%   data: {"dataType":"text","outputData":{"text":"K3 rotacional\n","truncated":false}}
%---
%[output:9500456d]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{I_{2,\\textrm{zz}} \\,{{\\dot{q} }_1 }^2 }{2}+I_{2,\\textrm{xx}} \\,{\\left({\\dot{q} }_2 \\,\\sin \\left(q_1 \\right)+{\\dot{q} }_3 \\,\\sin \\left(q_1 \\right)\\right)}\\,{\\left(\\frac{{\\dot{q} }_2 \\,\\sin \\left(q_1 \\right)}{2}+\\frac{{\\dot{q} }_3 \\,\\sin \\left(q_1 \\right)}{2}\\right)}+I_{2,\\textrm{yy}} \\,{\\left({\\dot{q} }_2 \\,\\cos \\left(q_1 \\right)+{\\dot{q} }_3 \\,\\cos \\left(q_1 \\right)\\right)}\\,{\\left(\\frac{{\\dot{q} }_2 \\,\\cos \\left(q_1 \\right)}{2}+\\frac{{\\dot{q} }_3 \\,\\cos \\left(q_1 \\right)}{2}\\right)}"}}
%---
%[output:6a9de17b]
%   data: {"dataType":"text","outputData":{"text":"K3 total\n","truncated":false}}
%---
%[output:3d5f4830]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\begin{array}{l}\n{\\left(\\frac{m_2 \\,{L_2 }^2 \\,{\\sin \\left(q_2 \\right)}^2 }{2}+m_2 \\,L_2 \\,l_{c,3} \\,{\\sin \\left(q_2 +\\frac{q_3 }{2}\\right)}^2 -m_2 \\,L_2 \\,l_{c,3} \\,\\sigma_2 +\\frac{m_2 \\,{l_{c,3} }^2 \\,{\\sin \\left(q_2 +q_3 \\right)}^2 }{2}+\\frac{I_{2,\\textrm{zz}} }{2}\\right)}\\,{{\\dot{q} }_1 }^2 +{\\left(\\frac{m_2 \\,{L_2 }^2 }{2}-m_2 \\,{\\left(2\\,\\sigma_2 -1\\right)}\\,L_2 \\,l_{c,3} +\\sigma_1 +\\frac{I_{2,\\textrm{xx}} }{2}+\\frac{I_{2,\\textrm{xx}} \\,\\sigma_3 }{2}-\\frac{I_{2,\\textrm{yy}} \\,\\sigma_3 }{2}\\right)}\\,{{\\dot{q} }_2 }^2 +{\\left({l_{c,3} }^2 \\,m_2 -I_{2,\\textrm{yy}} \\,\\sigma_3 +I_{2,\\textrm{xx}} \\,{\\sin \\left(q_1 \\right)}^2 -L_2 \\,l_{c,3} \\,m_2 \\,{\\left(2\\,\\sigma_2 -1\\right)}\\right)}\\,{\\dot{q} }_2 \\,{\\dot{q} }_3 +{\\left(\\sigma_1 +\\frac{I_{2,\\textrm{xx}} }{2}+\\frac{I_{2,\\textrm{xx}} \\,\\sigma_3 }{2}-\\frac{I_{2,\\textrm{yy}} \\,\\sigma_3 }{2}\\right)}\\,{{\\dot{q} }_3 }^2 \\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\frac{m_2 \\,{l_{c,3} }^2 }{2}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 ={\\sin \\left(\\frac{q_3 }{2}\\right)}^2 \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_3 ={\\sin \\left(q_1 \\right)}^2 -1\n\\end{array}"}}
%---
%[output:7ef652c3]
%   data: {"dataType":"text","outputData":{"text":"Energia Cinética K\n","truncated":false}}
%---
%[output:142c571f]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\begin{array}{l}\n{\\left(\\frac{I_{1,\\textrm{zz}} }{2}+I_{2,\\textrm{zz}} +\\sigma_6 +\\sigma_5 +\\frac{{l_{c,3} }^2 \\,m_2 }{4}-\\frac{{L_2 }^2 \\,m_2 \\,{\\cos \\left(q_2 \\right)}^2 }{2}-\\frac{{l_{c,2} }^2 \\,m_2 \\,{\\cos \\left(q_2 \\right)}^2 }{2}-\\frac{{l_{c,3} }^2 \\,m_2 \\,\\cos \\left(2\\,q_2 +2\\,q_3 \\right)}{4}+\\frac{\\sigma_1 }{2}-\\frac{L_2 \\,l_{c,3} \\,m_2 \\,\\cos \\left(2\\,q_2 +q_3 \\right)}{2}\\right)}\\,{{\\dot{q} }_1 }^2 +{\\left(I_{2,\\textrm{xx}} +\\sigma_6 +\\sigma_5 +\\sigma_4 -\\sigma_3 +\\sigma_2 +\\sigma_1 \\right)}\\,{{\\dot{q} }_2 }^2 +{\\left(I_{2,\\textrm{xx}} +{l_{c,3} }^2 \\,m_2 -\\sigma_3 +\\sigma_2 +\\sigma_1 \\right)}\\,{\\dot{q} }_2 \\,{\\dot{q} }_3 +{\\left(\\frac{I_{2,\\textrm{xx}} }{2}+\\sigma_4 -\\frac{\\sigma_3 }{2}+\\frac{\\sigma_2 }{2}\\right)}\\,{{\\dot{q} }_3 }^2 \\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =L_2 \\,l_{c,3} \\,m_2 \\,\\cos \\left(q_3 \\right)\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 =I_{2,\\textrm{yy}} \\,{\\cos \\left(q_1 \\right)}^2 \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_3 =I_{2,\\textrm{xx}} \\,{\\cos \\left(q_1 \\right)}^2 \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_4 =\\frac{{l_{c,3} }^2 \\,m_2 }{2}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_5 =\\frac{{l_{c,2} }^2 \\,m_2 }{2}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_6 =\\frac{{L_2 }^2 \\,m_2 }{2}\n\\end{array}"}}
%---
%[output:040d6e47]
%   data: {"dataType":"text","outputData":{"text":"\\left(\\frac{I_{1,\\mathrm{zz}}}{2}+I_{2,\\mathrm{zz}}+\\frac{{L_{2}}^2\\,m_{2}}{2}+\\frac{{l_{c,2}}^2\\,m_{2}}{2}+\\frac{{l_{c,3}}^2\\,m_{2}}{4}-\\frac{{L_{2}}^2\\,m_{2}\\,{\\cos\\left(q_{2}\\right)}^2}{2}-\\frac{{l_{c,2}}^2\\,m_{2}\\,{\\cos\\left(q_{2}\\right)}^2}{2}-\\frac{{l_{c,3}}^2\\,m_{2}\\,\\cos\\left(2\\,q_{2}+2\\,q_{3}\\right)}{4}+\\frac{L_{2}\\,l_{c,3}\\,m_{2}\\,\\cos\\left(q_{3}\\right)}{2}-\\frac{L_{2}\\,l_{c,3}\\,m_{2}\\,\\cos\\left(2\\,q_{2}+q_{3}\\right)}{2}\\right)\\,{\\dot{q}_{1}}^2+\\left(I_{2,\\mathrm{xx}}+\\frac{{L_{2}}^2\\,m_{2}}{2}+\\frac{{l_{c,2}}^2\\,m_{2}}{2}+\\frac{{l_{c,3}}^2\\,m_{2}}{2}-I_{2,\\mathrm{xx}}\\,{\\cos\\left(q_{1}\\right)}^2+I_{2,\\mathrm{yy}}\\,{\\cos\\left(q_{1}\\right)}^2+L_{2}\\,l_{c,3}\\,m_{2}\\,\\cos\\left(q_{3}\\right)\\right)\\,{\\dot{q}_{2}}^2+\\left(I_{2,\\mathrm{xx}}+{l_{c,3}}^2\\,m_{2}-I_{2,\\mathrm{xx}}\\,{\\cos\\left(q_{1}\\right)}^2+I_{2,\\mathrm{yy}}\\,{\\cos\\left(q_{1}\\right)}^2+L_{2}\\,l_{c,3}\\,m_{2}\\,\\cos\\left(q_{3}\\right)\\right)\\,\\dot{q}_{2}\\,\\dot{q}_{3}+\\left(\\frac{I_{2,\\mathrm{xx}}}{2}+\\frac{{l_{c,3}}^2\\,m_{2}}{2}-\\frac{I_{2,\\mathrm{xx}}\\,{\\cos\\left(q_{1}\\right)}^2}{2}+\\frac{I_{2,\\mathrm{yy}}\\,{\\cos\\left(q_{1}\\right)}^2}{2}\\right)\\,{\\dot{q}_{3}}^2\n","truncated":false}}
%---
%[output:7f39a6fc]
%   data: {"dataType":"text","outputData":{"text":"\nCoeficiente de q_1^2:\n","truncated":false}}
%---
%[output:077dc5c1]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{I_{1,\\textrm{zz}} }{2}+I_{2,\\textrm{zz}} +\\frac{{L_2 }^2 \\,m_2 }{2}+\\frac{{l_{c,2} }^2 \\,m_2 }{2}+\\frac{{l_{c,3} }^2 \\,m_2 }{4}-\\frac{{L_2 }^2 \\,m_2 \\,{\\cos \\left(q_2 \\right)}^2 }{2}-\\frac{{l_{c,2} }^2 \\,m_2 \\,{\\cos \\left(q_2 \\right)}^2 }{2}-\\frac{{l_{c,3} }^2 \\,m_2 \\,\\cos \\left(2\\,q_2 +2\\,q_3 \\right)}{4}+\\frac{L_2 \\,l_{c,3} \\,m_2 \\,\\cos \\left(q_3 \\right)}{2}-\\frac{L_2 \\,l_{c,3} \\,m_2 \\,\\cos \\left(2\\,q_2 +q_3 \\right)}{2}"}}
%---
%[output:704206c4]
%   data: {"dataType":"text","outputData":{"text":"LaTeX: \\frac{I_{1,\\mathrm{zz}}}{2}+I_{2,\\mathrm{zz}}+\\frac{{L_{2}}^2\\,m_{2}}{2}+\\frac{{l_{c,2}}^2\\,m_{2}}{2}+\\frac{{l_{c,3}}^2\\,m_{2}}{4}-\\frac{{L_{2}}^2\\,m_{2}\\,{\\cos\\left(q_{2}\\right)}^2}{2}-\\frac{{l_{c,2}}^2\\,m_{2}\\,{\\cos\\left(q_{2}\\right)}^2}{2}-\\frac{{l_{c,3}}^2\\,m_{2}\\,\\cos\\left(2\\,q_{2}+2\\,q_{3}\\right)}{4}+\\frac{L_{2}\\,l_{c,3}\\,m_{2}\\,\\cos\\left(q_{3}\\right)}{2}-\\frac{L_{2}\\,l_{c,3}\\,m_{2}\\,\\cos\\left(2\\,q_{2}+q_{3}\\right)}{2}\n","truncated":false}}
%---
%[output:2d6a7ce0]
%   data: {"dataType":"text","outputData":{"text":"\nCoeficiente de q_2^2:\n","truncated":false}}
%---
%[output:2e952745]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"I_{2,\\textrm{xx}} +\\frac{{L_2 }^2 \\,m_2 }{2}+\\frac{{l_{c,2} }^2 \\,m_2 }{2}+\\frac{{l_{c,3} }^2 \\,m_2 }{2}-I_{2,\\textrm{xx}} \\,{\\cos \\left(q_1 \\right)}^2 +I_{2,\\textrm{yy}} \\,{\\cos \\left(q_1 \\right)}^2 +L_2 \\,l_{c,3} \\,m_2 \\,\\cos \\left(q_3 \\right)"}}
%---
%[output:36b56336]
%   data: {"dataType":"text","outputData":{"text":"LaTeX: I_{2,\\mathrm{xx}}+\\frac{{L_{2}}^2\\,m_{2}}{2}+\\frac{{l_{c,2}}^2\\,m_{2}}{2}+\\frac{{l_{c,3}}^2\\,m_{2}}{2}-I_{2,\\mathrm{xx}}\\,{\\cos\\left(q_{1}\\right)}^2+I_{2,\\mathrm{yy}}\\,{\\cos\\left(q_{1}\\right)}^2+L_{2}\\,l_{c,3}\\,m_{2}\\,\\cos\\left(q_{3}\\right)\n","truncated":false}}
%---
%[output:712422de]
%   data: {"dataType":"text","outputData":{"text":"\nCoeficiente de q_3^2:\n","truncated":false}}
%---
%[output:3513f52d]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{I_{2,\\textrm{xx}} }{2}+\\frac{{l_{c,3} }^2 \\,m_2 }{2}-\\frac{I_{2,\\textrm{xx}} \\,{\\cos \\left(q_1 \\right)}^2 }{2}+\\frac{I_{2,\\textrm{yy}} \\,{\\cos \\left(q_1 \\right)}^2 }{2}"}}
%---
%[output:6b8a4ea2]
%   data: {"dataType":"text","outputData":{"text":"LaTeX: \\frac{I_{2,\\mathrm{xx}}}{2}+\\frac{{l_{c,3}}^2\\,m_{2}}{2}-\\frac{I_{2,\\mathrm{xx}}\\,{\\cos\\left(q_{1}\\right)}^2}{2}+\\frac{I_{2,\\mathrm{yy}}\\,{\\cos\\left(q_{1}\\right)}^2}{2}\n","truncated":false}}
%---
%[output:951e5899]
%   data: {"dataType":"text","outputData":{"text":"\nCoeficiente de q_2 q_3:\n","truncated":false}}
%---
%[output:9378bc7c]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"I_{2,\\textrm{xx}} +{l_{c,3} }^2 \\,m_2 -I_{2,\\textrm{xx}} \\,{\\cos \\left(q_1 \\right)}^2 +I_{2,\\textrm{yy}} \\,{\\cos \\left(q_1 \\right)}^2 +L_2 \\,l_{c,3} \\,m_2 \\,\\cos \\left(q_3 \\right)"}}
%---
%[output:1fb8c5dd]
%   data: {"dataType":"text","outputData":{"text":"LaTeX: I_{2,\\mathrm{xx}}+{l_{c,3}}^2\\,m_{2}-I_{2,\\mathrm{xx}}\\,{\\cos\\left(q_{1}\\right)}^2+I_{2,\\mathrm{yy}}\\,{\\cos\\left(q_{1}\\right)}^2+L_{2}\\,l_{c,3}\\,m_{2}\\,\\cos\\left(q_{3}\\right)\n","truncated":false}}
%---
