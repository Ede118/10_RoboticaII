function q_solucion = FncCinematicaInversa(Robot, P)
% FNCCINEMATICAINVERSA Calcula la cinemática inversa analítica para un robot de 3-GDL.
%
% Inputs:
%   Robot: Objeto SerialLink de 3-GDL (contiene las longitudes a2 y a3)
%   P: Vector de posición [Px, Py, Pz]
%
% Output:
%   q_solucion: Configuración articular [q1, q2, q3] (Codo Arriba)

Px = P(1);
Py = P(2);
Pz = P(3);

a2 = Robot.links(2).a;
a3 = Robot.links(3).a;

q1 = atan2(Py, Px);

r = sqrt(Px^2 + Py^2);

cos_q3 = (r^2 + Pz^2 - a2^2 - a3^2) / (2 * a2 * a3);

% Verificacion del punto del robot
if abs(cos_q3) > 1
    error('El punto objetivo esta fuera del espacio de trabajo del robot.');
end

% Solucion Codo arriba
sin_q3 = -sqrt(1 - cos_q3^2); 
q3 = atan2(sin_q3, cos_q3);

% Cálculo de q2
q2 = atan2(Pz, r) - atan2(a3 * sin(q3), a2 + a3 * cos(q3));

q_solucion = [q1, q2, q3];

end
