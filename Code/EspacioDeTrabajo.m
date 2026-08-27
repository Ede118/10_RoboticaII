%% Ley temporal
clc; clear; close all;
% Cargar el modelo del robot
RobotSerie;

q1_range = robot.qlim(1, :);
q2_range = robot.qlim(2, :);
q3_range = robot.qlim(3, :);

paso = 0.08;
q1_vec = q1_range(1):paso:q1_range(2);
q2_vec = q2_range(1):paso:q2_range(2);

P_max = []; 
P_min = []; 

q3_max = 0;

for q2 = q2_vec
    for q1 = q1_vec
        T = robot.fkine([q1, q2, q3_max]);

        if isobject(T)
            pos = T.t;        
        else
            pos = T(1:3, 4);  
        end
        
        P_max = [P_max; pos(:)'];
    end
end

q3_min = q3_range(1); 

for q2 = q2_vec
    for q1 = q1_vec
        T = robot.fkine([q1, q2, q3_min]);
        
        if isobject(T)
            pos = T.t;       
        else
            pos = T(1:3, 4);  
        end
        
        P_min = [P_min; pos(:)'];
    end
end


figure('Name', 'Espacio de Trabajo - Máximos y Mínimos', 'NumberTitle', 'off');

%plano X-Y
subplot(1,3,1);
plot(P_max(:,1), P_max(:,2), 'r.', 'MarkerSize', 2); hold on;
plot(P_min(:,1), P_min(:,2), 'b.', 'MarkerSize', 2);
grid on; axis equal;
title('Plano X - Y ');
xlabel('X'); ylabel('Y');
legend('Límite Máximo', 'Límite Mínimo');

%plano X-Z
subplot(1,3,2);
plot(P_max(:,1), P_max(:,3), 'r.', 'MarkerSize', 2); hold on;
plot(P_min(:,1), P_min(:,3), 'b.', 'MarkerSize', 2);
grid on; axis equal;
title('Plano X - Z');
xlabel('X'); zlabel('Z');

%graficos
subplot(1,3,3);
plot3(P_max(:,1), P_max(:,2), P_max(:,3), 'r.', 'MarkerSize', 1); hold on;
plot3(P_min(:,1), P_min(:,2), P_min(:,3), 'b.', 'MarkerSize', 1);
grid on; axis equal;
title('Superficies 3D');
xlabel('X'); ylabel('Y'); zlabel('Z');