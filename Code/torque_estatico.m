%%Torque estático

m_es = 0.1; %kg
m_tr = 0.05;
m_gr = 0.05;
m_v = 0.2;

d1 = 0.15; %m
d2 = 0.15;
d3 = 0.05;

g = 9.81;

T2 = ((m_v + m_gr))*g*((d2 + d3)) + m_es*g*d2/2; 

T1 = (m_v + m_gr)*g*(d2 + d3 + d1) + m_es*g*(d2/2+d1+d1/2) + m_tr*g*d1;

%[N*m]*[100 cm/1m] 

T_2 = T2 * 100 %N*cm
T_1 = T1 * 100