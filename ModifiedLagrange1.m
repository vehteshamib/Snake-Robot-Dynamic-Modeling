function dz = ModifiedLagrange1(t,z)

global m1 m2 m3 L1 L2 L3 I1 I2 I3

Tau1 = 0.01*sin(t);
Tau2 = 0.03*cos(10*t+pi/4);

dz = zeros(7,1);
q1 = z(1); q2 = z(2); q3 = z(3); q4 = z(4); q5 = z(5);
u1 = z(6); u2 = z(7);

M = [                m1 + m2 + m3,                          0,          -(L1*sin(q3)*(m2 + m3))/2,        -(L2*sin(q4)*(m2 + 2*m3))/2,        -(L3*m3*sin(q5))/2
    0,               m1 + m2 + m3,           (L1*cos(q3)*(m2 + m3))/2,         (L2*cos(q4)*(m2 + 2*m3))/2,         (L3*m3*cos(q5))/2
    -(L1*sin(q3)*(m2 + m3))/2,   (L1*cos(q3)*(m2 + m3))/2,     I1 + (L1^2*m2)/4 + (L1^2*m3)/4, (L1*L2*cos(q3 - q4)*(m2 + 2*m3))/4, (L1*L3*m3*cos(q3 - q5))/4
    -(L2*sin(q4)*(m2 + 2*m3))/2, (L2*cos(q4)*(m2 + 2*m3))/2, (L1*L2*cos(q3 - q4)*(m2 + 2*m3))/4,         I2 + (L2^2*m2)/4 + L2^2*m3, (L2*L3*m3*cos(q4 - q5))/2
    -(L3*m3*sin(q5))/2,          (L3*m3*cos(q5))/2,          (L1*L3*m3*cos(q3 - q5))/4,          (L2*L3*m3*cos(q4 - q5))/2,          (m3*L3^2)/4 + I3];

W1 = [                    cos(q3),                           0
    sin(q3),                           0
    0,                           1
    -(2*sin(q3 - q4))/L2,       -(L1*cos(q3 - q4))/L2
    (2*sin(q3 - 2*q4 + q5))/L3, (L1*cos(q3 - 2*q4 + q5))/L3];

dq = W1*[u1;u2];

dq1 = dq(1);    dq2 = dq(2);    dq3 = dq(3);     dq4 = dq(4);     dq5 = dq(5); 

F = [(L3*m3*cos(q5)*dq5^2)/2 + dq3*((L1*dq3*m2*cos(q3))/2 + (L1*dq3*m3*cos(q3))/2) + dq4*((L2*dq4*m2*cos(q4))/2 + L2*dq4*m3*cos(q4))
    (L3*m3*sin(q5)*dq5^2)/2 + dq3*((L1*dq3*m2*sin(q3))/2 + (L1*dq3*m3*sin(q3))/2) + dq4*((L2*dq4*m2*sin(q4))/2 + L2*dq4*m3*sin(q4))
    - Tau1 - (L1*L2*dq4^2*m2*sin(q3 - q4))/4 - (L1*L2*dq4^2*m3*sin(q3 - q4))/2 - (L1*L3*dq5^2*m3*sin(q3 - q5))/4
    Tau1 - Tau2 + (L1*L2*dq3^2*m2*sin(q3 - q4))/4 + (L1*L2*dq3^2*m3*sin(q3 - q4))/2 - (L2*L3*dq5^2*m3*sin(q4 - q5))/2
    Tau2 + (L3*m3*(L1*sin(q3 - q5)*dq3^2 + 2*L2*sin(q4 - q5)*dq4^2))/4];

W1_dot_u = [                                                         -dq3*u1*sin(q3)
    dq3*u1*cos(q3)
    0
    -((2*u1*cos(q3 - q4) - L1*u2*sin(q3 - q4))*(dq3 - dq4))/L2
    ((2*u1*cos(q3 - 2*q4 + q5) - L1*u2*sin(q3 - 2*q4 + q5))*(dq3 - 2*dq4 + dq5))/L3];


M_ML = W1' * M * W1;
F_ML = W1' * (F - M * W1_dot_u);

dz(1:5) = dq;
dz(6:7) = M_ML\F_ML;