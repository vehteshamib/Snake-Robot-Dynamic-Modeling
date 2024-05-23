function dz = ModifiedLagrange(t,z)

global m1 m2 m3 L1 L2 L3 I1 I2 I3

Tau1 = 0.01*sin(t);
Tau2 = 0.03*cos(10*t+pi/4);

dz = zeros(7,1);
q1 = z(1); q2 = z(2); q3 = z(3); q4 = z(4); q5 = z(5);
u1 = z(6); u2 = z(7);



a = [ -sin(q3), cos(q3),                                               0,                                       0,                                   0
    -sin(q4), cos(q4), (L1*cos(q3)*cos(q4))/2 + (L1*sin(q3)*sin(q4))/2,     (L2*cos(q4)^2)/2 + (L2*sin(q4)^2)/2,                                   0
    -sin(q5), cos(q5), (L1*cos(q3)*cos(q5))/2 + (L1*sin(q3)*sin(q5))/2, L2*cos(q4)*cos(q5) + L2*sin(q4)*sin(q5), (L3*cos(q5)^2)/2 + (L3*sin(q5)^2)/2];


M = [                m1 + m2 + m3,                          0,          -(L1*sin(q3)*(m2 + m3))/2,        -(L2*sin(q4)*(m2 + 2*m3))/2,        -(L3*m3*sin(q5))/2
    0,               m1 + m2 + m3,           (L1*cos(q3)*(m2 + m3))/2,         (L2*cos(q4)*(m2 + 2*m3))/2,         (L3*m3*cos(q5))/2
    -(L1*sin(q3)*(m2 + m3))/2,   (L1*cos(q3)*(m2 + m3))/2,     I1 + (L1^2*m2)/4 + (L1^2*m3)/4, (L1*L2*cos(q3 - q4)*(m2 + 2*m3))/4, (L1*L3*m3*cos(q3 - q5))/4
    -(L2*sin(q4)*(m2 + 2*m3))/2, (L2*cos(q4)*(m2 + 2*m3))/2, (L1*L2*cos(q3 - q4)*(m2 + 2*m3))/4,         I2 + (L2^2*m2)/4 + L2^2*m3, (L2*L3*m3*cos(q4 - q5))/2
    -(L3*m3*sin(q5))/2,          (L3*m3*cos(q5))/2,          (L1*L3*m3*cos(q3 - q5))/4,          (L2*L3*m3*cos(q4 - q5))/2,          (m3*L3^2)/4 + I3];

a1 = a(:,1:3);
a2 = a(:,4:5);

W = [-a1\a2;eye(2,2)];

dq = W*[dq4;dq5];

dq1=dq(1); dq2=dq(2); dq3=dq(3);

a_dot_dq = [                                                                                                                -dq3*(dq1*cos(q3) + dq2*sin(q3))
    - dq3*((L1*dq3*sin(q3 - q4))/2 - (L1*dq4*sin(q3 - q4))/2) - dq2*dq4*sin(q4) - dq1*dq4*cos(q4)
    - dq3*((L1*dq3*sin(q3 - q5))/2 - (L1*dq5*sin(q3 - q5))/2) - dq4*(L2*dq4*sin(q4 - q5) - L2*dq5*sin(q4 - q5)) - dq2*dq5*sin(q5) - dq1*dq5*cos(q5)];

F = [(L3*m3*cos(q5)*dq5^2)/2 + dq3*((L1*dq3*m2*cos(q3))/2 + (L1*dq3*m3*cos(q3))/2) + dq4*((L2*dq4*m2*cos(q4))/2 + L2*dq4*m3*cos(q4))
    (L3*m3*sin(q5)*dq5^2)/2 + dq3*((L1*dq3*m2*sin(q3))/2 + (L1*dq3*m3*sin(q3))/2) + dq4*((L2*dq4*m2*sin(q4))/2 + L2*dq4*m3*sin(q4))
    - Tau1 - (L1*L2*dq4^2*m2*sin(q3 - q4))/4 - (L1*L2*dq4^2*m3*sin(q3 - q4))/2 - (L1*L3*dq5^2*m3*sin(q3 - q5))/4
    Tau1 - Tau2 + (L1*L2*dq3^2*m2*sin(q3 - q4))/4 + (L1*L2*dq3^2*m3*sin(q3 - q4))/2 - (L2*L3*dq5^2*m3*sin(q4 - q5))/2
    Tau2 + (L3*m3*(L1*sin(q3 - q5)*dq3^2 + 2*L2*sin(q4 - q5)*dq4^2))/4];

M1 = M(1:3,:);
M2 = M(4:5,:);

F1 = F(1:3);
F2 = F(4:5);


W_dot_dq = [  -(3*L2*dq4^2*sin(q5) + 4*L2*dq4^2*sin(2*q3 - 2*q4 + q5) + L3*dq5^2*sin(2*q3 + q4 - 2*q5) + L3*dq5^2*sin(q4 - 2*q5) - L2*dq4^2*sin(2*q3 - 4*q4 + 3*q5) - L3*dq5^2*sin(2*q3 - 3*q4 + 2*q5) - 3*L2*dq4^2*sin(2*q3 - q5) - 4*L2*dq4^2*sin(2*q4 - q5) + L2*dq4^2*sin(4*q4 - 3*q5) + L3*dq5^2*sin(3*q4 - 2*q5) - 2*L2*dq4*dq5*sin(q5) - 2*L3*dq4*dq5*sin(q4) - 4*L2*dq3*dq4*sin(2*q3 - 2*q4 + q5) - 2*L2*dq4*dq5*sin(2*q3 - 2*q4 + q5) + 2*L3*dq3*dq5*sin(2*q3 + q4 - 2*q5) - 2*L3*dq4*dq5*sin(2*q3 + q4 - 2*q5) - 2*L3*dq4*dq5*sin(q4 - 2*q5) + 2*L2*dq3*dq4*sin(2*q3 - 4*q4 + 3*q5) + 2*L3*dq3*dq5*sin(2*q3 - 3*q4 + 2*q5) + 2*L2*dq3*dq4*sin(2*q3 - q5) - 4*L3*dq3*dq5*sin(2*q3 - q4) + 2*L2*dq4*dq5*sin(2*q3 - q5) + 2*L3*dq4*dq5*sin(2*q3 - q4) + 2*L2*dq4*dq5*sin(2*q4 - q5))/(12*sin(q4 - q5) - 4*sin(3*q4 - 3*q5))
    -(8*L3*dq5^2*sin(q5) - 24*L2*dq4^2*sin(q4) + 8*L3*dq5^2*sin(2*q3 - 2*q4 + q5) + L2*dq4^2*sin(q4 - 2*q5) + 8*L2*dq4^2*sin(2*q3 - 3*q4 + 2*q5) - 24*L2*dq4^2*sin(2*q3 - q4) + L2*dq4^2*sin(2*q5 - q4) + 8*L2*dq4^2*sin(3*q4 - 2*q5) + 8*L3*dq5^2*sin(2*q3 - q5) + 8*L3*dq5^2*sin(2*q4 - q5) + 16*L2*dq4*dq5*sin(q4) - 16*L3*dq4*dq5*sin(q5) - 16*L3*dq3*dq5*sin(2*q3 - 2*q4 + q5) - L2*dq4*dq5*sin(q4 - 2*q5) - 16*L2*dq3*dq4*sin(2*q3 - 3*q4 + 2*q5) + 16*L2*dq3*dq4*sin(2*q3 - q4) + 16*L2*dq4*dq5*sin(2*q3 - q4) + 16*L3*dq3*dq5*sin(2*q3 - q5) - L2*dq4*dq5*sin(2*q5 - q4) - 16*L3*dq4*dq5*sin(2*q3 - q5))/(64*sin(q4 - q5)^2)
    (L2*dq4^2*sin(q3 - 3*q4 + 2*q5) - 3*L2*dq4^2*sin(q3 - q4) + L3*dq5^2*sin(q3 - q5) + L3*dq5^2*sin(q3 - 2*q4 + q5) - L2*dq3*dq4*sin(q3 - 3*q4 + 2*q5) + L2*dq3*dq4*sin(q3 - q4) + 2*L2*dq4*dq5*sin(q3 - q4) + L3*dq3*dq5*sin(q3 - q5) - 2*L3*dq4*dq5*sin(q3 - q5) - L3*dq3*dq5*sin(q3 - 2*q4 + q5))/(L1 - L1*cos(2*q4 - 2*q5))
    0
    0];

M_Em = W' * M * W;
F_Em = W' * (F - M * W_dot_dq);

dz(1:5) = dq;
dz(6:7) = M_Em\F_Em;