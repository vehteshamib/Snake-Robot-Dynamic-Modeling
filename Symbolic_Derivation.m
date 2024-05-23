clc, clear all

syms q1 q2 q3 q4 q5 % q1=x, q2=y, q3=theta1, q4=theta2, q5=theta3
syms dq1 dq2 dq3 dq4 dq5 % dq/dt
syms Tau1 Tau2
syms m1 m2 m3 L1 L2 L3 I1 I2 I3

V1 = [dq1;dq2];
V2 = V1 + dq3*L1/2*[-sin(q3);cos(q3)] + dq4*L2/2*[-sin(q4);cos(q4)];
V3 = V2 + dq4*L2/2*[-sin(q4);cos(q4)] + dq5*L3/2*[-sin(q5);cos(q5)];

T = 1/2*m1*V1.'*V1 + 1/2 * I1 * dq3^2 + 1/2*m2*V2.'*V2 + 1/2 * I2 * dq4^2 + 1/2*m3*V3.'*V3 + 1/2 * I3 * dq5^2;
V = 0;
Q = [0;0;-Tau1;Tau1-Tau2;Tau2];

e1 = [-sin(q3);cos(q3)];
e2 = [-sin(q4);cos(q4)];
e3 = [-sin(q5);cos(q5)];

q = [q1; q2; q3 ; q4; q5];
dq = [dq1; dq2; dq3 ; dq4; dq5];

constraints = simplify([V1.'*e1 ; V2.'*e2 ; V3.'*e3])

a = simplify(jacobian(constraints,dq))
a_dot = zeros(3,5);
for i=1:5
    a_dot=a_dot+diff(a,q(i))*dq(i);
end
a_dot_dq = simplify(a_dot*dq)

a1 = a(:,1:3);
a2 = a(:,4:5);

W = [-a1^(-1)*a2;eye(2,2)];

W_dot = zeros(5,2);
for i=1:5
    W_dot = W_dot+diff(W,q(i))*dq(i);
end

% Lagrangian
L = T - V;

dL_dq = jacobian(L,q);
dL_ddq = jacobian(L,dq);

% Mass Matrix
M = simplify(jacobian (dL_ddq , dq))
B = simplify(jacobian (dL_ddq , q)*dq - dL_dq.');

F = Q - B





