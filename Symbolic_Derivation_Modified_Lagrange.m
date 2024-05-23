clc, clear all

syms q1 q2 q3 q4 q5 % q1=x, q2=y, q3=theta1, q4=theta2, q5=theta3
syms dq1 dq2 dq3 dq4 dq5 % dq/dt
syms Tau1 Tau2
syms m1 m2 m3 L1 L2 L3 I1 I2 I3
syms u1 u2
syms ddq1 ddq2 ddq3 ddq4 ddq5 

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
ddq = [ddq1; ddq2; ddq3 ; ddq4; ddq5];


dT_dt=0;
for i=1:5
    dT_dt = dT_dt + diff(T,q(i))*dq(i) + diff(T,dq(i))*ddq(i);
end
InputPower = Q.'*dq;
PowErr = simplify(dT_dt - InputPower)


constraints = simplify([V1.'*e1 ; V2.'*e2 ; V3.'*e3])

ConsErr = simplify(constraints.'*constraints) 

a = simplify(jacobian(constraints,dq))
a_dot = zeros(3,5);
for i=1:5
    a_dot=a_dot+diff(a,q(i))*dq(i);
end
a_dot_dq = simplify(a_dot*dq)

a1 = a(:,[2,4,5]);
a2 = a(:,[1,3]);

aa = -a1^(-1)*a2;
W = simplify([1 0;aa(1,:);0 1;aa(2:3,:)])

W_dot = zeros(5,2);
for i=1:5
    W_dot = W_dot+diff(W,q(i))*dq(i);
end
W_dot_dq = simplify(W_dot*[dq1;dq3])

% Lagrangian
L = T - V;


dL_dq = jacobian(L,q);
dL_ddq = jacobian(L,dq);

% Mass Matrix
M = simplify(jacobian (dL_ddq , dq))
B = simplify(jacobian (dL_ddq , q)*dq - dL_dq.');


F = Q - B

% Embedding
a1 = a(:,1:3);
a2 = a(:,4:5);

W = [-a1\a2;eye(2,2)];
W_dot = zeros(5,2);
for i=1:5
    W_dot = W_dot+diff(W,q(i))*dq(i);
end



% % First Modified Lagrange Method u1 = dx/cos(theta1); u2 = dtheta1;
% Y = [1/cos(q3) 0 0 0 0;0 0 1 0 0];
% A = [Y;a];
% B = [eye(2);zeros(3,2)];
% 
% W1 = simplify(A\B)
% W1_dot = zeros(5,2);
% for i=1:5
%     W1_dot = W1_dot+diff(W1,q(i))*dq(i);
% end
% u = [u1;u2];
% 
% W1_dot_u = simplify(W1_dot*u)
% 
% dq = W1*u;
% dq1 = dq(1); dq2 = dq(2);  dq3 = dq(3);  dq4 = dq(4); dq5 = dq(5); 
% F = simplify(eval(F));
% 
% 
% % Second Modified Lagrange Method 
% Y = [1/cos(q3) 0 0 0 0;0 0 0 1 0];
% A = [Y;a];
% B = [eye(2);zeros(3,2)];
% 
% W2 = simplify(A\B)
% W2_dot = zeros(5,2);
% for i=1:5
%     W2_dot = W2_dot+diff(W2,q(i))*dq(i);
% end
% u = [u1;u2];
% 
% W2_dot_u = simplify(W2_dot*u)
% 
% dq = W2*u;
% dq1 = dq(1); dq2 = dq(2);  dq3 = dq(3);  dq4 = dq(4); dq5 = dq(5); 
% F = simplify(eval(F))
% 
% % Third Modified Lagrange Method 
% Y = [1/cos(q3) 0 0 0 0;0 0 L1/2/sin(q4-q3) L2/2/sin(q4-q3)*cos(q4-q3) 0];
% A = [Y;a];
% B = [eye(2);zeros(3,2)];
% 
% W3 = simplify(A\B)
% W3_dot = zeros(5,2);
% for i=1:5
%     W3_dot = W3_dot+diff(W3,q(i))*dq(i);
% end
% u = [u1;u2];
% 
% W3_dot_u = simplify(W3_dot*u)
% 
% dq = W3*u;
% dq1 = dq(1); dq2 = dq(2);  dq3 = dq(3);  dq4 = dq(4); dq5 = dq(5); 
% F = simplify(eval(F))
% 
