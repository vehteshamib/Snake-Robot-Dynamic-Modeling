clc, clear all, close all

global m1 m2 m3 L1 L2 L3 I1 I2 I3
m1 = 1; m2 = 1.5; m3 = 1.3;
L1 = 1; L2 = 0.7; L3 = 2;
I1 = 0.1; I2 = 0.2; I3 = 0.3;


% options=odeset('maxstep',10^-2);
options=[];
timespan = [0:0.01:20];

t0=clock;
z0 = [0;0;pi/10;pi/5;pi/4;0;0;0;0;0;0;0;0];
[t,z] = ode45(@Int_Mul,timespan,z0,options);
q1IM = z(:,1); q2IM = z(:,2); q3IM = z(:,3); q4IM = z(:,4); q5IM = z(:,5);
dq1IM = z(:,6); dq2IM = z(:,7); dq3IM = z(:,8); dq4IM = z(:,9); dq5IM = z(:,10);
t1 = clock;
timeIM=etime(t1,t0)

Tau1 = 0.01*sin(t);
Tau2 = 0.03*cos(10*t+pi/4);

Len = length(t);

PowerErrorIM = t;
ConstraintErrorIM = t;
for i = 1:Len
    q1=z(i,1); q2=z(i,2); q3=z(i,3); q4=z(i,4);  q5=z(i,5);  dq1=z(i,6); dq2=z(i,7); dq3=z(i,8); dq4=z(i,9); dq5=z(i,10);
    A = Int_Mul(t(i),z(i,:)');
    ddq1 = A(6); ddq2 = A(7); ddq3 = A(8); ddq4 = A(9); ddq5 = A(10);
    
    ConstraintErrorIM(i) = ((L2*dq4)/2 + dq2*cos(q4) - dq1*sin(q4) + (L1*dq3*cos(q3 - q4))/2)^2 + ((L3*dq5)/2 + dq2*cos(q5) - dq1*sin(q5) + (L1*dq3*cos(q3 - q5))/2 + L2*dq4*cos(q4 - q5))^2 + (dq2*cos(q3) - dq1*sin(q3))^2;
    PowerErrorIM(i) = Tau1(i)*dq3 - Tau1(i)*dq4 + Tau2(i)*dq4 - Tau2(i)*dq5 + I1*ddq3*dq3 + I2*ddq4*dq4 + I3*ddq5*dq5 + ddq1*dq1*m1 + ddq2*dq2*m1 + ddq1*dq1*m2 + ddq2*dq2*m2 + ddq1*dq1*m3 + ddq2*dq2*m3 + (L1^2*ddq3*dq3*m2)/4 + (L1^2*ddq3*dq3*m3)/4 + (L2^2*ddq4*dq4*m2)/4 + L2^2*ddq4*dq4*m3 + (L3^2*ddq5*dq5*m3)/4 - (L1*dq1*dq3^2*m2*cos(q3))/2 - (L1*dq1*dq3^2*m3*cos(q3))/2 - (L2*dq1*dq4^2*m2*cos(q4))/2 - L2*dq1*dq4^2*m3*cos(q4) - (L3*dq1*dq5^2*m3*cos(q5))/2 - (L1*dq2*dq3^2*m2*sin(q3))/2 - (L1*dq2*dq3^2*m3*sin(q3))/2 - (L2*dq2*dq4^2*m2*sin(q4))/2 - L2*dq2*dq4^2*m3*sin(q4) - (L3*dq2*dq5^2*m3*sin(q5))/2 + (L1*ddq3*dq2*m2*cos(q3))/2 + (L1*ddq2*dq3*m2*cos(q3))/2 + (L1*ddq3*dq2*m3*cos(q3))/2 + (L1*ddq2*dq3*m3*cos(q3))/2 + (L2*ddq4*dq2*m2*cos(q4))/2 + (L2*ddq2*dq4*m2*cos(q4))/2 + L2*ddq4*dq2*m3*cos(q4) + L2*ddq2*dq4*m3*cos(q4) + (L3*ddq5*dq2*m3*cos(q5))/2 + (L3*ddq2*dq5*m3*cos(q5))/2 - (L1*ddq3*dq1*m2*sin(q3))/2 - (L1*ddq1*dq3*m2*sin(q3))/2 - (L1*ddq3*dq1*m3*sin(q3))/2 - (L1*ddq1*dq3*m3*sin(q3))/2 - (L2*ddq4*dq1*m2*sin(q4))/2 - (L2*ddq1*dq4*m2*sin(q4))/2 - L2*ddq4*dq1*m3*sin(q4) - L2*ddq1*dq4*m3*sin(q4) - (L3*ddq5*dq1*m3*sin(q5))/2 - (L3*ddq1*dq5*m3*sin(q5))/2 + (L1*L2*ddq4*dq3*m2*cos(q3 - q4))/4 + (L1*L2*ddq3*dq4*m2*cos(q3 - q4))/4 + (L1*L2*ddq4*dq3*m3*cos(q3 - q4))/2 + (L1*L2*ddq3*dq4*m3*cos(q3 - q4))/2 + (L1*L3*ddq5*dq3*m3*cos(q3 - q5))/4 + (L1*L3*ddq3*dq5*m3*cos(q3 - q5))/4 + (L2*L3*ddq5*dq4*m3*cos(q4 - q5))/2 + (L2*L3*ddq4*dq5*m3*cos(q4 - q5))/2 + (L1*L2*dq3*dq4^2*m2*sin(q3 - q4))/4 - (L1*L2*dq3^2*dq4*m2*sin(q3 - q4))/4 + (L1*L2*dq3*dq4^2*m3*sin(q3 - q4))/2 - (L1*L2*dq3^2*dq4*m3*sin(q3 - q4))/2 + (L1*L3*dq3*dq5^2*m3*sin(q3 - q5))/4 - (L1*L3*dq3^2*dq5*m3*sin(q3 - q5))/4 + (L2*L3*dq4*dq5^2*m3*sin(q4 - q5))/2 - (L2*L3*dq4^2*dq5*m3*sin(q4 - q5))/2;
end

TotalConstraintErrorIM = sum(ConstraintErrorIM)
TotalPowerErrorIM = sum(abs(PowerErrorIM))
disp('---------------------------------------------------------------------')

t0=clock;
z0 = [0;0;pi/10;pi/5;pi/4;0;0;0;0;0];
[t,z] = ode45(@Augmented,timespan,z0,options);
q1Au = z(:,1); q2Au = z(:,2); q3Au = z(:,3); q4Au = z(:,4); q5Au = z(:,5);
dq1Au = z(:,6); dq2Au = z(:,7); dq3Au = z(:,8); dq4Au = z(:,9); dq5Au = z(:,10);
t1 = clock;
timeAu=etime(t1,t0)

PowerErrorAu = t;
ConstraintErrorAu = t;
for i = 1:Len
    q1=z(i,1); q2=z(i,2); q3=z(i,3); q4=z(i,4);  q5=z(i,5);  dq1=z(i,6); dq2=z(i,7); dq3=z(i,8); dq4=z(i,9); dq5=z(i,10);
    A = Augmented(t(i),z(i,:)');
    ddq1 = A(6); ddq2 = A(7); ddq3 = A(8); ddq4 = A(9); ddq5 = A(10);
    
    ConstraintErrorAu(i)=((L2*dq4)/2 + dq2*cos(q4) - dq1*sin(q4) + (L1*dq3*cos(q3 - q4))/2)^2 + ((L3*dq5)/2 + dq2*cos(q5) - dq1*sin(q5) + (L1*dq3*cos(q3 - q5))/2 + L2*dq4*cos(q4 - q5))^2 + (dq2*cos(q3) - dq1*sin(q3))^2;
    PowerErrorAu(i) = Tau1(i)*dq3 - Tau1(i)*dq4 + Tau2(i)*dq4 - Tau2(i)*dq5 + I1*ddq3*dq3 + I2*ddq4*dq4 + I3*ddq5*dq5 + ddq1*dq1*m1 + ddq2*dq2*m1 + ddq1*dq1*m2 + ddq2*dq2*m2 + ddq1*dq1*m3 + ddq2*dq2*m3 + (L1^2*ddq3*dq3*m2)/4 + (L1^2*ddq3*dq3*m3)/4 + (L2^2*ddq4*dq4*m2)/4 + L2^2*ddq4*dq4*m3 + (L3^2*ddq5*dq5*m3)/4 - (L1*dq1*dq3^2*m2*cos(q3))/2 - (L1*dq1*dq3^2*m3*cos(q3))/2 - (L2*dq1*dq4^2*m2*cos(q4))/2 - L2*dq1*dq4^2*m3*cos(q4) - (L3*dq1*dq5^2*m3*cos(q5))/2 - (L1*dq2*dq3^2*m2*sin(q3))/2 - (L1*dq2*dq3^2*m3*sin(q3))/2 - (L2*dq2*dq4^2*m2*sin(q4))/2 - L2*dq2*dq4^2*m3*sin(q4) - (L3*dq2*dq5^2*m3*sin(q5))/2 + (L1*ddq3*dq2*m2*cos(q3))/2 + (L1*ddq2*dq3*m2*cos(q3))/2 + (L1*ddq3*dq2*m3*cos(q3))/2 + (L1*ddq2*dq3*m3*cos(q3))/2 + (L2*ddq4*dq2*m2*cos(q4))/2 + (L2*ddq2*dq4*m2*cos(q4))/2 + L2*ddq4*dq2*m3*cos(q4) + L2*ddq2*dq4*m3*cos(q4) + (L3*ddq5*dq2*m3*cos(q5))/2 + (L3*ddq2*dq5*m3*cos(q5))/2 - (L1*ddq3*dq1*m2*sin(q3))/2 - (L1*ddq1*dq3*m2*sin(q3))/2 - (L1*ddq3*dq1*m3*sin(q3))/2 - (L1*ddq1*dq3*m3*sin(q3))/2 - (L2*ddq4*dq1*m2*sin(q4))/2 - (L2*ddq1*dq4*m2*sin(q4))/2 - L2*ddq4*dq1*m3*sin(q4) - L2*ddq1*dq4*m3*sin(q4) - (L3*ddq5*dq1*m3*sin(q5))/2 - (L3*ddq1*dq5*m3*sin(q5))/2 + (L1*L2*ddq4*dq3*m2*cos(q3 - q4))/4 + (L1*L2*ddq3*dq4*m2*cos(q3 - q4))/4 + (L1*L2*ddq4*dq3*m3*cos(q3 - q4))/2 + (L1*L2*ddq3*dq4*m3*cos(q3 - q4))/2 + (L1*L3*ddq5*dq3*m3*cos(q3 - q5))/4 + (L1*L3*ddq3*dq5*m3*cos(q3 - q5))/4 + (L2*L3*ddq5*dq4*m3*cos(q4 - q5))/2 + (L2*L3*ddq4*dq5*m3*cos(q4 - q5))/2 + (L1*L2*dq3*dq4^2*m2*sin(q3 - q4))/4 - (L1*L2*dq3^2*dq4*m2*sin(q3 - q4))/4 + (L1*L2*dq3*dq4^2*m3*sin(q3 - q4))/2 - (L1*L2*dq3^2*dq4*m3*sin(q3 - q4))/2 + (L1*L3*dq3*dq5^2*m3*sin(q3 - q5))/4 - (L1*L3*dq3^2*dq5*m3*sin(q3 - q5))/4 + (L2*L3*dq4*dq5^2*m3*sin(q4 - q5))/2 - (L2*L3*dq4^2*dq5*m3*sin(q4 - q5))/2;
end
TotalConstraintErrorAu = sum(ConstraintErrorAu)
TotalPowerErrorAu = sum(abs(PowerErrorAu))
disp('---------------------------------------------------------------------')

t0=clock;
z0 = [0;0;pi/10;pi/5;pi/4;0;0;0;0;0];
[t,z] = ode45(@Elimination,timespan,z0,options);
q1El = z(:,1); q2El = z(:,2); q3El = z(:,3); q4El = z(:,4); q5El = z(:,5);
dq1El = z(:,6); dq2El = z(:,7); dq3El = z(:,8); dq4El = z(:,9); dq5El = z(:,10);
t1 = clock;
timeEl=etime(t1,t0)

PowerErrorEl = t;
ConstraintErrorEl = t;
for i = 1:Len
    q1=z(i,1); q2=z(i,2); q3=z(i,3); q4=z(i,4);  q5=z(i,5);  dq1=z(i,6); dq2=z(i,7); dq3=z(i,8); dq4=z(i,9); dq5=z(i,10);
    A = Elimination(t(i),z(i,:)');
    ddq1 = A(6); ddq2 = A(7); ddq3 = A(8); ddq4 = A(9); ddq5 = A(10);
    
    ConstraintErrorEl(i)=((L2*dq4)/2 + dq2*cos(q4) - dq1*sin(q4) + (L1*dq3*cos(q3 - q4))/2)^2 + ((L3*dq5)/2 + dq2*cos(q5) - dq1*sin(q5) + (L1*dq3*cos(q3 - q5))/2 + L2*dq4*cos(q4 - q5))^2 + (dq2*cos(q3) - dq1*sin(q3))^2;
    PowerErrorEl(i) = Tau1(i)*dq3 - Tau1(i)*dq4 + Tau2(i)*dq4 - Tau2(i)*dq5 + I1*ddq3*dq3 + I2*ddq4*dq4 + I3*ddq5*dq5 + ddq1*dq1*m1 + ddq2*dq2*m1 + ddq1*dq1*m2 + ddq2*dq2*m2 + ddq1*dq1*m3 + ddq2*dq2*m3 + (L1^2*ddq3*dq3*m2)/4 + (L1^2*ddq3*dq3*m3)/4 + (L2^2*ddq4*dq4*m2)/4 + L2^2*ddq4*dq4*m3 + (L3^2*ddq5*dq5*m3)/4 - (L1*dq1*dq3^2*m2*cos(q3))/2 - (L1*dq1*dq3^2*m3*cos(q3))/2 - (L2*dq1*dq4^2*m2*cos(q4))/2 - L2*dq1*dq4^2*m3*cos(q4) - (L3*dq1*dq5^2*m3*cos(q5))/2 - (L1*dq2*dq3^2*m2*sin(q3))/2 - (L1*dq2*dq3^2*m3*sin(q3))/2 - (L2*dq2*dq4^2*m2*sin(q4))/2 - L2*dq2*dq4^2*m3*sin(q4) - (L3*dq2*dq5^2*m3*sin(q5))/2 + (L1*ddq3*dq2*m2*cos(q3))/2 + (L1*ddq2*dq3*m2*cos(q3))/2 + (L1*ddq3*dq2*m3*cos(q3))/2 + (L1*ddq2*dq3*m3*cos(q3))/2 + (L2*ddq4*dq2*m2*cos(q4))/2 + (L2*ddq2*dq4*m2*cos(q4))/2 + L2*ddq4*dq2*m3*cos(q4) + L2*ddq2*dq4*m3*cos(q4) + (L3*ddq5*dq2*m3*cos(q5))/2 + (L3*ddq2*dq5*m3*cos(q5))/2 - (L1*ddq3*dq1*m2*sin(q3))/2 - (L1*ddq1*dq3*m2*sin(q3))/2 - (L1*ddq3*dq1*m3*sin(q3))/2 - (L1*ddq1*dq3*m3*sin(q3))/2 - (L2*ddq4*dq1*m2*sin(q4))/2 - (L2*ddq1*dq4*m2*sin(q4))/2 - L2*ddq4*dq1*m3*sin(q4) - L2*ddq1*dq4*m3*sin(q4) - (L3*ddq5*dq1*m3*sin(q5))/2 - (L3*ddq1*dq5*m3*sin(q5))/2 + (L1*L2*ddq4*dq3*m2*cos(q3 - q4))/4 + (L1*L2*ddq3*dq4*m2*cos(q3 - q4))/4 + (L1*L2*ddq4*dq3*m3*cos(q3 - q4))/2 + (L1*L2*ddq3*dq4*m3*cos(q3 - q4))/2 + (L1*L3*ddq5*dq3*m3*cos(q3 - q5))/4 + (L1*L3*ddq3*dq5*m3*cos(q3 - q5))/4 + (L2*L3*ddq5*dq4*m3*cos(q4 - q5))/2 + (L2*L3*ddq4*dq5*m3*cos(q4 - q5))/2 + (L1*L2*dq3*dq4^2*m2*sin(q3 - q4))/4 - (L1*L2*dq3^2*dq4*m2*sin(q3 - q4))/4 + (L1*L2*dq3*dq4^2*m3*sin(q3 - q4))/2 - (L1*L2*dq3^2*dq4*m3*sin(q3 - q4))/2 + (L1*L3*dq3*dq5^2*m3*sin(q3 - q5))/4 - (L1*L3*dq3^2*dq5*m3*sin(q3 - q5))/4 + (L2*L3*dq4*dq5^2*m3*sin(q4 - q5))/2 - (L2*L3*dq4^2*dq5*m3*sin(q4 - q5))/2;
end
TotalConstraintErrorEl = sum(ConstraintErrorEl)
TotalPowerErrorEl = sum(abs(PowerErrorEl))
disp('---------------------------------------------------------------------')

t0=clock;
z0 = [0;0;pi/10;pi/5;pi/4;0;0;0;0;0];
[t,z] = ode45(@GreenWood,timespan,z0,options);
q1GW = z(:,1); q2GW = z(:,2); q3GW = z(:,3); q4GW = z(:,4); q5GW = z(:,5);
dq1GW = z(:,6); dq2GW = z(:,7); dq3GW = z(:,8); dq4GW = z(:,9); dq5GW = z(:,10);
t1 = clock;
timeGW=etime(t1,t0)

PowerErrorGW = t;
ConstraintErrorGW=t;
for i = 1:Len
    q1=z(i,1); q2=z(i,2); q3=z(i,3); q4=z(i,4);  q5=z(i,5);  dq1=z(i,6); dq2=z(i,7); dq3=z(i,8); dq4=z(i,9); dq5=z(i,10);
    A = GreenWood(t(i),z(i,:)');
    ddq1 = A(6); ddq2 = A(7); ddq3 = A(8); ddq4 = A(9); ddq5 = A(10);
    
    ConstraintErrorGW(i)=((L2*dq4)/2 + dq2*cos(q4) - dq1*sin(q4) + (L1*dq3*cos(q3 - q4))/2)^2 + ((L3*dq5)/2 + dq2*cos(q5) - dq1*sin(q5) + (L1*dq3*cos(q3 - q5))/2 + L2*dq4*cos(q4 - q5))^2 + (dq2*cos(q3) - dq1*sin(q3))^2;
    PowerErrorGW(i) = Tau1(i)*dq3 - Tau1(i)*dq4 + Tau2(i)*dq4 - Tau2(i)*dq5 + I1*ddq3*dq3 + I2*ddq4*dq4 + I3*ddq5*dq5 + ddq1*dq1*m1 + ddq2*dq2*m1 + ddq1*dq1*m2 + ddq2*dq2*m2 + ddq1*dq1*m3 + ddq2*dq2*m3 + (L1^2*ddq3*dq3*m2)/4 + (L1^2*ddq3*dq3*m3)/4 + (L2^2*ddq4*dq4*m2)/4 + L2^2*ddq4*dq4*m3 + (L3^2*ddq5*dq5*m3)/4 - (L1*dq1*dq3^2*m2*cos(q3))/2 - (L1*dq1*dq3^2*m3*cos(q3))/2 - (L2*dq1*dq4^2*m2*cos(q4))/2 - L2*dq1*dq4^2*m3*cos(q4) - (L3*dq1*dq5^2*m3*cos(q5))/2 - (L1*dq2*dq3^2*m2*sin(q3))/2 - (L1*dq2*dq3^2*m3*sin(q3))/2 - (L2*dq2*dq4^2*m2*sin(q4))/2 - L2*dq2*dq4^2*m3*sin(q4) - (L3*dq2*dq5^2*m3*sin(q5))/2 + (L1*ddq3*dq2*m2*cos(q3))/2 + (L1*ddq2*dq3*m2*cos(q3))/2 + (L1*ddq3*dq2*m3*cos(q3))/2 + (L1*ddq2*dq3*m3*cos(q3))/2 + (L2*ddq4*dq2*m2*cos(q4))/2 + (L2*ddq2*dq4*m2*cos(q4))/2 + L2*ddq4*dq2*m3*cos(q4) + L2*ddq2*dq4*m3*cos(q4) + (L3*ddq5*dq2*m3*cos(q5))/2 + (L3*ddq2*dq5*m3*cos(q5))/2 - (L1*ddq3*dq1*m2*sin(q3))/2 - (L1*ddq1*dq3*m2*sin(q3))/2 - (L1*ddq3*dq1*m3*sin(q3))/2 - (L1*ddq1*dq3*m3*sin(q3))/2 - (L2*ddq4*dq1*m2*sin(q4))/2 - (L2*ddq1*dq4*m2*sin(q4))/2 - L2*ddq4*dq1*m3*sin(q4) - L2*ddq1*dq4*m3*sin(q4) - (L3*ddq5*dq1*m3*sin(q5))/2 - (L3*ddq1*dq5*m3*sin(q5))/2 + (L1*L2*ddq4*dq3*m2*cos(q3 - q4))/4 + (L1*L2*ddq3*dq4*m2*cos(q3 - q4))/4 + (L1*L2*ddq4*dq3*m3*cos(q3 - q4))/2 + (L1*L2*ddq3*dq4*m3*cos(q3 - q4))/2 + (L1*L3*ddq5*dq3*m3*cos(q3 - q5))/4 + (L1*L3*ddq3*dq5*m3*cos(q3 - q5))/4 + (L2*L3*ddq5*dq4*m3*cos(q4 - q5))/2 + (L2*L3*ddq4*dq5*m3*cos(q4 - q5))/2 + (L1*L2*dq3*dq4^2*m2*sin(q3 - q4))/4 - (L1*L2*dq3^2*dq4*m2*sin(q3 - q4))/4 + (L1*L2*dq3*dq4^2*m3*sin(q3 - q4))/2 - (L1*L2*dq3^2*dq4*m3*sin(q3 - q4))/2 + (L1*L3*dq3*dq5^2*m3*sin(q3 - q5))/4 - (L1*L3*dq3^2*dq5*m3*sin(q3 - q5))/4 + (L2*L3*dq4*dq5^2*m3*sin(q4 - q5))/2 - (L2*L3*dq4^2*dq5*m3*sin(q4 - q5))/2;
end
TotalConstraintErrorGW = sum(ConstraintErrorGW)
TotalPowerErrorGW = sum(abs(PowerErrorGW))
disp('---------------------------------------------------------------------')


t0=clock;
z0 = [0;0;pi/10;pi/5;pi/4;0;0];
[t,z] = ode45(@Embedding1,timespan,z0,options);
q1Em = z(:,1); q2Em = z(:,2); q3Em = z(:,3); q4Em = z(:,4); q5Em = z(:,5);
dq1Em = z(:,6); dq3Em = z(:,7);
t1 = clock;
timeEm=etime(t1,t0)

PowerErrorEm = t;
ConstraintErrorEm = t;
for i = 1:Len
    q1=z(i,1); q2=z(i,2); q3=z(i,3); q4=z(i,4);  q5=z(i,5);  dq1=z(i,6); dq3=z(i,7);
    
    W = [                                    1,                           0
        tan(q3),                           0
        0,                           1
        -(2*sin(q3 - q4))/(L2*cos(q3)),       -(L1*cos(q3 - q4))/L2
        (2*sin(q3 - 2*q4 + q5))/(L3*cos(q3)), (L1*cos(q3 - 2*q4 + q5))/L3];
    
    dq = W*[dq1;dq3];  dq2 = dq(2); dq4 = dq(4); dq5 = dq(5);
    
    A = Embedding1(t(i),z(i,:)');
    ddq1 = A(6); ddq3 = A(7);
    
    W_dot_dq = [                                                                                                                                                                                                                                                                       0
        dq1*dq3*(tan(q3)^2 + 1)
        0
        (L1*dq3^2*sin(q3 - q4) - L1*dq3*dq4*sin(q3 - q4))/L2 - (cos(q3)*(2*dq1*dq3*cos(q3 - q4) - 2*dq1*dq4*cos(q3 - q4)) + 2*dq1*dq3*sin(q3 - q4)*sin(q3))/(L2*cos(q3)^2)
        (cos(q3)*(2*dq1*dq3*cos(q3 - 2*q4 + q5) - 4*dq1*dq4*cos(q3 - 2*q4 + q5) + 2*dq1*dq5*cos(q3 - 2*q4 + q5)) + 2*dq1*dq3*sin(q3 - 2*q4 + q5)*sin(q3))/(L3*cos(q3)^2) - (L1*dq3^2*sin(q3 - 2*q4 + q5) - 2*L1*dq3*dq4*sin(q3 - 2*q4 + q5) + L1*dq3*dq5*sin(q3 - 2*q4 + q5))/L3];
    
    
    ddq = W*[ddq1;ddq3] + W_dot_dq;
    ddq2 = ddq(2); ddq4 = ddq(4);  ddq5 = ddq(5);
    
    ConstraintErrorEm(i)=((L2*dq4)/2 + dq2*cos(q4) - dq1*sin(q4) + (L1*dq3*cos(q3 - q4))/2)^2 + ((L3*dq5)/2 + dq2*cos(q5) - dq1*sin(q5) + (L1*dq3*cos(q3 - q5))/2 + L2*dq4*cos(q4 - q5))^2 + (dq2*cos(q3) - dq1*sin(q3))^2;
    PowerErrorEm(i) = Tau1(i)*dq3 - Tau1(i)*dq4 + Tau2(i)*dq4 - Tau2(i)*dq5 + I1*ddq3*dq3 + I2*ddq4*dq4 + I3*ddq5*dq5 + ddq1*dq1*m1 + ddq2*dq2*m1 + ddq1*dq1*m2 + ddq2*dq2*m2 + ddq1*dq1*m3 + ddq2*dq2*m3 + (L1^2*ddq3*dq3*m2)/4 + (L1^2*ddq3*dq3*m3)/4 + (L2^2*ddq4*dq4*m2)/4 + L2^2*ddq4*dq4*m3 + (L3^2*ddq5*dq5*m3)/4 - (L1*dq1*dq3^2*m2*cos(q3))/2 - (L1*dq1*dq3^2*m3*cos(q3))/2 - (L2*dq1*dq4^2*m2*cos(q4))/2 - L2*dq1*dq4^2*m3*cos(q4) - (L3*dq1*dq5^2*m3*cos(q5))/2 - (L1*dq2*dq3^2*m2*sin(q3))/2 - (L1*dq2*dq3^2*m3*sin(q3))/2 - (L2*dq2*dq4^2*m2*sin(q4))/2 - L2*dq2*dq4^2*m3*sin(q4) - (L3*dq2*dq5^2*m3*sin(q5))/2 + (L1*ddq3*dq2*m2*cos(q3))/2 + (L1*ddq2*dq3*m2*cos(q3))/2 + (L1*ddq3*dq2*m3*cos(q3))/2 + (L1*ddq2*dq3*m3*cos(q3))/2 + (L2*ddq4*dq2*m2*cos(q4))/2 + (L2*ddq2*dq4*m2*cos(q4))/2 + L2*ddq4*dq2*m3*cos(q4) + L2*ddq2*dq4*m3*cos(q4) + (L3*ddq5*dq2*m3*cos(q5))/2 + (L3*ddq2*dq5*m3*cos(q5))/2 - (L1*ddq3*dq1*m2*sin(q3))/2 - (L1*ddq1*dq3*m2*sin(q3))/2 - (L1*ddq3*dq1*m3*sin(q3))/2 - (L1*ddq1*dq3*m3*sin(q3))/2 - (L2*ddq4*dq1*m2*sin(q4))/2 - (L2*ddq1*dq4*m2*sin(q4))/2 - L2*ddq4*dq1*m3*sin(q4) - L2*ddq1*dq4*m3*sin(q4) - (L3*ddq5*dq1*m3*sin(q5))/2 - (L3*ddq1*dq5*m3*sin(q5))/2 + (L1*L2*ddq4*dq3*m2*cos(q3 - q4))/4 + (L1*L2*ddq3*dq4*m2*cos(q3 - q4))/4 + (L1*L2*ddq4*dq3*m3*cos(q3 - q4))/2 + (L1*L2*ddq3*dq4*m3*cos(q3 - q4))/2 + (L1*L3*ddq5*dq3*m3*cos(q3 - q5))/4 + (L1*L3*ddq3*dq5*m3*cos(q3 - q5))/4 + (L2*L3*ddq5*dq4*m3*cos(q4 - q5))/2 + (L2*L3*ddq4*dq5*m3*cos(q4 - q5))/2 + (L1*L2*dq3*dq4^2*m2*sin(q3 - q4))/4 - (L1*L2*dq3^2*dq4*m2*sin(q3 - q4))/4 + (L1*L2*dq3*dq4^2*m3*sin(q3 - q4))/2 - (L1*L2*dq3^2*dq4*m3*sin(q3 - q4))/2 + (L1*L3*dq3*dq5^2*m3*sin(q3 - q5))/4 - (L1*L3*dq3^2*dq5*m3*sin(q3 - q5))/4 + (L2*L3*dq4*dq5^2*m3*sin(q4 - q5))/2 - (L2*L3*dq4^2*dq5*m3*sin(q4 - q5))/2;
end
TotalConstraintErrorEm = sum(ConstraintErrorEm)
TotalPowerErrorEm = sum(abs(PowerErrorEm))
disp('---------------------------------------------------------------------')

t0=clock;
z0 = [0;0;pi/10;pi/5;pi/4;0;0];
[t,z] = ode45(@ModifiedLagrange1,timespan,z0,options);
q1ML1 = z(:,1); q2ML1 = z(:,2); q3ML1 = z(:,3); q4ML1 = z(:,4); q5ML1 = z(:,5);
u1 = z(:,6); u2 = z(:,7);
t1 = clock;
timeML1=etime(t1,t0)

PowerErrorML1 = t;
ConstraintErrorML1 = t;
for i = 1:Len
    q1=z(i,1); q2=z(i,2); q3=z(i,3); q4=z(i,4);  q5=z(i,5);  u1=z(i,6); u2=z(i,7);
    W1 = [                    cos(q3),                           0
        sin(q3),                           0
        0,                           1
        -(2*sin(q3 - q4))/L2,       -(L1*cos(q3 - q4))/L2
        (2*sin(q3 - 2*q4 + q5))/L3, (L1*cos(q3 - 2*q4 + q5))/L3];
    
    dq = W1 * [u1;u2];  dq1 = dq(1); dq2 = dq(2); dq3 = dq(3); dq4 = dq(4); dq5 = dq(5);
    
    A = ModifiedLagrange1(t(i),z(i,:)');
    du1 = A(6); du2 = A(7);
    
    W1_dot_u = [-dq3*u1*sin(q3)
        dq3*u1*cos(q3)
        0
        -((2*u1*cos(q3 - q4) - L1*u2*sin(q3 - q4))*(dq3 - dq4))/L2
        ((2*u1*cos(q3 - 2*q4 + q5) - L1*u2*sin(q3 - 2*q4 + q5))*(dq3 - 2*dq4 + dq5))/L3];
    
    ddq = W1*[du1;du2] + W1_dot_u;
    ddq1 = ddq(1); ddq2 = ddq(2);  ddq3 = ddq(3); ddq4 = ddq(4); ddq5 = ddq(5);
    
    ConstraintErrorML1(i)=((L2*dq4)/2 + dq2*cos(q4) - dq1*sin(q4) + (L1*dq3*cos(q3 - q4))/2)^2 + ((L3*dq5)/2 + dq2*cos(q5) - dq1*sin(q5) + (L1*dq3*cos(q3 - q5))/2 + L2*dq4*cos(q4 - q5))^2 + (dq2*cos(q3) - dq1*sin(q3))^2;
    PowerErrorML1(i) = Tau1(i)*dq3 - Tau1(i)*dq4 + Tau2(i)*dq4 - Tau2(i)*dq5 + I1*ddq3*dq3 + I2*ddq4*dq4 + I3*ddq5*dq5 + ddq1*dq1*m1 + ddq2*dq2*m1 + ddq1*dq1*m2 + ddq2*dq2*m2 + ddq1*dq1*m3 + ddq2*dq2*m3 + (L1^2*ddq3*dq3*m2)/4 + (L1^2*ddq3*dq3*m3)/4 + (L2^2*ddq4*dq4*m2)/4 + L2^2*ddq4*dq4*m3 + (L3^2*ddq5*dq5*m3)/4 - (L1*dq1*dq3^2*m2*cos(q3))/2 - (L1*dq1*dq3^2*m3*cos(q3))/2 - (L2*dq1*dq4^2*m2*cos(q4))/2 - L2*dq1*dq4^2*m3*cos(q4) - (L3*dq1*dq5^2*m3*cos(q5))/2 - (L1*dq2*dq3^2*m2*sin(q3))/2 - (L1*dq2*dq3^2*m3*sin(q3))/2 - (L2*dq2*dq4^2*m2*sin(q4))/2 - L2*dq2*dq4^2*m3*sin(q4) - (L3*dq2*dq5^2*m3*sin(q5))/2 + (L1*ddq3*dq2*m2*cos(q3))/2 + (L1*ddq2*dq3*m2*cos(q3))/2 + (L1*ddq3*dq2*m3*cos(q3))/2 + (L1*ddq2*dq3*m3*cos(q3))/2 + (L2*ddq4*dq2*m2*cos(q4))/2 + (L2*ddq2*dq4*m2*cos(q4))/2 + L2*ddq4*dq2*m3*cos(q4) + L2*ddq2*dq4*m3*cos(q4) + (L3*ddq5*dq2*m3*cos(q5))/2 + (L3*ddq2*dq5*m3*cos(q5))/2 - (L1*ddq3*dq1*m2*sin(q3))/2 - (L1*ddq1*dq3*m2*sin(q3))/2 - (L1*ddq3*dq1*m3*sin(q3))/2 - (L1*ddq1*dq3*m3*sin(q3))/2 - (L2*ddq4*dq1*m2*sin(q4))/2 - (L2*ddq1*dq4*m2*sin(q4))/2 - L2*ddq4*dq1*m3*sin(q4) - L2*ddq1*dq4*m3*sin(q4) - (L3*ddq5*dq1*m3*sin(q5))/2 - (L3*ddq1*dq5*m3*sin(q5))/2 + (L1*L2*ddq4*dq3*m2*cos(q3 - q4))/4 + (L1*L2*ddq3*dq4*m2*cos(q3 - q4))/4 + (L1*L2*ddq4*dq3*m3*cos(q3 - q4))/2 + (L1*L2*ddq3*dq4*m3*cos(q3 - q4))/2 + (L1*L3*ddq5*dq3*m3*cos(q3 - q5))/4 + (L1*L3*ddq3*dq5*m3*cos(q3 - q5))/4 + (L2*L3*ddq5*dq4*m3*cos(q4 - q5))/2 + (L2*L3*ddq4*dq5*m3*cos(q4 - q5))/2 + (L1*L2*dq3*dq4^2*m2*sin(q3 - q4))/4 - (L1*L2*dq3^2*dq4*m2*sin(q3 - q4))/4 + (L1*L2*dq3*dq4^2*m3*sin(q3 - q4))/2 - (L1*L2*dq3^2*dq4*m3*sin(q3 - q4))/2 + (L1*L3*dq3*dq5^2*m3*sin(q3 - q5))/4 - (L1*L3*dq3^2*dq5*m3*sin(q3 - q5))/4 + (L2*L3*dq4*dq5^2*m3*sin(q4 - q5))/2 - (L2*L3*dq4^2*dq5*m3*sin(q4 - q5))/2;
end
TotalConstraintErrorML1 = sum(ConstraintErrorML1)
TotalPowerErrorML1 = sum(abs(PowerErrorML1))
disp('---------------------------------------------------------------------')


figure(1)
hold on
plot(t,q1IM,'y-','linewidth',10)
plot(t,q1Au,'g--','linewidth',8)
plot(t,q1El,'r-','linewidth',6)
plot(t,q1GW,'b--','linewidth',4)
plot(t,q1Em,'m-','linewidth',2)
plot(t,q1ML1,'k--','linewidth',1)

legend('Integrated Multiplier Method','Augmented Method','Elimination Method','GreenWood','Embedding Method','Modified Lagrange')

set(gca,'fontsize',18,'fontweight','bold');
xlabel('Time (s)','fontsize',25,'fontweight','bold');
ylabel('X (m)','fontsize',25,'fontweight','bold');


figure(2)
hold on
plot(t,ConstraintErrorIM,'y-','linewidth',3)
plot(t,ConstraintErrorAu,'g--','linewidth',2)
plot(t,ConstraintErrorEl,'r-','linewidth',2)
plot(t,ConstraintErrorGW,'b--','linewidth',2)
plot(t,ConstraintErrorEm,'m-','linewidth',2)
plot(t,ConstraintErrorML1,'k--','linewidth',2)

legend('Integrated Multiplier Method','Augmented Method','Elimination Method','GreenWood','Embedding Method','Modified Lagrange')

set(gca,'fontsize',18,'fontweight','bold');
xlabel('Time (s)','fontsize',25,'fontweight','bold');
ylabel('Constraints Error','fontsize',25,'fontweight','bold')


figure(3)
hold on
plot(t,PowerErrorIM,'y-','linewidth',3)
plot(t,PowerErrorAu,'g--','linewidth',2)
plot(t,PowerErrorEl,'r-','linewidth',2)
plot(t,PowerErrorGW,'b--','linewidth',2)
plot(t,PowerErrorEm,'m-','linewidth',2)
plot(t,PowerErrorML1,'k--','linewidth',2)

legend('Integrated Multiplier Method','Augmented Method','Elimination Method','GreenWood','Embedding Method','Modified Lagrange')

set(gca,'fontsize',18,'fontweight','bold');
xlabel('Time (s)','fontsize',25,'fontweight','bold');
ylabel('Power Error (Watt)','fontsize',25,'fontweight','bold');

SumPowerErrorIM = PowerErrorIM;
SumPowerErrorAu = PowerErrorAu;
SumPowerErrorEl = PowerErrorEl;
SumPowerErrorGW = PowerErrorGW;
SumPowerErrorEm = PowerErrorEm;
SumPowerErrorML1 = PowerErrorML1;

for i = 1:Len
    SumPowerErrorIM(i) = sum(abs(PowerErrorIM(1:i)))/i;
    SumPowerErrorAu(i) = sum(abs(PowerErrorAu(1:i)))/i;
    SumPowerErrorEl(i) = sum(abs(PowerErrorEl(1:i)))/i;
    SumPowerErrorGW(i) = sum(abs(PowerErrorGW(1:i)))/i;
    SumPowerErrorEm(i) = sum(abs(PowerErrorEm(1:i)))/i;
    SumPowerErrorML1(i) = sum(abs(PowerErrorML1(1:i)))/i;
end

figure(4)
hold on
plot(t,SumPowerErrorIM,'y-','linewidth',3)
plot(t,SumPowerErrorAu,'g--','linewidth',2)
plot(t,SumPowerErrorEl,'r-','linewidth',2)
plot(t,SumPowerErrorGW,'b--','linewidth',2)
plot(t,SumPowerErrorEm,'m-','linewidth',2)
plot(t,SumPowerErrorML1,'k--','linewidth',2)

legend('Integrated Multiplier Method','Augmented Method','Elimination Method','GreenWood','Embedding Method','Modified Lagrange')

set(gca,'fontsize',18,'fontweight','bold');
xlabel('Time (s)','fontsize',25,'fontweight','bold');
ylabel('Power Error (Watt)','fontsize',25,'fontweight','bold');

% figure(5)
% hold on
% plot(t,q1ML1 - q1IM,'y-','linewidth',2)
% plot(t,q1ML1 - q1Au,'g--','linewidth',2)
% plot(t,q1ML1 - q1El,'r-','linewidth',2)
% plot(t,q1ML1 - q1GW,'b--','linewidth',2)
% plot(t,q1ML1 - q1Em,'m-','linewidth',2)
% 
% legend('Integrated Multiplier Method','Augmented Method','Elimination Method','GreenWood','Embedding Method')
% set(gca,'fontsize',18,'fontweight','bold');
% xlabel('Time (s)','fontsize',25,'fontweight','bold');
% ylabel('Difference With Modified Lagrange Method','fontsize',25,'fontweight','bold');
% 
