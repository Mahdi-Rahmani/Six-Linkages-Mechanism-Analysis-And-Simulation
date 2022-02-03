clc;
clear all;
close all;
% first of all we should write our data
w2 = 2*pi;
a2 = 0;
a = 0.6290;
% size of the linkages
AB = 0.0234094; BC = 0.04830114; BD = 0.02871411; CD = 0.04315669; EG = 0.0135;

%%%%%%%%%  PART1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% our unknown variables are AC, ED, teta3 , teta4
%first loop functions:
    %AB * cos(teta2) + BC * cos(180 - a -teta3) = AC;
    %AB * sin(teta2) = BC * sin(180 - a -teta3);
%second loop function:
    %CD * cos(teta3) = ED * cos(teta4) + CG;
    %CD * sin(teta3) = EG + ED * sin(teta4);

%we concider (AC, ED, teta3 , teta4) as (x(1),x(2),x(3),x(4))
%the initial value of teta2
teta2 = 1.2220;
%the initial values are achieved by solving equation with (teta2 = 70.0169 deg)
firstGuess = [0.051;0.04674666;2.0396;2.5773];
t = linspace(0,1,361);
ans1 = [0;0;0;0];    
i = 1;
while teta2 <= 7.5052
    F1 = @(x)[ AB * cos(teta2) + BC * cos(pi - a -x(3)) - x(1);
               AB * sin(teta2) - BC * sin(pi - a -x(3));
               CD * cos(x(3)) - x(2) * cos(x(4)) - 0.071 + x(1);
               CD * sin(x(3)) - EG - x(2) * sin(x(4))];
     x = fsolve(F1,firstGuess);
     ans1(:,i) = x;
     firstGuess = x;
     i = i + 1;
     teta2 = teta2 + 2*pi/360;
end

%now we should calculate the values related to first order diffrence
%equation(the velocity and angular velocity
%unknown variables are  Vac , Ved ,w3 , w4 ( v(1) v(2) v(3) v(4))
teta2 = 1.2220;
firstGuess = [-0.16394;0.11399;-1.169;0.9489];
ans2 = [0;0;0;0];    
k = 1;
while teta2 <= 7.5052
    F2 = @(v)[ -AB * w2 * sin(teta2) + BC * v(3) * sin(pi - a - ans1(3,k)) - v(1);
          AB * w2 * cos(teta2) + BC * v(3) * cos(pi - a -ans1(3,k));
          -CD * v(3) * sin(ans1(3,k)) - v(2) * cos(ans1(4,k)) + ans1(2,k)* v(4) *sin(ans1(4,k))+v(1) ;
          CD * v(3) * cos(ans1(3,k)) - v(2) * sin(ans1(4,k)) - ans1(2,k) * v(4) * cos(ans1(4,k))];
     v = fsolve(F2,firstGuess);
     ans2(:,k) = v;
     firstGuess = v;
     k = k + 1;
     teta2 = teta2 + 2*pi/360;
          
end

%now we should calculate the values related to second order diffrence
%equation(the acceleration and angular accelration)
%unknown variables are  Aac , Aed ,a3 , a4 
teta2 = 1.2220;
firstGuess = [0.07269;0.31979;19.4993;10.9065];
ans3 = [0;0;0;0];    
k = 1;
while teta2 <= 7.5052
    F3 = @(A)[ -AB * ((w2)^2) * cos(teta2) + BC *( A(3) * sin(pi - a - ans1(3,k)) - (ans2(3,k)^2)*cos(pi - a - ans1(3,k))) - A(1);
          -AB * ((w2)^2) * sin(teta2) + BC *( A(3) * cos(pi - a - ans1(3,k)) + (ans2(3,k)^2)*sin(pi - a - ans1(3,k)));
          -CD *(A(3) * sin(ans1(3,k)) + (ans2(3,k)^2) * cos(ans1(3,k))) - A(2)* cos(ans1(4,k))+ 2*ans2(2,k)*ans2(4,k) *sin(ans1(4,k))+ans1(2,k)*A(4)*sin(ans1(4,k)) + ans1(2,k)*(ans2(4,k)^2)*cos(ans1(4,k)) + A(1);
           CD * (A(3) * cos(ans1(3,k)) - (ans2(3,k)^2) * sin(ans1(3,k))) - A(2)* sin(ans1(4,k))- 2*ans2(2,k)*ans2(4,k) *cos(ans1(4,k))-ans1(2,k)*A(4)*cos(ans1(4,k)) + ans1(2,k)*(ans2(4,k)^2)*sin(ans1(4,k))];
     A = fsolve(F3,firstGuess);
     ans3(:,k) = A;
     firstGuess = A;
     k = k + 1;
     teta2 = teta2 + 2*pi/360;          
end

% plot linkage 2 (velocity and acceleration)
V_linkage2 = w2 *ones(1,361);
A_linkage2 = a2 *ones(1,361);
figure (1)
plot(t,A_linkage2)
hold on
plot(t,V_linkage2)
ylim([-1 8])
legend('acceleration','velocity');
title('linkage2')
xlabel('t')
grid on

% plot linkage 3 (velocity and acceleration)
figure (2)
plot(t,ans2(3,:))
hold on
plot(t,ans3(3,:))
legend('velocity','acceleration');
title('linkage3')
xlabel('t')
grid on

% plot linkage 4 (velocity and acceleration)
figure (3)
plot(t,ans2(4,:))
hold on
plot(t,ans3(4,:))
legend('velocity','acceleration');
title('linkage4')
xlabel('t')
grid on

% plot linkage 6(velocity and acceleration){AC}
figure (4)
plot(t,ans2(1,:))
hold on
plot(t,ans3(1,:))
legend('velocity','acceleration');
title('linkage6')
xlabel('t')
grid on

%%%%%%%%%  PART2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now we want to calculate the velocity and accelration of node 5 on
%floating linkage
%first we should find the Vx and Vy and also Ax and Ay
%then with pythagoras formula we can find the V and A for this floating
%point

teta2 = 1.2220;
VFloatingAns = zeros(1,361);
AFloatingAns = zeros(1,361);

VFloatingAnsX = zeros(1,361);
VFloatingAnsY = zeros(1,361);
AFloatingAnsX = zeros(1,361);
AFloatingAnsY = zeros(1,361);
k = 1;
while teta2 <= 7.5052
     Vx = ans2(2,k) * cos(ans1(4,k)) - ans1(2,k) * ans2(4,k) * sin(ans1(4,k));
     Vy = ans2(2,k) * sin(ans1(4,k)) + ans1(2,k) * ans2(4,k) * cos(ans1(4,k));
     VFloatingAnsX(k) = Vx;
     VFloatingAnsY(k) = Vy;
     VFloatingAns(k) = sqrt((Vx)^2 + (Vy)^2);
     
     Ax = ans3(2,k) * cos(ans1(4,k)) - 2*ans2(2,k) * ans2(4,k) * sin(ans1(4,k))- ans1(2,k) * ans3(4,k) * sin(ans1(4,k))- ans1(2,k) * ((ans2(4,k))^2) * cos(ans1(4,k));
     Ay = ans3(2,k) * sin(ans1(4,k)) + 2*ans2(2,k) * ans2(4,k) * cos(ans1(4,k))+ ans1(2,k) * ans3(4,k) * cos(ans1(4,k))- ans1(2,k) * ((ans2(4,k))^2) * sin(ans1(4,k));
     AFloatingAnsX(k) = Ax;
     AFloatingAnsY(k) = Ay;
     AFloatingAns(k) = sqrt((Ax)^2 + (Ay)^2);
     
     k = k + 1;
     teta2 = teta2 + 2*pi/360;          
end

% the acceleration and velocity in x and y directions
figure (5)
plot(t,VFloatingAnsX)
hold on
plot(t,VFloatingAnsY)
hold on
plot(t,AFloatingAnsX)
hold on
plot(t,AFloatingAnsY)
legend('velocity in x direction','velocity in y direction','acceleration in x direction','acceleration in y direction');
title('linkage5')
xlabel('t')
grid on

% the absolute acceleration and velocity 
figure (6)
plot(t,VFloatingAns)
hold on
plot(t,AFloatingAns)
legend('velocity','acceleration');
title('linkage5')
xlabel('t')
grid on


%%%%%%%%%  PART3  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now we want to do the force analysis
g = 9.81;
p = 100;
teta2 = 1.2220;
k = 1;
forceAnalysis = zeros(11,1);
%data for second linkage
f0 = 2 * (AB/2) * (w2)^2;
while teta2 <= 7.5052
    %data for third linkage
    m = 8;
    I3 = 0.02536;
    BG = 0.029;
    f1 = m * AB * (w2)^2;
    f2 = m * BG * (ans2(3,k))^2;
    f3 = m * BG * (ans3(3,k));
    gama = pi - ans1(3,k) - 0.8552;%B= 49 deg
    phy = 0.7854; % 45 deg
    %data for forth linkage
    l = 0.04674;
    I4 = 0.045;
    m4 = 4;  
%%%%%%%%%%%%%%%%%%%%%%  Ax  Ay   Bx  By   T   Fd  Ex   Ey   N   Cy  Cx 
    CoeficientMatrix = [1   ,0  ,-1 ,0   ,0  ,0  ,0   ,0   ,0  ,0  ,0;%
                        0   ,1  ,0  ,-1  ,0  ,0  ,0   ,0   ,0  ,0  ,0;%
                        0   ,0  ,-AB*sin(teta2),AB*cos(teta2),-1,0,0,0,0,0,0;%
                        0   ,0  ,1  ,0 ,0  ,-cos(ans1(4,k)-(pi/2)), 0, 0, 0, 0, -1;
                        0   ,0  ,0  ,1 ,0  ,-sin(ans1(4,k)-(pi/2)),0,0,0,1,0;
                        0   ,0  ,0  ,0 ,0  ,cos(ans1(4,k)-(pi/2))*BD*sin(pi-gama)-sin(ans1(4,k)-(pi/2))*BD*cos(phy-gama),0,0,0,BC*cos(pi-a-ans1(3,k)),-BC*sin(pi-a-ans1(3,k));
                        0   ,0  ,0  ,0 ,0  ,-cos(ans1(4,k)-(pi/2)),1,0,0,0,0;
                        0   ,0  ,0  ,0 ,0  ,sin(ans1(4,k)-(pi/2)),0,1,0,0,0;
                        0   ,0  ,0  ,0 ,0  ,ans1(2,k),0,0,0,0,0;
                        0   ,0  ,0  ,0 ,0  ,0,0,0,1,-1,0;
                        0   ,0  ,0  ,0 ,0  ,0,0,0,0,0,1];
     AnswerMatrix = [-(f0)*cos(teta2);
                     -(f0)*sin(teta2);
                      0;
                      (f3)*cos((pi/2)-gama) - (f2)*cos(gama)-(f1)*cos(teta2);
                      (f2)*sin(gama) - (f1)*sin(teta2) + (f3) * sin((pi/2)-gama);
                      -(f1)*sin(teta2)*BG + (f3)*BG + (I3)*ans3(3,k);
                      (m4)*(l/2)*ans3(4,k)*cos(ans1(4,k)-(pi/2)) - (m4)*(l/2)*(ans2(4,k)^2)*cos(pi-ans1(4,k));
                      -(m4)*(l/2)*ans3(4,k)*sin(ans1(4,k)-(pi/2)) - (m4)*(l/2)*(ans2(4,k)^2)*sin(pi-ans1(4,k));
                      -(l/2)*(m4)*(l/2)*ans3(4,k) - I4 * ans3(4,k);
                      0;
                      p];
    ANS = linsolve(CoeficientMatrix,AnswerMatrix);
    forceAnalysis(:,k) = ANS;
    k = k + 1;
    teta2 = teta2 + 2*pi/360;          
end
% plot Ax , Ay
figure (7)
plot(t,forceAnalysis(1,:))
hold on
plot(t,forceAnalysis(2,:))
legend('Ax','Ay');
title('Force analysis in A')
xlabel('t')
ylabel('N')
grid on

% plot Bx , By
figure (8)
plot(t,forceAnalysis(3,:))
hold on
plot(t,forceAnalysis(4,:))
legend('Bx','By');
title('Force analysis in B')
xlabel('t')
ylabel('N')
grid on

% plot T
figure (9)
plot(t,forceAnalysis(5,:))
title('Tork')
xlabel('t')
ylabel('N.m')
grid on

% plot FD
figure (10)
plot(t,forceAnalysis(6,:))
title('Force analysis in D')
xlabel('t')
ylabel('N')
grid on

% plot Ex , Ey
figure (11)
plot(t,forceAnalysis(7,:))
hold on
plot(t,forceAnalysis(8,:))
legend('Ex','Ey');
title('Force analysis in E')
xlabel('t')
ylabel('N')
grid on

% plot N
figure (12)
plot(t,forceAnalysis(9,:))
title('Force analysis (N)')
xlabel('t')
ylabel('N')
grid on

% plot Cy, Cx
figure (13)
plot(t,forceAnalysis(10,:))
hold on
plot(t,forceAnalysis(11,:))
ylim([-80 110])
legend('Cy','Cx');
title('Force analysis in C')
xlabel('t')
ylabel('N')
grid on



