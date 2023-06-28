

clc;
clear;
close all;

% Vi =[];
% Vi(1)=0.3;
% Vi(2)=-0.3;
% Vi(3)=-0.4;
% for i = 4:1:27
%       Vi(i)=0.001;
% end

[t,x] = ode45(@NnAdaptive,[0,10],[0.3 -0.3 -0.4 0.001 0.001 0.001 0.001 0.001...
    0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001...
    0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001]);

figure(1)

plot(t,x(:,1:3));
grid on

gt = 1;
kt = 1;
s = 1;
m = 27;
fb = 0.26;

pmax = 20.185;
tb = sqrt(m)*fb/(2*pmax);

ka1 = 1.5; ka2 = 2; ku = 2.5;
mu = [-0.6 0 0.6]';

u = zeros (length(t),1);
fx = zeros(length(t),1);

for i = 1:length(t)
    phii = exp(-norm(x(i,1:3)'-mu)^2/(2*s^2));
    phi = phii*ones(27,1);
    pt = x(i,4:30)*phi;
    z1 = x(i,1);
    a1 = -ka1 * z1;
    z2 = x(i,2) - a1;
    a2 = -(ka1 + ka2)*z2 + z1*ka1^2;
    z3 = x(i,3)-a2;
    
    a2dot = -ka1^3*z1 + (ka1^2 + ka1*ka2 + ka2^2)*z2 - (ka1+ka2)*z3;
    zeta = 1e-6;
    db = 0.1;
     
    u(i) = -ku*z3 - tanh(z3*pt/zeta)*pt - tanh(z3*a2dot/zeta)*a2dot -...
        tanh(z3*fb/(2*zeta))*fb/2-tanh(z3*db/zeta)*db;
end

figure(2)
plot(t,u);
grid on

pt = zeros(length(t),1);
for i = 1:length(t)
    phii = exp(-norm(x(i,1:3)'-mu)^2/(2*s^2));
    phi = phii*ones(27,1);
    pt(i) = x(i,4:30)*phi;
    fx(i) = x(i,1)^2*x(i,3)^2 + x(i,2)^4;
end
figure(3)
plot(t,pt,t,fx);
grid on
 