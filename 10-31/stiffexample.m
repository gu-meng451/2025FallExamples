%% Ex: van der Pol's equation

clear all
close all
clc

%% Set the parameters
%
% mu is the only parameter of the system.  as mu increases interesting
% things happen.
mu = 5;

%%
% The equations of motion
f = @(t,y) [ y(2); mu*(1-y(1)^2)*y(2) - y(1)];
y0 = [2;0];
tf = mu*10;

%% ode45
sol1 = ode45(f,[0,tf],y0);
sol1.stats

plot(sol1.x,sol1.y(1,:), "DisplayName","ode45")
hold on

%% ode15s:
sol2 = ode15s(f,[0,tf],y0);
sol2.stats

plot(sol2.x,sol2.y(1,:), "DisplayName","ode15s")
