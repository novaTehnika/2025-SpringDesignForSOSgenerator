% springDesign.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/19/2025
%
% PURPOSE/DESCRIPTION:
% This script ...
%
% FILE DEPENDENCY:
%
%
% UPDATES:
% See git log.
%
% Copyright (C) 2025  Jeremy W. Simmons II
% 
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program. If not, see <https://www.gnu.org/licenses/>.
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% define parameters of the diffusion process
P_f = 100*101.3e3; % [Pa] target pressure for dissolution
V_e = 0.5e-3; % [m^3] final mixture volume
v_prime = 0.05; % change in volume as a fraction of final mixture volume
lambda = 0; % margin between preload and target pressure

% Vary the diameter of the piston
d = linspace(1,8,10000)*0.0254; % [in -> m]
A = pi/4*d.^2;

deltay = v_prime*V_e./A; % displacement over dissolution phase
F_f = P_f*A;
F_p = (1-lambda)*P_f*A; % preload force

figure
plot(d/0.0254,deltay/0.0254)
xlabel('piston diameter (in)')
ylabel('displacement (in)')
title('Piston displacement over dissolution phase vs. piston diameter')

figure
plot(d/0.0254,F_f/4.44822162)
xlabel('piston diameter (in)')
ylabel('force (lbf)')
title('Spring force at target pressure vs. piston diameter')

%% assuming a change in force for a linear spring
% Vary the diameter of the piston
d = [2,3,4,5]*0.0254; % [in -> m]
A = pi/4*d.^2;

deltay = v_prime*V_e./A; % displacement over dissolution phase
F_f = P_f*A;
F_p = (1-lambda)*P_f*A; % preload force

%
maxDifflection = 0.25;
f_prime = logspace(log10(0.001),log10(0.5),100); % target fractional change in force during dissolution
F_max = zeros(numel(d),numel(f_prime));
k_linSpring = zeros(numel(d),numel(f_prime));
y_p = zeros(numel(d),numel(f_prime));
L = zeros(numel(d),numel(f_prime));

for i = 1:numel(d)
    F_max(i,:) = (1+f_prime)*F_f(i);
    k_linSpring(i,:) = f_prime*F_f(i)/deltay(i); % linear spring constant
    y_p(i,:) = F_p(i)./k_linSpring(i,:); % preload displacement of linear spring
    L(i,:) = (deltay(i)+y_p(i,:))/maxDifflection;
end

figure
clear leg
for i = 1:numel(d)
    loglog(f_prime,y_p(i,:)/0.0254)
    hold on
    leg(i) = {[num2str(d(i)/0.0254),'in piston']};
end
plot(f_prime([1,end]),deltay(i)*[1,1])
leg(numel(d)+1) = {'disp. during diss.'};
grid on
xlabel('change in force, fraction')
ylabel('preload displacement (in)')
title('Preload displacement vs. fraction change in force during dissolution')
legend(leg)

figure
clear leg
for i = 1:numel(d)
    loglog(f_prime,L(i,:)/0.0254)
    hold on
    leg(i) = {[num2str(d(i)/0.0254),'in piston']};
end
grid on
xlabel('change in force, fraction')
ylabel('unloaded length (in)')
title('Unloaded spring length vs. fraction change in force during dissolution')
legend(leg)

figure
clear leg
for i = 1:numel(d)
    loglog(f_prime,k_linSpring(i,:)/4.44822162/0.0254)
    hold on
    leg(i) = {[num2str(d(i)/0.0254),'in piston']};
end
grid on
xlabel('change in force, fraction')
ylabel('spring stiffness (lbf/in)')
title('Spring stiffness vs. fraction change in force during dissolution')
legend(leg)

%% Solve for the change in force and failure

% define material
parMatl.E = 200e9; % [Pa]
parMatl.mu = 0.3;
parMatl.par.yeildStress = 250e6; % [Pa]

% Determine dominate spring designs

 % Belleville spring
fun = @(x) bellevilleOpFun(x,1);
nonlcon = @(x) bellevilleOpFun(x,2);
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);

% plot change in force over the extent of travel

function varargout = bellevilleOpFun(x,deltay,parMatl,mode)
    % parse design variables
    par.D = x(1);
    par.d = x(2);
    par.t = x(3);
    par.H = x(4);
    par.E = parMatl.E;
    par.mu = parMatl.mu;
    parMatl.par.yeildStress;

    % find deflection at preload
    [F,~,~] = belleville(y,par);

    switch mode
        case 1 % objective function
            % obj 1: change in force over stroke

            % obj 2: installed length at preload

            varargout{1} = [];
        case 2 % nonlinear constraints
            % constraint 1: yeild stress not exceeded
             % stress at maximum difflection

             % compare to yeild stress
             varargout{1} = stressMax - parMatl.yeildStress;
            % constraint 2: maximum deflection not exceeded
            varargout{2} = y_p + deltay - yMax; 
    end

end

function [F,stress,s_max] = belleville(s,par)
% Model source: https://www.engineersedge.com/belleville_spring.htm
% (Accessed 19 June 2025)

    % calculation coeffs.
    delta = par.D/par.d;
    alpha = 1/pi*((delta-1)/delta)^2/((delta+1)/(delta-1)-2/log(delta));
    beta = 1/pi*6/log(delta)*((delta-1)/log(delta)-1);
    gamma = (delta-1)/pi*3/log(delta);
    h = par.H-par.t;
    
    
    % Common calc
    A = 4*par.E*par.t^4/((1-par.mu^2)*alpha*par.D^2);
    
    % Force
    F = A*s/par.t*((h/par.t-s/par.t)*(h/par.t-s/2/par.t)+1);
    
    % Stress
    stress = A*(beta*(h/par.t-s/2/par.t)+gamma);
    
    % Max displacement
    s_max = h;

end
