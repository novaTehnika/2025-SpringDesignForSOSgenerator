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


% Vary the diameter of the piston
d = linspace(1,8,10)*0.0254; % [in -> m]
A = pi/4*d.^2;

% Solve for the change in force and failure

% define material
parMatl.E = 1e3; % [Pa]
parMatl.mu = 0.3;
parMatl.par.yeildStress; % [Pa]

% Determine dominate spring designs
fun = @(x) bellevilleOpFun(x,1);
nonlcon = @(x) bellevilleOpFun(x,2);
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

            varargout{1} = 
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
