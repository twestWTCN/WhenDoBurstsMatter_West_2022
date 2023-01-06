function f = ABC_fx_JR_mmc(x,P,D)
% state equations for a neural mass model of motor cortex
% Bhatt et al. 2016 Neuroimage
%
% FORMAT [f,J,D] = spm_fx_mmc(x,u,P,M)
% FORMAT [f,J]   = spm_fx_mmc(x,u,P,M)
% FORMAT [f]     = spm_fx_mmc(x,u,P,M)
% x      - state vector
%   x(:,1) - voltage     (middle pyramidal cells)
%   x(:,2) - conductance (middle pyramdidal cells)
%   x(:,3) - voltage     (superficial pyramidal cells)
%   x(:,4) - conductance (superficial pyramidal cells)
%   x(:,5) - current     (inhibitory interneurons)
%   x(:,6) - conductance (inhibitory interneurons)
%   x(:,7) - voltage     (deep pyramidal cells)
%   x(:,8) - conductance (deep pyramidal cells)
%
% f        - dx(t)/dt  = f(x(t))
% J        - df(t)/dx(t)
% D        - delay operator dx(t)/dt = f(x(t - d))
%                                    = D(d)*f(x(t))
%
% Prior fixed parameter scaling [Defaults]
%
% E  = (forward, backward, lateral) extrinsic rates 
% G  = intrinsic rates
% D  = propagation delays (intrinsic, extrinsic)
% T  = synaptic time constants
% S  = slope of sigmoid activation function
%
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging


% pre-synaptic inputs: s(V)
%--------------------------------------------------------------------------
R    = P.Rz(2:end);              % gain of activation function (1st is extrinsic- so remove)
S = sigmoidin(x,R,0);
% R    = (2/3);     %0.5.*                  % slope of sigmoid activation function
% B    = 0;                        % bias or background (sigmoid)
% R    = R.*exp(P.S);              % gain of activation function
% F    = 1./(1 + exp(-R*x + B));   % firing rate
% S    = F - 1/(1 + exp(B));       % deviation from baseline firing

% input
%==========================================================================

% time constants and intrinsic connections
%==========================================================================
% T    = ones(n,1)*T/1000;
% G    = ones(n,1)*G;
T = P.T;
G = P.G;
% C = P.C;
% extrinsic connections
%--------------------------------------------------------------------------
% forward  (i)   2  sp -> mp (+ve)
% forward  (ii)  1  sp -> dp (+ve)
% backward (i)   2  dp -> sp (-ve)
% backward (ii)  1  dp -> ii (-ve)
%--------------------------------------------------------------------------
% free parameters on time constants and intrinsic connections
%--------------------------------------------------------------------------
% G(:,1)  mp -> mp (-ve self)  4
% G(:,2)  mp -> sp (+ve rec )  4
% G(:,3)  ii -> mp (-ve rec )  4
% G(:,4)  ii -> ii (-ve self)  4
% G(:,5)  mp -> ii (+ve rec )  4
% G(:,6)  dp -> ii (+ve rec )  2
% G(:,7)  sp -> sp (-ve self)  4
% G(:,8)  sp -> mp (+ve rec )  4
% G(:,9)  ii -> dp (-ve rec )  2
% G(:,10) dp -> dp (-ve self)  1
% G(:,11) sp -> dp (+ve rec)  2
% G(:,12) ii -> sp (-ve rec)  4
% G(:,13) sp -> ii (+ve rec)  4
% G(:,14) dp -> sp (+ve rec)  2
%--------------------------------------------------------------------------
% Neuronal states (deviations from baseline firing)
%--------------------------------------------------------------------------
%   S(:,1) - voltage     (middle pyramidal cells)
%   S(:,2) - conductance (middle pyramidal cells)
%   S(:,3) - voltage     (superficial pyramidal cells)
%   S(:,4) - conductance (superficial pyramidal cells)
%   S(:,5) - current     (inhibitory interneurons)
%   S(:,6) - conductance (inhibitory interneurons)
%   S(:,7) - voltage     (deep pyramidal cells)
%   S(:,8) - conductance (deep pyramidal cells)
%--------------------------------------------------------------------------
 
% Motion of states: f(x)
%--------------------------------------------------------------------------
 
% Conductance
%==========================================================================
 
% Middle layer (middle pyramidal cells): Hidden causes
%--------------------------------------------------------------------------
u      =  0; 
u      = - G(1).*S(end-D(1),1) - G(:,3).*S(end-D(3),5) + G(:,8).*S(end-D(8),3) + u;
f(:,2) =  (u - 2*x(end,2) - x(end,1)./T(:,1))./T(:,1);
 
% Supra-granular layer (superficial pyramidal cells): Hidden causes - error
%--------------------------------------------------------------------------
u      = 0; 
u      =  - G(:,7).*S(end-D(7),3) + G(:,2).*S(end-D(2),1) - G(:,12).*S(end-D(12),5) + G(:,14).*S(end-D(14),7) + u;
f(:,4) =  (u - 2*x(end,4) - x(end,3)./T(:,2))./T(:,2);
 
% Supra-granular layer (inhibitory interneurons): Hidden states - error
%--------------------------------------------------------------------------
u      = 0; 
u      =  - G(:,4).*S(end-D(4),5) + G(:,5).*S(end-D(5),1) + G(:,6).*S(end-D(6),7) + G(:,13).*S(end-D(13),3) + u;
f(:,6) =  (u - 2*x(end,6) - x(end,5)./T(:,3))./T(:,3);
 
% Infra-granular layer (deep pyramidal cells): Hidden states
%--------------------------------------------------------------------------
u      =  0; %
u      = - G(:,10).*S(end-D(10),7) - G(:,9).*S(end-D(9),5) + G(:,11).*S(end-D(11),3) + u;
f(:,8) =  (u - 2*x(end,8) - x(end,7)./T(:,4))./T(:,4);
 
% Voltage
%==========================================================================
f(:,1) = x(end,2);
f(:,3) = x(end,4);
f(:,5) = x(end,6);
f(:,7) = x(end,8);
f = f';
% f      = spm_vec(f);
 
