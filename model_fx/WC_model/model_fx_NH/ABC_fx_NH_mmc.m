function dX = ABC_fx_NH_mmc(x,P,D)
% state equations for a neural mass model of motor cortex
% Bhatt et al. 2016 Neuroimage
%
% FORMAT [f,J,D] = spm_fx_mmc(x,u,P,M)
% FORMAT [f,J]   = spm_fx_mmc(x,u,P,M)
% FORMAT [f]     = spm_fx_mmc(x,u,P,M)
% x      - state vector
%   x(:,1) - rate     (middle pyramidal cells)
%   x(:,2) - rate     (superficial pyramidal cells)
%   x(:,3) - rate     (inhibitory interneurons)
%   x(:,4) - rate     (deep pyramidal cells)
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
% pre-synaptic inputs: s(V)
%--------------------------------------------------------------------------
% input
%==========================================================================

% time constants and intrinsic connections
%==========================================================================
T = P.T;
G = P.G;
Sn = P.Sn;
Mn = P.Mn;
Bn = P.Bn;
% Br = P.Br;
dX(:,1) = (-x(end,1) + sigNevHol(-G(1).*x(end-D(1),1) - G(3).*x(end-D(3),4) + G(8).*x(end-D(8),2),Mn(1),Sn(1),Bn(1)) )./T(1);
dX(:,2) = (-x(end,2) + sigNevHol(-G(7).*x(end-D(7),2) + G(2).*x(end-D(2),1) - G(12).*x(end-D(12),3) + G(14).*x(end-D(7),4),Mn(2),Sn(2),Bn(2)) )./T(2) ;
dX(:,3) = (-x(end,3) + sigNevHol(-G(4).*x(end-D(4),3) + G(5).*x(end-D(5),1) + G(6).*x(end-D(6),4) + G(13).*x(end-D(13),2),Mn(3),Sn(3),Bn(3)) )./T(3) ;
dX(:,4) = (-x(end,4) + sigNevHol(-G(10).*x(end-D(10),4) - G(9).*x(end-D(9),3) + G(11).*x(end-D(11),2),Mn(4),Sn(4),Bn(4)) )./T(4) ;

