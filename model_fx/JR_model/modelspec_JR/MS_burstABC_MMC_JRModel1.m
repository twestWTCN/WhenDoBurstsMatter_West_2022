function [R p m uc] = MS_burstABC_MMC_JRModel1(R)
% THIS IS THE STN/GPE
% Null-Model
[R,m] = getStateDetailsJR(R);
%   m.outstates{1} = [0 0 1 0 0 0 0 0];
% setup exogenous noise
% m.uset.p = DCM.Ep;
for nm = 1:m.m
    m.uset.p.covar{nm} = eye(m.Cint(nm));
    m.uset.p.scale(nm) = 1;
    m.uset.p.alpha{nm} = repmat(0.6,1,m.Cint(nm)); % this is the F exponent of the noise
end

% Leadfield
m.obs.LF = [0.10 0.3 0.10 0.5]; % contributions
m.obs.Cnoise = 0.2; % contributions
R.obs.Cnoise = m.obs.Cnoise; %legacy implementations

%% Prepare Priors
% 1 MMC
% 3 GPE
% 4 STN

% Excitatory connections
p.A{1} =  repmat(-32,m.m,m.m);
p.A_s{1} = repmat(0,m.m,m.m);

p.A{2} =  repmat(-32,m.m,m.m);
p.A_s{2} = repmat(0,m.m,m.m);

% Modulations
p.B{1} =  repmat(-32,m.m,m.m);
p.B_s{1} = repmat(0,m.m,m.m);
p.B{2} =  repmat(-32,m.m,m.m);
p.B_s{2} = repmat(0,m.m,m.m);


% Leadfield
p.obs.LF = [0];
p.obs.LF_s = repmat(1,size(p.obs.LF));

p.obs.Cnoise = [0];
p.obs.Cnoise_s = repmat(1/2,size(p.obs.Cnoise));

p.obs.mixing = [1]; %zeros(size(R.obs.mixing));
p.obs.mixing_s = repmat(0,size(p.obs.mixing));

p.obs.AlpNoise = 0;
p.obs.AlpNoise_s = repmat(1/8,size(p.obs.AlpNoise));

% Delays
p.DExt = repmat(-32,size(p.A{1})).*~((p.A{1}>-32) | (p.A{2}>-32)) ;
p.DExt_s = repmat(1/8,size(p.DExt));

% Sigmoid transfer for connections
% p.S = [0 0];
% p.S_s = [1/8 1/8];

% time constants and gains
for i = 1:m.m
    prec = 1/8;
    p.int{i}.T = zeros(1,m.Tint(i));
    p.int{i}.T_s = repmat(1/4,size(p.int{i}.T));
    p.int{i}.G = zeros(1,m.Gint(i));
    p.int{i}.G_s = repmat(1/4,size(p.int{i}.G));
    p.int{i}.DInt = zeros(1,m.Dint(i));
    p.int{i}.DInt_s = repmat(1/8,size(p.int{i}.DInt));

    p.int{i}.BG = zeros(1,m.Gint(i));
    p.int{i}.BG_s = repmat(prec,size(p.int{i}.BG));
    
    % Input strengths
    p.int{i}.C = zeros(1,m.Cint(i));
    p.int{i}.C_s = repmat(1/8,size(p.int{i}.C));
    p.int{i}.alpha = zeros(1,m.Cint(i));
    p.int{i}.alpha_s = repmat(1/8,size(p.int{i}.alpha));
    
    p.int{i}.S = zeros(1,m.Sint(i));
    p.int{i}.S_s = repmat(1/8,size(p.int{i}.S));    
end