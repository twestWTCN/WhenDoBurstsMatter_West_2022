function [xstore_cond,tvec,wflag,J,Es] = ABC_fx_compile_NH(R,x,uc,pc,m)
% To Do:
% 1)Precompute the expectations of the within source parameters and take
%   outside of the integration loop.
% If you want to estimate the noise floor - then use 'decon' to deconnect
% both intrinsic/extrinsic couplings.
cs = 0; % cond counter
wflag= 0; tvec = [];
for condsel = 1:numel(R.condnames)
    cs = cs+1;
    us = uc{cs};
    p = pc;
    % Compiles NMM functions with delays, seperates intrinsic and extrinsic
    % dynamics then summates
    xinds = m.xinds;
    % model-specific parameters
    %==========================================================================
    % model or node-specific state equations of motions
    %--------------------------------------------------------------------------
    fx{1} = @ABC_fx_NH_mmc;                                    % motor micro circuit
    
    %% Get Intrinsic parameters
    [DInt,DMat] = recallIntDelayTableWS_emp(R,p,m);
    %% Precompute parameter expectations
    % Parameter Priors
%         pQ = getModelPriorsNH(m);
%         pQ = getModelPriorsNH_emp(m);
    pQ = getModelPriorsNH_emp(m);
    n     = m.m;
    
    nbank = cell(1,n); qbank = cell(1,n);
    for i = 1:n
        N.x  = m.x{i};
        nbank{i} = N;
        Q = p.int{i};
        Q.T = pQ(i).T.*exp(Q.T);
        Q.G = pQ(i).G.*exp(Q.G);
        Q.Mn = pQ(i).Mn.*exp(Q.Mn);
        Q.Sn = pQ(i).Sn.*exp(Q.Sn);
        Q.Bn = pQ(i).Bn.*exp(Q.Bn);
        Q.C  = pQ(i).C.*exp(Q.C);
        qbank{i} = Q;
    end
    
    %% TIME INTEGRATION STARTS HERE ===========================================
    f = zeros(xinds(end),1); dt = R.IntP.dt;
    DDt = DMat/dt; 
    % set buffers to match maximum delay
    R.IntP.bufferExt = 1+max(DDt(:));
    R.IntP.bufferInt =  1+max(DDt(:));
    if iscell(x)
        xstore= full(repmat(spm_vec(x),1,max([R.IntP.bufferExt R.IntP.bufferInt])));
    else
        xstore = x;
    end
    xint = zeros(m.n,1);
    TOL = exp(-4);
    for tstep = max([R.IntP.bufferExt R.IntP.bufferInt]):R.IntP.nt
        % intrinsic flow at target
        %----------------------------------------------------------------------
        i = 1;
        xi = xstore(m.xinds(i,1):m.xinds(i,2),tstep-R.IntP.bufferInt+1:tstep)';
        [f(m.xinds(i,1):m.xinds(i,2))] = fx{i}(xi,qbank{i},DInt{i});
        dW = us{i}(tstep,:).*sqrt(dt).*qbank{i}.C;
        xint = xint + (f*dt) + dW'; % EM Step
        
        xstore = [xstore xint]; % This is done for speed reasons! Faster than indexing (!!)
        
    end
    if wflag == 1
        xstore_cond{condsel} = NaN;
    end
    xstore_cond{condsel} = xstore;
    if nargout>3
        [J{condsel},Es{condsel}] = findJacobian_immediateVicinity(R,xstore(:,end-R.IntP.bufferExt:end),uc,p,m,condsel);

%         [Q,Js] = adapted_spm_dcm_delay(DMat,J{condsel})
%         eig(Q*Js)
        Es{condsel} = nan;
    end
    a = 1;
end
