function [R,m] = getStateDetailsJR(R)
m.m = numel(R.nmsim_name); % # of sources
m.un = 0; %initialize
for nm = 1:m.m
        m.dipfit.model(nm).source = R.nmsim_name{nm};

switch R.nmsim_name{nm}
    case 'THAL'
        m.x{nm} = [0 0];
        m.Gint(nm) = 1;
        m.Tint(nm) = 1;
        m.Sint(nm) = 2;
        m.Dint(nm) = 1;
        m.outstates{nm} = [1 0]; % relay cell voltage
        m.Cint(nm) = 1;
    case 'MMC'
        m.x{nm} = [0 0 0 0 0 0 0 0];
        m.Gint(nm) = 14;
        m.Tint(nm) = 4;
        m.Sint(nm) = 9;
        m.Dint(nm) = 14; % 14 delays
        m.Cint(nm) = 8; % number of noisy inputs
        m.outstates{nm} = [1 0 1 0 1 0 1 0]; % Superficial layer voltage
        m.UInt{nm} = [0 1 0 1 0 1 0 1]; % states recieving input
       case 'STN'
        m.x{nm} = [0 0];
        m.Gint(nm) = 0;
        m.Tint(nm) = 1;
        m.Sint(nm) = 1;        
        m.outstates{nm} = [1 0]; %
       case 'GPE'
        m.x{nm} = [0 0];
        m.Gint(nm) = 0;
        m.Tint(nm) = 1;
        m.Sint(nm) = 1;        
        m.outstates{nm} = [1 0]; % 
       case 'STR'
        m.x{nm} = [0 0];
        m.Gint(nm) = 1;
        m.Tint(nm) = 1;
        m.Sint(nm) = 2;        
        m.outstates{nm} = [1 0]; % 
       case 'GPI'
        m.x{nm} = [0 0];
        m.Gint(nm) = 1;
        m.Tint(nm) = 1;
        m.Sint(nm) = 2;        
        m.outstates{nm} = [1 0]; % 
end
m.un = m.un + m.Cint(nm); % this is total number of noisy inputs
end

R.obs.outstates = find([m.outstates{:}]);
for i=1:numel(R.chdat_name)
    R.obs.obsstates(i) = find(strcmp(R.chdat_name{i},R.chsim_name));
end
m.n =  size([m.x{:}],2); % Number of states
% Precompute xinds to make things easier with indexing
% Compute X inds (removes need to spm_unvec which is slow)
xinds = zeros(size(m.x,2),2);
for i = 1:size(m.x,2)
    if i == 1
        xinds(i,1) = 1;
        xinds(i,2) = size(m.x{i},2);
    else
        xinds(i,1) = xinds(i-1,2)+1;
        xinds(i,2) = xinds(i,1) + (size(m.x{i},2)-1);
    end
end
m.xinds = xinds;