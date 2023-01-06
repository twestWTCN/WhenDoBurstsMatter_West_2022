function [st_pts,end_pts,isi_pts,flex,stidnSave] = gen_inv_pts_TW(R,subject)
% function modified from Kai Millers toolbox.

%% set parameters, load finger flexion
fsamp=1000;  %sampling frequency

pack, disp(subject)
load([R.path.KM_exp '/data/' subject '/' subject '_fingerflex'],'flex','cue'), clear data
dlength=size(flex,1);

% Smooth flex data
gwin = gausswin(0.3*fsamp);

for dig = 1:5
    flexS = conv(flex(:,dig),gwin,'same');
    X = find(cue==dig);
    XSeg = SplitVec(X,'consecutive'); % split up cue periods
    if numel(XSeg)<2
        warning('Cue trace is not available, try to reconstruct from movement')
        X = find(flexS>(0.65*max(flexS)));
        XSeg = SplitVec(X,'consecutive'); % split up cue periods
        xs = 0;
        for s = 1:numel(XSeg)
            X = XSeg{s};
            XSeg{s} = X-(0.5*fsamp);
        end
    end
    
    
    end_pts = []; st_pts = []; s = 0; ofs = 0;
    for seg = 1:numel(XSeg)
        if (XSeg{seg}(1)-0.5*fsamp)>1 && (XSeg{seg}(end)+2*fsamp)<numel(flexS)
            winInds = XSeg{seg}(1)-0.5*fsamp:XSeg{seg}(end)+1*fsamp;
            X = flexS(winInds);
            XD = diff(X); % find derivative
            % find inflection point
            stind = find(abs(XD)>100,1,'first'); % find when derivative is exceeded
            % correct for ofset ~100ms
            stind = stind-(0.1*fsamp); % account for rise time
            %% Find movement onset
            if stind>500
                s = s+1;
                stidnSave(s) = stind/fsamp; % reaction time
                st_pts(s) = stind+winInds(1);
            else
                a = 1;
                %less than cue (preemptive)
            end
%             plot(abs(XD)); hold on
%             if stind>1
%             scatter(stind,abs(XD(stind)))
%             end
%             yyaxis right
%             plot(cue(winInds))
%             clf
            %% Find movement offset
            winInds = XSeg{seg}(end)-0.5*fsamp:XSeg{seg}(end)+2*fsamp; % this is the cue end minus some time
            X = flexS(winInds);
            XD = diff(X);
            stind = find(abs(XD)>100,1,'last');
            % correct for ofset ~100ms
            stind = stind+(0.1*fsamp);
            if stind>500
                ofs = ofs+1;
                end_pts(ofs) = stind+winInds(1);
            else
                a = 1;
                %less than cue (preemptive)
            end
            %         plot(abs(XD)); hold on
            %         scatter(stind,abs(XD(stind)))
            %         yyaxis right
            %         plot(cue(winInds))
        end
    end
    %% Find Movement Ofset
    if dig == 1
        ch1_st_pts = st_pts; ch1_end_pts = end_pts;
    elseif dig == 2
        ch2_st_pts = st_pts; ch2_end_pts = end_pts;
    elseif dig == 3
        ch3_st_pts = st_pts; ch3_end_pts = end_pts;
    elseif dig == 4
        ch4_st_pts = st_pts; ch4_end_pts = end_pts;
    elseif dig == 5
        ch5_st_pts = st_pts; ch5_end_pts = end_pts;
    end
    
end

%% %%%%%generate random picks for "rest" condition
catEnds = [ch1_end_pts ch2_end_pts ch3_end_pts ch4_end_pts ch5_end_pts];
catStds = [ch1_st_pts  ch2_st_pts  ch3_st_pts  ch4_st_pts  ch5_st_pts];
isos = 0; isi_pts = [];
for seg = 1:numel(catEnds)
    x = catEnds(seg)+1*fsamp; % 1 second from end point
    xwin = x-1*fsamp:x+2*fsamp;
    xd = abs(x - catStds)<(1*fsamp); % check distance
    if ~any(xd)
        isos = isos+1;
        isi_pts(isos) = x;
    else
        a = 1;
        % too close
    end
    %     xc = sum(flex(xwin,:),2);
    %     plot(xc); hold on
    %     scatter(501,xc(501));
    %     clf
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% order points and keep label
%%%% pts=[movementstart movementinversion movementlabel]
st_pts = [...
    [ch1_st_pts' 1+0*ch1_st_pts'];...
    [ch2_st_pts' 2+0*ch2_st_pts'];...
    [ch3_st_pts' 3+0*ch3_st_pts'];...
    [ch4_st_pts' 4+0*ch4_st_pts'];...
    [ch5_st_pts' 5+0*ch5_st_pts'];...
    ];
end_pts = [...
    [ch1_end_pts' 1+0*ch1_end_pts'];...
    [ch2_end_pts' 2+0*ch2_end_pts'];...
    [ch3_end_pts' 3+0*ch3_end_pts'];...
    [ch4_end_pts' 4+0*ch4_end_pts'];...
    [ch5_end_pts' 5+0*ch5_end_pts'];...
    ];

isi_pts = [isi_pts' 1+0*isi_pts'];