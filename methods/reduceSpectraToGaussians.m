function ymodel = reduceSpectraToGaussians(X,Y,peakN,peakrange,minProm,plotop)
[hgt,Xpk,W,p] = findpeaks(Y,'MinPeakProminence',minProm);
pkLoc = X(Xpk);

% Restrict range
Xpk = Xpk(pkLoc>peakrange(1) & pkLoc<peakrange(2)); % only accept peaks within range
hgt = hgt(pkLoc>peakrange(1) & pkLoc<peakrange(2));

%initialize model
if numel(Xpk)>0
    ymodel = zeros(size(Y));
    % sort by height
    [hgt,ia] = sort(hgt,'descend');
    Xpk = Xpk(ia);
    
    % find difference in peaks
    prcdiff = abs(100*(diff(hgt)./hgt(1)));
    Xpk(find(prcdiff<15)+1) = [];
    hgt(find(prcdiff<15)+1) = [];
    
    % restrict to 3 peaks
    if numel(Xpk)>3
        Xpk = Xpk(1:3);
    end
    
    if plotop; plot(X,Y); hold on; end
    
    % now iteratively fit
    for i = 1:numel(Xpk) % first 3 peaks
        % find range before next inflection
        infPnts = find(abs(diff(sign(diff(Y))))>0)+ 1;
        peakPnt = find((Xpk(i)-infPnts)==0);
        try
            rangePk = infPnts(peakPnt-1):infPnts(peakPnt+1);
        catch
            warning('Not enough inflection points found,using max width')
            maxwid = 5/diff(X(1:2)); % max half width
            % establish which end is missing
            if (peakPnt-1)<=0
                rangePk = (infPnts(peakPnt)-maxwid):infPnts(peakPnt+1);
            elseif(peakPnt+1)>numel(peakPnt)
                rangePk = infPnts(peakPnt-1):(infPnts(peakPnt)+maxwid);
            else
                rangePk = (infPnts(peakPnt)-maxwid):(infPnts(peakPnt)+maxwid);
            end
            
            % truncate to make sure in range
            rangePk(rangePk>numel(X)) = [];
            rangePk = fix(rangePk);
        end
        gammaEst = (range(rangePk).*(X(2)-X(1)))/2;
        % fit cauchy over that range
        modelFun =  @(p,x)  p(1)./(1+((x-p(2))/p(3)).^2);
        startingVals = [Y(Xpk(i)) X(Xpk(i)) gammaEst]; % hgt centre gamma
        nlModel = fitnlm(X(rangePk),Y(rangePk),modelFun,startingVals);
        
        peakPrt(:,i) = predict(nlModel,X');
        if i>1
            peakPrt(peakPrt(:,i)<ymodel',i) = 0;
            peakPrt(peakPrt<peakPrt(:,i)) = 0;
        end
        ymodel = sum(peakPrt,2)';
        
        if plotop; plot(X,  peakPrt(:,i)); end
    end
    gwid = 5/diff(X(1:2));
    ymodel =  smoothdata(ymodel,'gaussian',gwid);
    ymodel = (ymodel-mean(ymodel))./std(ymodel);
    ymodel = ymodel-min(ymodel);
    
else % if no peaks then flat
    gwid = 5/diff(X(1:2));
    ymodel =  smoothdata(Y,'gaussian',gwid);
    ymodel = (ymodel-mean(ymodel))./std(ymodel);
    ymodel = ymodel-min(ymodel);
end

if plotop; plot(X,ymodel); end
