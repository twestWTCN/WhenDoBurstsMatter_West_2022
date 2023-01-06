function [data,chloc,indM1,dist] = findClosest2ROI(data,loc)
        % taken from https://www.researchgate.net/publication/273471194_The_Frontal_Control_of_Stopping
        roi = [37 -25 62; -37 -25 62];
        dist = [];
        for lr = 1:2
            V = loc-roi(lr,:);
            dist(:,lr) = arrayfun(@(ROWIDX) norm(V(ROWIDX,:)), (1:size(V,1)).'); % Euclidean distance from ROI in mm
        end
        indM1 = find(sum(dist<30,2)); % find channels that are less than 50mm from ROI
        % find side
        side = find(sum(dist<30,1));
        data = data(:,indM1);
        chloc = loc(indM1,:);
        dist = dist(indM1,side);