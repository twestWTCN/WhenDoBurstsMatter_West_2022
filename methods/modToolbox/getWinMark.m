function winmark = getWinMark(X)
mark = zeros(size(X));
mark([1 end],:) = 1;
mark = mark(:);
winmark{1} = mark; % this marks the edge of the windows
