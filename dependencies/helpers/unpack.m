function x= unpack(y,pts)
% unpacks data Y(1,:) into set of N epochs X(:,N) with start and indices in pts(N,1) and
% pts(N,2)
for i = 1:size(pts,1)
    if pts(i,2)<=size(y,1)
        x(:,i) = y(pts(i,1):pts(i,2));
    end
end

