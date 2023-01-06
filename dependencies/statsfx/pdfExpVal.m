function [EV,SIGMA] = pdfExpVal(pdf,sup)
% Equations yielding the expected value and variance of arbitrary
% probability density function
%https://ltcconline.net/greenl/courses/117/DoubIntProb/ExpValVariance.htm

% Expected value
EV =  sum(pdf.*sup*diff(sup(2:3)));
% Standard deviation
SIGMA = sqrt(sum((pdf.*diff(sup(2:3))).*((sup-EV).^2)));
