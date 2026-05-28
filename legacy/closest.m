function [X, Diff] = closest(A,B)
%% CLOSEST.M
% Find the closest element to B in A

%% Find the closest element
X = A( abs(A-B) == min( abs(A - B) ) );
Diff = min(abs( A - B ));