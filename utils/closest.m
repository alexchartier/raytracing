function [X, Diff, idx] = closest(A,B)
%% CLOSEST.M
% % Find the closest element to B in A
% [X, Diff, idx] = closest(A,B)

%% Find the closest element

idx = find(abs(A-B) == min( abs(A - B) ))
X = A(idx(1));
Diff = min(abs( A - B ));
idx = A == X;