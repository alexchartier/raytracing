function array = arange(start, inc, npts)
%% arange.m
% produce array of length 'npts' with spacing 'inc' and starting point
% 'start'
% array = arange(start, inc, npts)

array = start:inc:(start + inc * (npts - 1));