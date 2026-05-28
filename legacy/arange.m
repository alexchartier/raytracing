function array = arange(start, inc, npts)
%% arange.m
% produce array of length 'npts' with spacing 'inc' and starting point
% 'start'
array = start:inc:(start + inc * npts);