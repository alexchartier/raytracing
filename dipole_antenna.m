%% simulate offset dipole antenna performance

d=dipole ('length', 12, 'width', 0.02); %, 'FeedOffset', 2);
d_off=dipole ('length', 12, 'width', 0.02, 'FeedOffset', 2);

%% run analysis
clf
freq = linspace(1E6, 20E6, 19);
Z = 200;
%impedance(d, freq)

s = sparameters(d,freq, Z);
s_off = sparameters(d_off,freq, Z);
hold on
rfplot(s, 'k')
rfplot(s_off, 'r')
legend({'12-m dipole', '4/8-m offset dipole'})
hold off

set(gca, 'LineWidth', 2)

%patternAzimuth(d, 2E6)



%%
subplot(1, 2, 1)
subplot(1, 2, 2)
hg = smithplot(s,1,1,'GridType','ZY');
hg.LineStyle = '--';
hg = smithplot(s_off,1,1,'GridType','ZY');
hg.LineStyle = '--';



%% Ground plane whip
    % Example 1:
    % Create a reflector backed dipole with length of 30cm, width of 25cm
    % and spaced 7.5cm from the dipole for operation at 1 GHz.
clf
m = monopole('GroundPlaneLength',10,'GroundPlaneWidth',10);

%show(m)
%pattern(m, 5E6)
freq = linspace(5E6, 30E6, 6);
Z = 50:100:10000;

s = sparameters(m,freq, Z);

hold on
rfplot(s, 'k')

%legend({'12-m dipole', '4/8-m offset dipole'})
hold off

set(gca, 'LineWidth', 2)



 