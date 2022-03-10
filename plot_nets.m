%%
clf
N = 1:100;
m = 10;
hold on
plot(N, N, 'xr', 'DisplayName', 'a) Unconnected Group', 'MarkerSize', m)
plot(N, (N./2).^2 , 'og', 'DisplayName', 'b) Heterogeneous Network', 'MarkerSize', m)
plot(N, N + N.*(N -1)./2, '*m', 'DisplayName', 'c) Homogeneous Network', 'MarkerSize', m)
hold off
legend
set(gca, 'FontSize', 20)
xlabel('Number of nodes')
ylabel('Max Number of Soundings')
grid on
grid minor