%% plot_iri_profs

fname = 'auroral_e_dens.txt';
%fname = 'daytime_dens.txt';

vals = str2num(asciiread(fname));

dene = vals(:, 2);
dene(dene <= 0) = 1;
plasmafreq = sqrt(80.6 * dene ./ 1E12);
plot(plasmafreq, vals(:, 1), '-k', 'LineWidth', 4)
xlim([0, 14])
xlabel('IRI Plasma freq. (MHz)')
ylabel('Alt. (km)')
set(gca, 'FontSize', 20)
grid on
grid minor
