function D2 = interp_sami(D, alts, time)
% 
% D = load_sami('test.nc');
% alts = 90:2:400;
% time = datetime(2014, 5, 23, 6, 0, 0);
%  
% % plot
% ai = D.alt <= max(alts);
% subplot(2, 1, 1)
% plot(D.alt(ai), D.dene0(ai, 10, 10, tind))
% 
% subplot(2, 1, 2)
% plot(D2.alt, D2.dene0(:, 10, 10))

%%
D2 = D;
size_D = size(D.dene0);
size_D(1) = length(alts);
D2.dene = nan(size_D(1:end - 1));
D2.alt = alts;
tind = D.time == time;

dene_in = D.dene0(:, :, :, tind);

for i1 = 1:length(D.lon)
    for i2 = 1:length(D.lat)
        D2.dene(:, i1, i2) = interp1(D.alt, dene_in(:, i1, i2), alts);
    end
end


