function wspr_out = declump_wspr(wspr_t, npts_wanted)
%% % Example inputs
% times = datenum(2019, 3, 15, 0:1, 0, 0);
% wspr_fn = 'data/wspr/wsprspots-2019-03.csv';
% fmin = 2;
% fmax = 30;
% snrmin = -25;
% distmin = 500;
% distmax = 2000;
% npts_wanted = 100;
% 
% % Load WSPR data as example
% wspr = load_wspr_csv(wspr_fn, min(times), max(times), fmin, fmax, ...
%     snrmin, distmin, distmax);
% 
% 
% %% Select WSPR at just one time
% time = times(1);
% wspr_ti = abs(wspr.times - datenum(time)) < 5/60/24;
% fn = fieldnames(wspr);
% wspr_t = [];
% for k = 1:numel(fn)
%     wspr_t.(fn{k}) = wspr.(fn{k})(wspr_ti);
% end
% 
% 
% %%
% 
% 
% hold on
% scatter(wspr_t.lonI, wspr_t.latI, 'b')
% scatter(wspr_out.lonI, wspr_out.latI, '.k')
% xlim([-180, 180])
% ylim([-90, 90])
% hold off


%% Calculate the great-circle midpoints
npts = length(wspr_t.times);
wspr_t.latI = nan([npts, 1]);
wspr_t.lonI = nan([npts, 1]);

for l = 1:npts
    [wspr_t.latI(l), wspr_t.lonI(l)] = meanm([wspr_t.txlat(l), wspr_t.rxlat(l)], ...
        [wspr_t.txlon(l), wspr_t.rxlon(l)]);
end

% Kick out the non-unique elements
[~, IA, ~] = unique(wspr_t.latI + wspr_t.lonI);
fn = fieldnames(wspr_t);
for k = 1:numel(fn)
    wspr_t.(fn{k}) = wspr_t.(fn{k})(IA);
end

if length(IA) <= npts_wanted
    wspr_out = wspr_t;
else
    % declump
    [latI, lonI, idx] = declump(wspr_t.latI', wspr_t.lonI', npts_wanted);
    wspr_out = [];
    fn = fieldnames(wspr_t);
    for k = 1:numel(fn)
        wspr_out.(fn{k}) = wspr_t.(fn{k})(idx);
    end
end
  

%%
%    clf
%
%    hold on
%    m_proj('Miller Cylindrical', 'latitude', [-80 80], 'longitude', [-180, 180])
%
%    m_scatter(wspr_t.lonI, wspr_t.latI, 'b')
%    m_scatter(wspr_out.lonI, wspr_out.latI, '.m')
%
%
%        m_grid('color', 'k', 'FontSize', 20)
%        m_coast('color', 'k');
%    % xlim([-180, 180])
%    % ylim([-90, 90])
%    hold off
%
%    legend({"", "WSPR link midpoint", 'Selected by declumping algorithm'}, 'FontSize', 20)

















