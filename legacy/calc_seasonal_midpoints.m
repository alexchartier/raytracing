
%%

sols_eq = datenum(2021, 3:3:15, 21);

for t = 1:length(sols_eq) - 1
    midpt = (sols_eq(t+1) - sols_eq(t)) / 2 + sols_eq(t);
    disp(datestr(midpt))
end
