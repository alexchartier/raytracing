function UT = get_pharlap_time(time)
%% Specify time in PHaRLaP format
year = str2double(datestr(time, 'yyyy'));
month = str2double(datestr(time, 'mm'));
day  = str2double(datestr(time, 'dd'));
hour = str2double(datestr(time, 'HH'));
minute = str2double(datestr(time, 'MM'));
UT = [year month day hour minute];
