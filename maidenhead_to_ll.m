function [lat, lon] = maidenhead_to_ll(code)

Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
kv = 0:25;


cr = cellstr(upper(code)')';


% Lon
Field = kv(Alphabet == cr{1}) * 20;
Square = str2double(cr{3}) * 2;
SubSquareLow = kv(Alphabet == cr{5}) * (2/24);
SubSquareHigh = SubSquareLow + (2/24);

StartLon = Field + Square + SubSquareLow - 180;
EndLon = Field + Square + SubSquareHigh - 180;

lon = mean([StartLon, EndLon]);


% def GetLat(TWO, FOUR, SIX):

Field = kv(Alphabet == cr{2}) * 10;
Square = str2double(cr{4});
SubSquareLow = kv(Alphabet == cr{6}) * (1/24);
SubSquareHigh = SubSquareLow + (1/24);

StartLat = Field + Square + SubSquareLow - 90;
EndLat = Field + Square + SubSquareHigh - 90;
lat = mean([StartLat, EndLat]);


% def main(strMaidenHead = __MH__):
