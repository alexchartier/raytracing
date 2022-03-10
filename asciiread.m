function D = asciiread( File, Size )
% Function to read ASCII data into a character matrix
%    D = asciiread( File, Size )
% Input Arguments:
%    File = File name
% Optional Input Arguments:
%    Size = Bytes to read (defaults to inf) or [Rows,Cols]
% Output Arguments:
%    D     = ASCII file stored as character matrix
%
% Example. Read first four columns of ASCII file
%    Sites = asciiread('./data/sitelists/IGS_NORTH.sites',[Inf,4]);
% ----------------------------------------------------------------------------
D = [];

if nargin == 1; Size = inf; end

try
   fp = fopen(File,'r');
   D  = fread(fp,prod(Size),'uint8=>uint8');
   fclose(fp);

   if isempty(D); return; end
   D  = char(strread(char(D),'%s','delimiter','\n','whitespace',''));

   if length(Size)==2
      Size = min([size(D);Size]);
      D = D(1:Size(1),1:Size(2));
   end
catch
   fprintf('Warning from asciiread.m: ');
   if ~exist(File,'file'); fprintf('%s not found\n',File); return; end
   Err = lasterror; fprintf('%s\n',Err.message);
end
