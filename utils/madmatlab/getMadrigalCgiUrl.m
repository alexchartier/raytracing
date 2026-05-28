function cgiUrl = getMadrigalCgiUrl(url)
%  getMadrigalCgiUrl  	parse the main madrigal page to get the cgi url
% 
%  With Madrigal 3, this method simply returns the original url.
%
%  input: url to Madrigal
%
%  output: cgi url for that Madrigal Site
%
%  Note: parses the homepage for the accessData link

% get main page
  
if url(end) ~= '/'
    result = findstr(url, 'index.html');
    if length(result) == 0
        url = strcat(url,'/');
    end
end
these_options = weboptions('Timeout',300, 'ContentType', 'text');
pagedata = webread(url, these_options);

% get host name
[proto,ppath]=strtok(url,':');
ppath=ppath(4:end);
[host,page]=strtok(ppath,'/');
[host,port]=strtok(host,':');
if length(port) > 1
    port = port(2:end);
end
    
% get cgi-bin path:
% This searches for strings like
% <!-- This html comment exists simply to support the old remote API's: "/accessData.cgi -->
% or <A HREF="/madrigal/cgi-bin/accessData.cgi">
index1 = regexp(pagedata, '[^"]*accessData.cgi', 'once');
% check for error
if length(index1) == 0
    err.message = 'No Madrigal home page found at given url';
    err.identifier = 'madmatlab:badArguments';
    rethrow(err);
end
index2 = regexp(pagedata, 'accessData.cgi', 'once');
if index2 - index1 > 1
    % longer than just "/"
    cgibin=pagedata(index1:index2-1);
else
    cgibin=page;
end
    
% Build URL
if port
    cgiUrl = strcat(proto,'://',host,':',port,cgibin);
else
    cgiUrl = strcat(proto,'://',host,cgibin);
end

