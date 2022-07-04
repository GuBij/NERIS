function params=updateConstDir(maxline,option)
 global bool

 bool=0;
% fdir='./constant/';

 if nargin == 0
   params=readParams();
 else
   if isstr(maxline)
	maxline=str2num(maxline);
   end

   if nargin>1 && strcmp(option,'stations')
	bool=1;
%   elseif nargin>1
%	fdir=option;
   end
   params=readParams(maxline);
 end

 if bool==0
   params=calcParams(params);
 end

 createFile(params);
 if nargout==0
   quit;
 end
end
