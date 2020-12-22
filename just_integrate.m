function I=just_integrate(X,tstart,tstop,tinterval)
% 5 September 2020
% note becuase of rounding, below is not consistnet
%  [~,~,j]=intersect((i:tinterval:(i+1)),t);
%  use below instead
%  j=find(t==i):find(t==(i+1));
% Integrates the compartments into yearly totals
t=tstart:tinterval:tstop;
I=zeros(1,tstop-tstart);
for i = tstart:(tstop-1)
    %j=find(t==i):find(t==(i+1));
    j=find(t>=i & t<=(i+1));
    x=t(j);
    y=X(j);
    I(i-tstart+1) = trapz(x,y);
end
