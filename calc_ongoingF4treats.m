function py=calc_ongoingF4treats(t,outyF4z,year_start,year_end)
% year_start and year_end are in the scale of 0:81
% so for 2015 to 2030, 66 to 81
%load('testy')
%[z,j]=sort(outyt);
%t=outyt(j);
%Y=outyF4z(:,j);
%Iout=zeros(9,year_end-year_start);
for i=1:9 % ages
    z=outyF4z(i,:);
    I=integrate_over_year(t,z,year_start,year_end);
    Iout(i,:)=I;
end
py=sum(Iout);



function I=integrate_over_year(t,X,tstart,tstop)
% 5 September 2020
% note becuase of rounding, below is not consistnet
%  [~,~,j]=intersect((i:tinterval:(i+1)),t);
%  use below instead
%  j=find(t==i):find(t==(i+1));
% Integrates the compartments into yearly totals
t=1950+t;
I=zeros(1,tstop-tstart);
for i = tstart:(tstop-1)
    j=find(t>=i & t<=(i+1));
    x=t(j);
    y=X(j);
    I(i-tstart+1) = trapz(x,y);
end