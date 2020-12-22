function py=calc_ongoingF4costs(outyt,outyF4z,year_start,year_end,mortality_rates,discount_rate)
% year_start and year_end are in the scale of 0:81
% so for 2015 to 2030, 66 to 81
%load('testy')
if discount_rate>=1
    discount_rate = discount_rate/100;
end
t=outyt;
Y=outyF4z;
Iout=zeros(9,year_end-year_start);
for i=1:9 % ages
    z=Y(i,:);
    I=integrate_over_year(t,z,year_start,year_end);
    Iout(i,:)=I;
end
py=0;
for i= 1:9
    for j = 1:(year_end-year_start)
        mortality_probs=1-exp(-mortality_rates(i)*(j:(year_end-year_start-j+1)));
        if length(mortality_probs)>1
            d=[mortality_probs(1) diff(mortality_probs)];
        else
            d=mortality_probs;
        end
        k=0:(length(mortality_probs)-1);
        py=py+sum(Iout(i,j)*(1-d).*(1./((1+discount_rate).^k)));
    end
end




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