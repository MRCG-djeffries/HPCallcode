function [qaly]=calc_qaly(X,tstart,tstop,tinterval,q_weight,discount_rate)
% X is the vector of person years same length of tstart:tinterval:tstop
% discount rate is the discounting for future qalys
% q_weight is the utility rate for the stage
% tstart is the first year
% tstop is the last year

% testing qlay1, qaly2 and qaly3 should be identical
% rng(227083)
% X=1000*rand(16,1);
% newX=interp1(2015:2030,X,2015:1/12:2030);
% [qaly1]=calc_qaly(X,66,81,1,0.5,3);
% [qaly2]=calc_qaly(X,2015,2030,1,0.5,3);
% [qaly3]=calc_qaly(newX,2015,2030,1/12,0.5,3);
t=tstart:tinterval:tstop;
if discount_rate>=1
    discount_rate = discount_rate/100;
end
r = 1 - discount_rate;
qaly = 0;
k = 0; % year discount counter
for i = tstart:(tstop-1)
    [~,~,j]=intersect((i:tinterval:(i+1)),t);
    x=t(j);
    y=X(j);
    I = trapz(x,y);
    qaly=qaly+I*1/((1+discount_rate)^k);%r^k;
    k = k + 1;
end
qaly=q_weight*qaly;