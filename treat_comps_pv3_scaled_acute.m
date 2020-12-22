function [phi,phidash,moved_to_T4_former,moved_to_T4_current,total_not_T4,total_T4]=treat_comps_pv3_scaled_acute(treat_per,y1,X,age_scale00,age_scale10,t,pcom,alpha,time2015_inODEunits,pwid_weights,former_weights)
% This allows simultaneous applications of p=0 and p=1
% p=0 is PWID 
% p=1 is adanced for PWID and former
% y0 and y1 are the number of years p=0 and p = 1 are applied for
% y0 not used so dropped from argument list
% X are the compartments
% structure
% 180 (00), 180 (01) ,180 (10),180 (11)
% each 9*20
% scale00 and scale10 ar ethe former and current sacling factors
% n is the number of treatments per year
% t is the time of the slice
% NOTE it is assumed to be actual year, i.e. not 0:1/12:80
% but 1950:1/12:2030
% pcom is the probability of pwid completing treatement
% alpha is the SVR probability
% time2015_inODEunits is the year 2015 in the ODE timesacle
% so for example 0:1/12:80 for 1950:2030
% then time2015_inODEunits=781

% also outputs the number of subjects going to T4
if t>=y1
    treaty = treat_per;
else
    treaty = 0;
end

phi=zeros(756,1);
A_pwid = 384:21:567;

phi(A_pwid)=treaty*X(A_pwid);
% now move the treated to the T compartments

phidash=zeros(756,1);
A_pwid_treat = 399:21:567;
phidash(A_pwid_treat)=alpha*pcom*phi(A_pwid);
moved_to_T4_former=0;
moved_to_T4_current=0;
total_not_T4=0;
total_T4=0;
% test
%moved_to_T4=scale00*alpha*(phi(F3_former)+phi(F4_former)+phi(LT1_former)+phi(LT2_former))+...
%    scale10*alpha*pcom*(phi(F3_pwider)+phi(F4_pwider)+phi(LT1_pwider)+phi(LT2_pwider));
% total_not_T4 = age_scale10*pwid_weights'.*(phi(F0_pwider)+phi(F1_pwider)+phi(F2_pwider)+phi(F3_pwider)+...
%                phi(LT1_pwider)+phi(LT2_pwider))+...
%                age_scale00*former_weights'.*(phi(F0_former)+phi(F1_former)+phi(F2_former)+phi(F3_former)+...
%                phi(LT1_former)+phi(LT2_former));
% total_T4=age_scale00*former_weights'.*phi(F4_former)+ age_scale10*pwid_weights'.*phi(F4_pwider);       
% 
% moved_to_T4_former=age_scale00*alpha*former_weights'.*phi(F4_former);
% moved_to_T4_current=age_scale10*alpha*pcom*pwid_weights'.*phi(F4_pwider); % 9 by 1 for ages
function sumy=scale_sumphi_p_scaled(phi,age_scale00,age_scale10,jAdformer,jAdpwider,jMipwider,pwid_weights,former_weights)
% phi conatins the number of subjects treated
% it needs tobe scaled to compare to total treatments
% This is for any p
% jAdformer is the former advanced compartments
% jAdpwider is the current advanced compartments
% jMipwider is the current mild compartments
% sumy below is the f comps + h comps + g compartments
Apwid_A=[pwid_weights pwid_weights pwid_weights pwid_weights]';
Apwid_M=[pwid_weights pwid_weights pwid_weights]';
Aformer=[former_weights former_weights former_weights former_weights]';
sumy = age_scale10*sum(Apwid_A.*phi(jAdpwider))+...
       age_scale10*sum(Apwid_M.*phi(jMipwider))+...
       age_scale00*sum(Aformer.*phi(jAdformer));





