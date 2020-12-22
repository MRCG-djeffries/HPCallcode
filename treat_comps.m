function [phi,phidash]=treat_comps(X,scale00,scale10,n,t,pcom,alpha,time2015_inODEunits)
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

F3_former = 10:20:180;
F4_former = 11:20:180;
LT1_former = 14:20:180;
LT2_former = 15:20:180;

F3_pwider = 370:20:540;
F4_pwider = 371:20:540;
LT1_pwider =374:20:540;
LT2_pwider =375:20:540;

F0_former = 7:20:180;
F1_former = 8:20:180;
F2_former = 9:20:180;
F0_pwider = 367:20:540;
F1_pwider = 368:20:540;
F2_pwider = 369:20:540;


r_10=[F3_pwider,F4_pwider,LT1_pwider,LT2_pwider]; % advanced stages current

r_00=[F3_former,F4_former,LT1_former,LT2_former]; % advanced stages former

Adformer = scale00*sum(X(r_00)); % advanced stage totals
Adpwider = scale10*sum(X(r_10));

prop_treated = n/(Adformer+Adpwider);

phi=zeros(720,1);


j=[F3_former F4_former LT1_former LT2_former F3_pwider F4_pwider LT1_pwider LT2_pwider ];


if t>=time2015_inODEunits
    % how any treated in F3,F4,LT1,LT2
    phi(j)=min(prop_treated*X(j),X(j));
    % if there are more treats than subjects then there will be surplus
    ndash=n-scale_sumphi(phi,scale00,scale10,[F3_former F4_former LT1_former LT2_former],...
                         [F3_pwider F4_pwider LT1_pwider LT2_pwider]);
    if ndash > 1e-5
        % n was bigger than the denominatr of prop_treated
        % add extra allocation to phi overall all chronic for treat
        % i.e F0,F1,F2,F3,F4,LT1,LT2
        jformer=[F0_former F1_former F2_former F3_former F4_former LT1_former LT2_former];
        jpwider=[F0_pwider F1_pwider F2_pwider F3_pwider F4_pwider LT1_pwider LT2_pwider];
        j=[jformer jpwider];
        % needs to be scaled as ndash is scaled
        scaledbot=scale00*sum(X(jformer))+scale10*sum(X(jpwider));
        phi(j)=phi(j)+ndash*X(j)/scaledbot; % extra allocation
    end
end


% now move the treated to the T compartments
phidash=zeros(720,1);
T0_former = 16:20:180;
T1_former = 17:20:180;
T2_former = 18:20:180;
T3_former = 19:20:180;
T4_former = 20:20:180;
T0_pwider = 376:20:540;
T1_pwider = 377:20:540;
T2_pwider = 378:20:540;
T3_pwider = 379:20:540;
T4_pwider = 380:20:540;

phidash(T0_former)=alpha*phi(F0_former);
phidash(T1_former)=alpha*phi(F1_former);
phidash(T2_former)=alpha*phi(F2_former);
phidash(T3_former)=alpha*phi(F3_former);
phidash(T4_former)=alpha*(phi(F4_former)+phi(LT1_former)+phi(LT2_former));

phidash(T0_pwider)=alpha*pcom*phi(F0_pwider);
phidash(T1_pwider)=alpha*pcom*phi(F1_pwider);
phidash(T2_pwider)=alpha*pcom*phi(F2_pwider);
phidash(T3_pwider)=alpha*pcom*phi(F3_pwider);
phidash(T4_pwider)=alpha*pcom*(phi(F4_pwider)+phi(LT1_pwider)+phi(LT2_pwider));


function sumy=scale_sumphi(phi,scale00,scale10,formerj,pwiderj)
% phi conatins the number of subjects treated
% it needs tobe scaled to compare to total treatments
% This is for p=1
% formerj and pwiderj are the index of
% F3, F4, LT1 and LT2 compartments
sumy = scale00*sum(phi(formerj))+scale10*sum(phi(pwiderj));





