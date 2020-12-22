function [phi,phidash]=treat_comps_pv2(p,X,scale00,scale10,n,t,pcom,alpha,time2015_inODEunits)
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

c_10=[F0_pwider F1_pwider F2_pwider F3_pwider,F4_pwider,LT1_pwider,LT2_pwider]; % all chronic stages current

Topwider=scale10*sum(X(c_10));
Adformer = scale00*sum(X(r_00)); % advanced stage totals
Adpwider = scale10*sum(X(r_10));
Adtotal=Adformer+Adpwider;
omega = Adpwider/(p*Adtotal+(1-p)*Topwider);
% This automates on p
fprop =  n*omega/Adpwider;
gprop = n*p*(1-omega)/Adformer;
hprop = n*(1-p)*(1-omega)/(Topwider-Adpwider);

phi=zeros(720,1);

jAdpwider=[F3_pwider F4_pwider LT1_pwider LT2_pwider ]; % Advanced pwid
jMipwider=[F0_pwider F1_pwider F2_pwider]; % Mild pwid
jAdformer=[F3_former F4_former LT1_former LT2_former ]; % Advanced former


if t>=time2015_inODEunits
    % how any treated in F3,F4,LT1,LT2
    % f function
    phi(jAdpwider)=min(fprop*X(jAdpwider),X(jAdpwider));
    phi(jAdformer)=min(gprop*X(jAdformer),X(jAdformer));
    phi(jMipwider)=min(hprop*X(jMipwider),X(jMipwider));
    
    
    
    % if there are more treats than subjects then there will be surplus
    ndash=n-scale_sumphi_p(phi,scale00,scale10,jAdformer,jAdpwider,jMipwider);
                         
    if ndash > 1e-5
        if p==1
            %Treatments were allocated proportionally across injecting drug use status 
            %(i.e. patients were treated irrespective of whether they were current of former PWID)
            %, and if the number of treatments available was greater than the number of 
            %HCV-infected people with advanced liver disease, 
            %the remaining treatments were allocated to HCV-infected people with early liver disease, 
            %proportionally across injecting drug use status

            % n was bigger than the denominatr of prop_treated
            % add extra allocation to phi overall all chronic for treat
            % i.e F0,F1,F2,F3,F4,LT1,LT2
            jformer=[F0_former F1_former F2_former ];
            jpwider=[F0_pwider F1_pwider F2_pwider ];
            j=[jformer jpwider];
            % needs to be scaled as ndash is scaled
            scaledbot=scale00*sum(X(jformer))+scale10*sum(X(jpwider));
            phi(j)=phi(j)+ndash*X(j)/scaledbot; % extra allocation
        else % p=0
            %Treatments were allocated proportionally across liver disease stages, 
            %and if the number of treatments available was greater than the number 
            %of current PWID eligible for treatment, 
            %remaining treatments were allocated to former PWID, 
            %proportionally across disease stages.
            jformer=[F0_former F1_former F2_former F3_former F4_former LT1_former LT2_former];
            j=jformer;
            scaledbot=scale00*sum(X(jformer));
            phi(j)=phi(j)+ndash*X(j)/scaledbot; % extra allocation
        end
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


function sumy=scale_sumphi_p(phi,scale00,scale10,jAdformer,jAdpwider,jMipwider)
% phi conatins the number of subjects treated
% it needs tobe scaled to compare to total treatments
% This is for any p
% jAdformer is the former advanced compartments
% jAdpwider is the current advanced compartments
% jMipwider is the current mild compartments
% sumy below is the f comps + h comps + g compartments
sumy = scale10*sum(phi(jAdpwider))+...
       scale10*sum(phi(jMipwider))+...
       scale00*sum(phi(jAdformer));





