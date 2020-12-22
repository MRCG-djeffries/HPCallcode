function fitode_v6_with_abc
% 30 August 2020
% 31 August 2020 chronic HCV prevlance calibrated at 50% for current PWID in 2015.
% 1 September this use the abc values run using abc_ode.R
%                  P1     P2     P3     P4     P5     P6
% Min.:        0.2207 0.0590 0.0940 0.0668 0.0643 0.0613
% 2.5% Perc.:  0.2223 0.0591 0.0948 0.0672 0.0643 0.0622
% Median:      0.2570 0.0632 0.1041 0.0778 0.0681 0.0785
% Mean:        0.2573 0.0653 0.1071 0.0790 0.0691 0.0775
% Mode:        0.2561 0.0614 0.1048 0.0771 0.0663 0.0800
% 97.5% Perc.: 0.2855 0.0791 0.1348 0.0982 0.0800 0.0915
% Max.:        0.2888 0.0794 0.1368 0.0997 0.0842 0.0957

% P1 is delta parm_current(1)=r1 = 0.2573, was 0.26
% P2 is parm_current(6) = -log(1- r2 ) = 0.0653, was 0.117  F0 to F1 curr
% P3 is parm_former(6) = -log(1-r3) = 0.1071, was 0.112     F0 to F1 form
% P4 is parm_current(7) = -log(1- r4 ) = 0.079, was 0.089  F1 to F2 curr
% P5 is parm_former(7) = -log(1- r5);  = 0.0691, was 0.077 F1 to F2 form
% P6 is piv r6=0.01+0.19*rand; % = 0.0775 was 1.58264*0.056=0.089 curr

% delete figures
%************************************
% **** note M matrices 15,15 should be minus - it is a LT2 death
% also set the svr4 to dc and hcc to zero, not correctly implemented

% Saturday Sept 12 added 1 year
delete(get(0,'children'))

%%
[compnames,chronic_nams,...
    chronic_nums00,chronic_nums01,chronic_nums10,chronic_nums11,...
    S0_comps,S1_comps,S2_comps,S3_comps,S4_comps,A_comps,...
    F0_comps,F1_comps,F2_comps,F3_comps,F4_comps,...
    DC_comps,HCC_comps,LT1_comps,LT2_comps,...
    T0_comps,T1_comps,T2_comps,T3_comps,T4_comps,...
    S0_comps_00,S0_comps_10,S0_comps_01,S0_comps_11,...
    S1_comps_00,S1_comps_10,S1_comps_01,S1_comps_11,...
    S2_comps_00,S2_comps_10,S2_comps_01,S2_comps_11,...
    S3_comps_00,S3_comps_10,S3_comps_01,S3_comps_11,...
    S4_comps_00,S4_comps_10,S4_comps_01,S4_comps_11,...
    A_comps_00,A_comps_10,A_comps_01,A_comps_11,...
    F0_comps_00,F0_comps_10,F0_comps_01,F0_comps_11,...
    F1_comps_00,F1_comps_10,F1_comps_01,F1_comps_11,...
    F2_comps_00,F2_comps_10,F2_comps_01,F2_comps_11,...
    F3_comps_00,F3_comps_10,F3_comps_01,F3_comps_11,...
    F4_comps_00,F4_comps_10,F4_comps_01,F4_comps_11,...
    LT1_comps_00,LT1_comps_10,LT1_comps_01,LT1_comps_11,...
    LT2_comps_00,LT2_comps_10,LT2_comps_01,LT2_comps_11,...
    DC_comps_00,DC_comps_10,DC_comps_01,DC_comps_11,...
    HCC_comps_00,HCC_comps_10,HCC_comps_01,HCC_comps_11,...
    T0_comps_00,T0_comps_10,T0_comps_01,T0_comps_11,...
    T1_comps_00,T1_comps_10,T1_comps_01,T1_comps_11,...
    T2_comps_00,T2_comps_10,T2_comps_01,T2_comps_11,...
    T3_comps_00,T3_comps_10,T3_comps_01,T3_comps_11,...
    T4_comps_00,T4_comps_10,T4_comps_01,T4_comps_11]=coeffsof_v5;
%%
% get the data
[N0,M,Mdash,mort_current,mort_former,extra_parms_vals, phiminusvals, phiplussvals,extra_parms_nams,age_matrix,parm_names,parm_current,parm_former]= create_params_v6atabc;
% solve the diff equation
[t1,x1] = ode15s(@odeeq_v6,0:81,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix);
X=x1';
% chronic HCV prevlance calibrated at 50% for current PWID in 2015.
chronic_prev_current = [chronic_nums10]; % with no treatments chronic_nums11 all zero
total_current = [361:540];
figure;plot(t1,100*sum(X(chronic_prev_current,:))./sum(X(total_current,:)),'r--','LineWidth',2)
set(gca,'ylim',[0 70])
xticks(0:5:80)
set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
title('% of current PWID with HCV from 1950 to 2030')
grid()

chronic_prev = [chronic_nums00 chronic_nums10]; 
total_all = 1:720;
figure;plot(t1,100*sum(X(chronic_prev,:))./sum(X(total_all,:)),'r--','LineWidth',2)
set(gca,'ylim',[0 70])
xticks(0:5:80)
set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
title('% of all PWID with HCV from 1950 to 2030')
grid()

scale00=144000/sum(X([chronic_nums00],66)); % former
scale10=40000/sum(X([chronic_nums10],66)); % current

figure;plot(t1,scale00*sum(X([chronic_nums00],:))+scale10*sum(X([chronic_nums10],:)))
xticks(0:5:80)
set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
title('Absolute numbers of PWID with chronic HCV, scaled to the 2015 estimates')
grid()
figure; plot(t1,scale10*sum(X([chronic_nums10],:))+scale10*sum(X(S0_comps(19:27),:)))
lt2=scale00*sum(X(LT2_comps_00,:))+scale10*sum(X(LT2_comps_10,:));
lt1=scale00*sum(X(LT1_comps_00,:))+scale10*sum(X(LT1_comps_10,:));
lt=lt1+lt2;
dc=scale00*sum(X(DC_comps_00,:))+scale10*sum(X(DC_comps_10,:));
hcc=scale00*sum(X(HCC_comps_00,:))+scale10*sum(X(HCC_comps_10,:));
f4=scale00*sum(X(F4_comps_00,:))+scale10*sum(X(F4_comps_10,:));
f3=scale00*sum(X(F3_comps_00,:))+scale10*sum(X(F3_comps_10,:));
f2=scale00*sum(X(F2_comps_00,:))+scale10*sum(X(F2_comps_10,:));
f1=scale00*sum(X(F1_comps_00,:))+scale10*sum(X(F1_comps_10,:));
f0=scale00*sum(X(F0_comps_00,:))+scale10*sum(X(F0_comps_10,:));


figure;area(t1,[lt' hcc' dc' f4' f3' f2' f1' f0'])
xticks(0:5:80)
set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
legend({'lt' 'hcc' 'dc' 'f4' 'f3' 'f2' 'f1' 'f0'},'location','best');
grid()
title('Number of PWID at each stage of chronic HCV')

f=sum(X(chronic_prev_current,:))./sum(X(total_current,:));
sum2=f(53); % target 0.6
sum3=sum(X(chronic_prev_current,66))/sum(X(total_current,66)); % target 0.5sum4= (f0(66)+f1(66))/(f0(66)+f1(66)+f2(66)+f3(66)+f4(66)+hcc(66)+dc(66)+lt(66))
sum4= (f0(66)+f1(66))/(f0(66)+f1(66)+f2(66)+f3(66)+f4(66)+hcc(66)+dc(66)+lt(66)); % target = 0.66

figure;plot(t1,sum(X,1))
set(gca,'ylim',[900 1100])
title('Should be constant - tests for leaks')


% % test the qaly calculation, resolution
% [t1,x1] = ode45(@odeeq_v4,0:80,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix);
% X=x1';
% Xtest=scale10*sum(X([chronic_nums10],:));
% [qaly1]=calc_qaly(Xtest(66:81),66,81,1,0.5,3);
%
% [t1,x1] = ode45(@odeeq_v4,0:1/2:80,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix);
% X=x1';
% Xtest=scale10*sum(X([chronic_nums10],:));
% [qaly2]=calc_qaly(Xtest(133:161),66,81,1,0.5,3);
%
% [t1,x1] = ode45(@odeeq_v4,0:1/10:80,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix);
% X=x1';
% Xtest=scale10*sum(X([chronic_nums10],:));
% [qaly4]=calc_qaly(Xtest(661:801),66,81,1,0.5,3);
%
% [t1,x1] = ode45(@odeeq_v4,0:1/12:80,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix);
% X=x1';
% Xtest=scale10*sum(X(chronic_nums10,:));
% [qaly5]=calc_qaly(Xtest(793:961),66,81,1,0.5,3);


% MUST DO seperately as the utility weights are all different
% chronic_nums10 is amde up of
compvec10=zeros(11,9); % 9 compartments * 9 age groups
% susceptible current PWID, then acute
% then 9 chronic compartments
compvec10= [S0_comps_10;A_comps_10;F0_comps_10;F1_comps_10;F2_comps_10;F3_comps_10;F4_comps_10;HCC_comps_10;DC_comps_10;LT1_comps_10;LT2_comps_10];
% same for 00, i.e. former PWID
% *** note S00 does not effect cost effectiveness ratio
compvec00=zeros(11,9); % 9 compartments * 9 age groups
% susceptible current PWID, then acute
% then 9 chronic compartments
compvec00= [S0_comps_00;A_comps_00;F0_comps_00;F1_comps_00;F2_comps_00;F3_comps_00;F4_comps_00;HCC_comps_00;DC_comps_00;LT1_comps_00;LT2_comps_00];

utility_vec=[0.93 0.77 0.77 0.77 0.77 0.66 0.55 0.45 0.45 0.45 0.67]; % These are the utiltiy multipliers for

[t1,x1] = ode45(@odeeq_v5,0:1/12:80,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix);
X=x1';
% for 2000 to 2030 with 0:1/12:80 this is 601:961
% in this scale time2015 is 2015 whichis index 781 in above line
time2015 = 781;
% the compartments ordered in compvec, Table B3 Scott supplement
% scaline is scale10=40000/sum(X([chronic_nums10],66)); % current
% get individual multipliers
% ********** This is for current
qualy_vals10=zeros(11,1); % This conatins qualy for each chronic department
qualy_vals00=zeros(11,1); % This conatins qualy for each chronic department
testy = 0;
for i = 1 : 11
    Xtest=scale10*sum(X(compvec10(i,:),:));
    Ytest=scale00*sum(X(compvec00(i,:),:));
    %testy=testy+Xtest(time2015);
    qualy_vals10(i)=calc_qalyv2(Xtest(time2015:961),2015,2030,1/12,utility_vec(i),3); % 2000 to 2030 in steps of 1 month
    qualy_vals00(i)=calc_qalyv2(Ytest(time2015:961),2015,2030,1/12,utility_vec(i),3); % 2000 to 2030 in steps of 1 month
end
total_qual_current = sum(qualy_vals10)+sum(qualy_vals00);
death_vec= [-log(1-0.138) -log(1-0.605) -log(1-0.169) -log(1-0.034)]; %DC,HCC,LT1,LT2
%death_vec= [0.138 0.605 0.169 0.034]; %DC,HCC,LT1,LT2 these are already rates
% total deaths

% test for 12,13,14 and 15
testcompvec10= [DC_comps_10;HCC_comps_10;LT1_comps_10;LT2_comps_10];
testcompvec00= [DC_comps_00;HCC_comps_00;LT1_comps_00;LT2_comps_00];
death_vals10=zeros(4,1);
death_vals00=zeros(4,1);
Xdeaths = zeros(1,2030-2015);
for i = 1:4
    Xtest=scale10*sum(X(testcompvec10(i,:),:));
    Ytest=scale00*sum(X(testcompvec00(i,:),:));
    Xdeaths=Xdeaths+...
        just_integrate(Xtest(time2015:961),2015,2030,1/12)*death_vec(i)+...
        just_integrate(Ytest(time2015:961),2015,2030,1/12)*death_vec(i);
    death_vals10(i)=calc_qalyv2(Xtest(time2015:961),2015,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
    death_vals00(i)=calc_qalyv2(Ytest(time2015:961),2015,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
end
tot_liver_deaths=sum(death_vals00.*death_vec')+sum(death_vals10.*death_vec');
tot_liver_deaths_vec=death_vals10+death_vals10;

% incidence for compartments
% for acute duration is 12/52 years
% prev = inci*duration
% incidence = prev/duration
% incidence = prev * rate
rate=parm_current(2);
T00=rate*scale00*sum(X(A_comps_00,:));
T10=rate*scale10*sum(X(A_comps_10,:));
% integrate
inci_acute00=calc_qalyv2(T00(time2015:961),2015,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
inci_acute10=calc_qalyv2(T10(time2015:961),2015,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month

inci_acute=inci_acute00+inci_acute10;

% test to see if integration makes any differebce
% yes above is more accurate
[~,x1] = ode45(@odeeq_v5,0:80,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix);
X2=x1';
T002=rate*scale00*sum(X2(A_comps_00,:));
T102=rate*scale10*sum(X2(A_comps_10,:));
% integrate
inci_acute002=calc_qalyv2(T002(65:81),2015,2030,1,1,0); % 2000 to 2030 in steps of 1 month
inci_acute102=calc_qalyv2(T102(65:81),2015,2030,1,1,0); % 2000 to 2030 in steps of 1 month
inci_acute2=inci_acute002+inci_acute102;

% look at number of liver transplants
% DC to LT the rate is parm_names(13)
r_DCLT=parm_current(13);
T00=r_DCLT*scale00*sum(X(DC_comps_00,:));
T10=r_DCLT*scale10*sum(X(DC_comps_10,:));
% integrate
DCLT00=calc_qalyv2(T00(time2015:961),2015,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
DCLT10=calc_qalyv2(T10(time2015:961),2015,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month

% HCC to LT the rate is parm_names(15)
r_HCCLT=parm_current(15);
T00=r_HCCLT*scale00*sum(X(HCC_comps_00,:));
T10=r_HCCLT*scale10*sum(X(HCC_comps_10,:));
% integrate
HCCLT00=calc_qalyv2(T00(time2015:961),2015,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
HCCLT10=calc_qalyv2(T10(time2015:961),2015,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
total_liver_transplants = DCLT00 + DCLT10 + HCCLT00 + HCCLT10;

% compare with the prev of LT1
% rate is 1 as LT1 is the first year
T00=scale00*sum(X(LT1_comps_00,:));
T10=scale10*sum(X(LT1_comps_10,:));
% integrate
LT1T00=calc_qalyv2(T00(time2015:961),2015,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
LT1T10=calc_qalyv2(T10(time2015:961),2015,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
one_yearLT_compartment = LT1T00+LT1T10;
% will be less than above due to deaths

% calculate total costs
[costanual,costoneoff,namesof_anualcost,namesof_oneoff]=cost_params();
% anual costs
compvec10= [F0_comps_10;F1_comps_10;F2_comps_10;F3_comps_10;F4_comps_10;HCC_comps_10;DC_comps_10];
compvec00= [F0_comps_00;F1_comps_00;F2_comps_00;F3_comps_00;F4_comps_00;HCC_comps_00;DC_comps_00];
cost_vals10=zeros(7,1); % This conatins qualy for each chronic department
cost_vals00=zeros(7,1); % This conatins qualy for each chronic department
for i = 1 : 7
    Xtest=scale10*sum(X(compvec10(i,:),:));
    Ytest=scale00*sum(X(compvec00(i,:),:));
    %testy=testy+Xtest(time2015);
    cost_vals10(i)=calc_qalyv2(Xtest(time2015:961),2015,2030,1/12,costanual(i),3); % 2000 to 2030 in steps of 1 month
    cost_vals00(i)=calc_qalyv2(Ytest(time2015:961),2015,2030,1/12,costanual(i),3); % 2000 to 2030 in steps of 1 month
end
total_annual_costs = sum(cost_vals10)+sum(cost_vals00);

% total liver costs
% look at number of liver transplants
% DC to LT the rate is parm_names(13)
LTcost=costoneoff(3);
r_DCLT=parm_current(13);
T00=r_DCLT*scale00*sum(X(DC_comps_00,:));
T10=r_DCLT*scale10*sum(X(DC_comps_10,:));
% integrate
DCLT00=calc_qalyv2(T00(time2015:961),2015,2030,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
DCLT10=calc_qalyv2(T10(time2015:961),2015,2030,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month

% HCC to LT the rate is parm_names(15)
r_HCCLT=parm_current(15);
T00=r_HCCLT*scale00*sum(X(HCC_comps_00,:));
T10=r_HCCLT*scale10*sum(X(HCC_comps_10,:));
% integrate
HCCLT00=calc_qalyv2(T00(time2015:961),2015,2030,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
HCCLT10=calc_qalyv2(T10(time2015:961),2015,2030,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
total_cost_liver_transplants = DCLT00 + DCLT10 + HCCLT00 + HCCLT10;

% one off costs for diagnosis of incident chronic cases
rate=parm_current(2);
delta=parm_current(1);
T00=(1-delta)*rate*scale00*sum(X(A_comps_00,:));
T10=(1-delta)*rate*scale10*sum(X(A_comps_10,:));
% integrate
inci_acute00=calc_qalyv2(T00(time2015:961),2015,2030,1/12,costoneoff(4),3); % 2000 to 2030 in steps of 1 month
inci_acute10=calc_qalyv2(T10(time2015:961),2015,2030,1/12,costoneoff(4),3); % 2000 to 2030 in steps of 1 month
yearly_diagnosed_cases_in_chronic_stages=inci_acute00+inci_acute10;

% one off costs for transition incidence
% 'r_F3F4','r_F4DC','r_F4HCC'
% F3 to F4
r_F3F4=parm_current(9);
T00=r_F3F4*scale00*sum(X(F3_comps_00,:));
T10=r_F3F4*scale10*sum(X(F3_comps_10,:));
% integrate
F3F400=calc_qalyv2(T00(time2015:961),2015,2030,1/12,costoneoff(1),3); % 2000 to 2030 in steps of 1 month
F3F410=calc_qalyv2(T10(time2015:961),2015,2030,1/12,costoneoff(1),3); % 2000 to 2030 in steps of 1 month

% F4 to HCC
r_F4HCC=parm_current(11);
T00=r_F4HCC*scale00*sum(X(F4_comps_00,:));
T10=r_F4HCC*scale10*sum(X(F4_comps_10,:));
% integrate
F4HCC00=calc_qalyv2(T00(time2015:961),2015,2030,1/12,costoneoff(2),3); % 2000 to 2030 in steps of 1 month
F4HCC10=calc_qalyv2(T10(time2015:961),2015,2030,1/12,costoneoff(2),3); % 2000 to 2030 in steps of 1 month


total_transistion_costs = F3F400 + F3F410 + F4HCC00 + F4HCC10;


total_costs = total_annual_costs+total_cost_liver_transplants+ ...
    yearly_diagnosed_cases_in_chronic_stages+total_transistion_costs;


%% run with treatemnt
% This assumes the treatemnt is applied to current and former
% at advanced stages

ntreat=5662; % Number of treatments per year, 2015 and onwards
pcom=0.892; % prob of current PWID completing treatment
alpha=0.95; % prob of acheiving SVR
%[t1,x1t] = ode45(@odeeq_v5T,0:1/12:80,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix,...
%                scale00,scale10,ntreat,pcom,alpha,65); % note 2015 marker is 65 years
%XT=x1t';

% check that with p = 1 it is the same
p=1;
[t1,x1tp] = ode45(@odeeq_v5T_p,0:1/12:80,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix,...
                scale00,scale10,ntreat,pcom,alpha,65,p); % note 2015 marker is 65 years
XT=x1tp';
%sum(abs(XT(:)-XTp(:)))
% ** same
% make some checks when ntreat = 0
figure;plot(t1,100*sum(XT(chronic_prev_current,:))./sum(XT(total_current,:)),'r--','LineWidth',2)
set(gca,'ylim',[0 70])
xticks(0:5:80)
set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
title('% of current PWID with HCV from 1950 to 2030')
grid()

rate=parm_current(2);
delta=parm_current(1);
T00=(1-delta)*rate*scale00*sum(XT(A_comps_00,:));
T01=(1-delta)*rate*scale00*sum(XT(A_comps_01,:));
T11=(1-delta)*rate*scale10*sum(XT(A_comps_11,:));
T10=(1-delta)*rate*scale10*sum(XT(A_comps_10,:));
% integrate
yinci_PWID=just_integrate(T10(time2015:961),2015,2030,1/12)+...
           just_integrate(T00(time2015:961),2015,2030,1/12)+...
           just_integrate(T01(time2015:961),2015,2030,1/12)+...
           just_integrate(T11(time2015:961),2015,2030,1/12);
figure;plot(2015:2029,100*((yinci_PWID/yinci_PWID(1)) - 1),'r-','LineWidth',2);
grid()
title('% incidence change for treating advanced disease')
% deaths after treatment
% test for 12,13,14 and 15
death_vals10=zeros(4,1);
death_vals00=zeros(4,1);
death_vals11=zeros(4,1);
death_vals01=zeros(4,1);
XTdeaths=zeros(1,2030-2015);
% Now can have deaths in the failed treatment compartments
testcompvec01_failed=[192:20:360;193:20:360;194:20:360;195:20:360];
testcompvec11_failed=[552:20:720;553:20:720;554:20:720;555:20:720];
% death_vec= order is %DC,HCC,LT1,LT2 these are the probablities per year
for i = 1:4
    Xtest=scale10*sum(XT(testcompvec10(i,:),:));
    Ytest=scale00*sum(XT(testcompvec00(i,:),:));
    X1test=scale10*sum(XT(testcompvec11_failed(i,:),:));
    Y1test=scale00*sum(XT(testcompvec01_failed(i,:),:));
    XTdeaths=XTdeaths+...
        just_integrate(Xtest(time2015:961),2015,2030,1/12)*death_vec(i)+...
        just_integrate(Ytest(time2015:961),2015,2030,1/12)*death_vec(i)+...
        just_integrate(X1test(time2015:961),2015,2030,1/12)*death_vec(i)+...
        just_integrate(Y1test(time2015:961),2015,2030,1/12)*death_vec(i);
    death_vals10(i)=calc_qalyv2(Xtest(time2015:961),2015,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
    death_vals00(i)=calc_qalyv2(Ytest(time2015:961),2015,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
    death_vals11(i)=calc_qalyv2(X1test(time2015:961),2015,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
    death_vals01(i)=calc_qalyv2(Y1test(time2015:961),2015,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
end
tot_liver_deaths_trt=sum(death_vals00.*death_vec')+sum(death_vals10.*death_vec')+...
                     sum(death_vals01.*death_vec')+sum(death_vals11.*death_vec');
tot_liver_deaths_vec_trt=death_vals10+death_vals10;
%figure;plot(2015:2029,Xdeaths,'g-',2015:2029,XTdeaths,'r-','LineWidth',2);
%grid()
figure;plot(2015:2029,100*((XTdeaths/XTdeaths(1))),'r-','LineWidth',2);
grid()
title('% deaths change for treating advanced disease')


chronic_prev_current = [chronic_nums10]; % with no treatments chronic_nums11 all zero
total_current = [361:540];
figure;plot(t1,100*sum(XT(chronic_prev_current,:))./sum(XT(total_current,:)),'r--','LineWidth',2)
set(gca,'ylim',[0 70])
xticks(0:5:80)
set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
title('% of current PWID with HCV from 1950 to 2030')
grid()

%% now run with p = 0

ntreat=5725; % Number of treatments per year, 2015 and onwards
pcom=0.892; % prob of current PWID completing treatment
alpha=0.95; % prob of acheiving SVR
p=0;
[t1,x1t] = ode45(@odeeq_v5T_p,0:1/12:80,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix,...
                scale00,scale10,ntreat,pcom,alpha,65,p); % note 2015 marker is 65 years
XT=x1t';

% This reduces the number of incident cases

% one off costs for diagnosis of incident chronic cases



rate=parm_current(2);
delta=parm_current(1);
T00=(1-delta)*rate*scale00*sum(XT(A_comps_00,:));
T01=(1-delta)*rate*scale00*sum(XT(A_comps_01,:));
T11=(1-delta)*rate*scale10*sum(XT(A_comps_11,:));
T10=(1-delta)*rate*scale10*sum(XT(A_comps_10,:));
% integrate
yinci_PWID=just_integrate(T10(time2015:961),2015,2030,1/12)+...
           just_integrate(T00(time2015:961),2015,2030,1/12)+...
           just_integrate(T01(time2015:961),2015,2030,1/12)+...
           just_integrate(T11(time2015:961),2015,2030,1/12);
figure;plot(2015:2029,100*((yinci_PWID/yinci_PWID(1))),'r-','LineWidth',2);
grid()
title('% incidence change for treating current PWID')

% deaths
% deaths after treatment
% test for 12,13,14 and 15
death_vals10=zeros(4,1);
death_vals00=zeros(4,1);
death_vals11=zeros(4,1);
death_vals01=zeros(4,1);
XTdeaths=zeros(1,2030-2015);
% Now can have deaths in the failed treatment compartments
testcompvec01_failed=[192:20:360;193:20:360;194:20:360;195:20:360];
testcompvec11_failed=[552:20:720;553:20:720;554:20:720;555:20:720];
% death_vec= order is %DC,HCC,LT1,LT2 these are the probablities per year
for i = 1:4
    Xtest=scale10*sum(XT(testcompvec10(i,:),:));
    Ytest=scale00*sum(XT(testcompvec00(i,:),:));
    X1test=scale10*sum(XT(testcompvec11_failed(i,:),:));
    Y1test=scale00*sum(XT(testcompvec01_failed(i,:),:));
    XTdeaths=XTdeaths+...
        just_integrate(Xtest(time2015:961),2015,2030,1/12)*death_vec(i)+...
        just_integrate(Ytest(time2015:961),2015,2030,1/12)*death_vec(i)+...
        just_integrate(X1test(time2015:961),2015,2030,1/12)*death_vec(i)+...
        just_integrate(Y1test(time2015:961),2015,2030,1/12)*death_vec(i);
%        XTdeaths=XTdeaths+...
%        just_integrate(Xtest(time2015:961),2015,2030,1/12)*death_vec(i)+...
%        just_integrate(Ytest(time2015:961),2015,2030,1/12)*death_vec(i);
    death_vals10(i)=calc_qalyv2(Xtest(time2015:961),2015,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
    death_vals00(i)=calc_qalyv2(Ytest(time2015:961),2015,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
    death_vals11(i)=calc_qalyv2(X1test(time2015:961),2015,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
    death_vals01(i)=calc_qalyv2(Y1test(time2015:961),2015,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
end
tot_liver_deaths_trt=sum(death_vals00.*death_vec')+sum(death_vals10.*death_vec')+...
                     sum(death_vals01.*death_vec')+sum(death_vals11.*death_vec');
%figure;plot(2015:2029,Xdeaths,'g-',2015:2029,XTdeaths,'r-','LineWidth',2);
%grid()

figure;plot(2015:2029,100*((XTdeaths/XTdeaths(1))),'r-','LineWidth',2);
grid()
title('% deaths change for treating current PWID')








