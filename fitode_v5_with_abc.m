function fitode_v5_with_abc
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
% create_params_v5atabc sets max at 60% in 2002
% 50% in 2013

% create_params_v6atabc sets max at 56% in 2003
% 50% in 2013
% F0 to F1 61% in 2015 of PWID
% F3 to F4 11% in 2015 of PWID

[N0,M,Mdash,mort_current,mort_former,extra_parms_vals, phiminusvals, phiplussvals,extra_parms_nams,age_matrix,parm_names,parm_current,parm_former]= create_params_v5atabc;
% solve the diff equation
N0(1)=560;N0(361)=435.6;N0(367)=4.4;%4.36;
extra_parms_vals(1)=0.113;
[t1,x1] = ode45(@odeeq_v5,0:81,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix);
X=x1';
% chronic HCV prevlance calibrated at 50% for current PWID in 2015.
chronic_prev_current = [chronic_nums10]; % with no treatments chronic_nums11 all zero
total_current = [361:540];
figure;plot(t1,100*sum(X(chronic_prev_current,:))./sum(X(total_current,:)),'r--','LineWidth',2)
set(gca,'ylim',[0 70])
xticks(0:5:80)
set(gca,'xticklabel',1950:5:2031,'xticklabelrotation',90)
title('% of current PWID with HCV from 1950 to 2031')
set(gca,'xlim',[0 81])
grid()
% cols are f0,f1,f2,f3,f4,dc,hcc,lt1,lt2
% rows are the age
bot = sum(X(total_current,:));
agemat=[...
361 367   368   369   370   371   372   373   374   375	
381 387   388   389   390   391   392   393   394   395
401 407   408   409   410   411   412   413   414   415
421 427   428   429   430   431   432   433   434   435
441 447   448   449   450   451   452   453   454   455
461 467   468   469   470   471   472   473   474   475
481 487   488   489   490   491   492   493   494   495
501 507   508   509   510   511   512   513   514   515
521 527   528   529   530   531   532   533   534   535];
mn=size(X);
ageprop=zeros(9,mn(2));
cols  = lines(9);  
H = zeros(9,1);  % Store handles
figure;
for i = 1 : 9
    ageprop(i,:) = sum(X(agemat(i,:),:));
    H(i)=plot(t1,100*ageprop(i,:)./bot,'-','LineWidth',2,'Color',cols(i,:))
    hold on
end
set(gca,'ylim',[0 25])
xticks(0:5:80)
set(gca,'xticklabel',1950:5:2031,'xticklabelrotation',90)
title('% of current PWID with HCV from 1950 to 2031 by age')
set(gca,'xlim',[0 81])
grid()
legend(H, {'20-24','25-29','30-34','35-44','45-54','55-64','65-74','75-84','85+'},'location','best');

H = zeros(9,1);  % Store handles
age_weights=zeros(9,1);

figure;
for i = 1 : 9
    if i == 1
        age_weights(i) = 0.5*0.25*bot(66)/ageprop(1,66);
    else
        age_weights(i) = 0.5*(0.75/8)*bot(66)/ageprop(i,66);
    end
    ageprop(i,:) = sum(X(agemat(i,:),:));
    H(i)=plot(t1,100*age_weights(i)*ageprop(i,:)./bot,'-','LineWidth',2,'Color',cols(i,:))
    hold on
end
A=repmat(repelem(age_weights,9),1,82);
set(gca,'ylim',[0 30])
xticks(0:5:80)
set(gca,'xticklabel',1950:5:2031,'xticklabelrotation',90)
title('% of current PWID with HCV from 1950 to 2031 age weighted')
set(gca,'xlim',[0 81])
grid()
legend(H, {'20-24','25-29','30-34','35-44','45-54','55-64','65-74','75-84','85+'},'location','best');
%age_weights_matrix
figure;plot(t1,100*sum(A.*X(chronic_prev_current,:))./sum(X(total_current,:)),'r--','LineWidth',2)
set(gca,'ylim',[0 70])
xticks(0:5:80)
set(gca,'xticklabel',1950:5:2031,'xticklabelrotation',90)
title('% of current PWID with HCV from 1950 to 2031')
set(gca,'xlim',[0 81])
grid()

% what do the distribution for former look like

agemat_former=[...
361 367   368   369   370   371   372   373   374   375	
381 387   388   389   390   391   392   393   394   395
401 407   408   409   410   411   412   413   414   415
421 427   428   429   430   431   432   433   434   435
441 447   448   449   450   451   452   453   454   455
461 467   468   469   470   471   472   473   474   475
481 487   488   489   490   491   492   493   494   495
501 507   508   509   510   511   512   513   514   515
521 527   528   529   530   531   532   533   534   535]-360;
total_former=1:360;
bot = sum(X(total_former,:));
ageprop=zeros(9,mn(2));
cols  = lines(9);  
H = zeros(9,1);  % Store handles
figure;
for i = 1 : 9
    ageprop(i,:) = sum(X(agemat_former(i,:),:));
    H(i)=plot(t1,100*ageprop(i,:)./bot,'-','LineWidth',2,'Color',cols(i,:))
    hold on
end
set(gca,'ylim',[0 25])
xticks(0:5:80)
set(gca,'xticklabel',1950:5:2031,'xticklabelrotation',90)
title('% of former PWID with HCV from 1950 to 2031 by age')
set(gca,'xlim',[0 81])
grid()
legend(H, {'20-24','25-29','30-34','35-44','45-54','55-64','65-74','75-84','85+'},'location','best');

% age_weights for former
H = zeros(9,1);  % Store handles
age_weights_former=zeros(9,1);
for i = 1 : 9
    if i == 1
        age_weights_former(i) = 0.125*bot(66)/ageprop(1,66);
    else
        age_weights_former(i) = (0.875/8)*bot(66)/ageprop(i,66);
    end
end
A_former=repmat(repelem(age_weights_former,9),1,82);



%chronic_prev = [chronic_nums00 chronic_nums10];
%total_all = 1:720;
%figure;plot(t1,100*sum(X(chronic_prev,:))./sum(X(total_all,:)),'r--','LineWidth',2)
%set(gca,'ylim',[0 70])
%xticks(0:5:80)
%set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
%title('% of all PWID with HCV from 1950 to 2030')
%grid()

scale00=144000/sum(X([chronic_nums00],66)); % former
scale10=40000/sum(X([chronic_nums10],66)); % current
%scale10=80000/sum(X([361:540],66));


figure;plot(t1,scale00*sum(X([chronic_nums00],:))+scale10*sum(X([chronic_nums10],:)))
xticks(0:5:81)
set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
title('Absolute numbers of PWID with chronic HCV, scaled to the 2015 estimates')
set(gca,'xlim',[0 81])
grid()
% scaled
age_scale00=144000/sum(A_former(:,66).*X([chronic_nums00],66));
age_scale10=40000/sum(A(:,66).*X([chronic_nums10],66)); % current
figure;plot(t1,age_scale00*sum(A_former.*X(chronic_nums00,:))+age_scale10*sum(A.*X(chronic_nums10,:)))
xticks(0:5:81)
set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
title('Absolute numbers of PWID with chronic HCV, scaled to the 2015 estimates and age weighted')
set(gca,'xlim',[0 81])
grid()

%figure; plot(t1,scale10*sum(X([chronic_nums10],:))+scale10*sum(X(S0_comps(19:27),:)))
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
xticks(0:5:81)
set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
legend({'lt' 'hcc' 'dc' 'f4' 'f3' 'f2' 'f1' 'f0'},'location','best');
grid()
title('Number of PWID at each stage of chronic HCV')
set(gca,'xlim',[0 81])

% scaled version
Apwid=repmat(age_weights,1,82);
Aformer=repmat(age_weights_former,1,82);
lt2=age_scale00*sum(Aformer.*X(LT2_comps_00,:))+age_scale10*sum(Apwid.*X(LT2_comps_10,:));
lt1=age_scale00*sum(Aformer.*X(LT1_comps_00,:))+age_scale10*sum(Apwid.*X(LT1_comps_10,:));
lt=lt1+lt2;
dc=age_scale00*sum(Aformer.*X(DC_comps_00,:))+age_scale10*sum(Apwid.*X(DC_comps_10,:));
hcc=age_scale00*sum(Aformer.*X(HCC_comps_00,:))+age_scale10*sum(Apwid.*X(HCC_comps_10,:));
f4=age_scale00*sum(Aformer.*X(F4_comps_00,:))+age_scale10*sum(Apwid.*X(F4_comps_10,:));
f3=age_scale00*sum(Aformer.*X(F3_comps_00,:))+age_scale10*sum(Apwid.*X(F3_comps_10,:));
f2=age_scale00*sum(Aformer.*X(F2_comps_00,:))+age_scale10*sum(Apwid.*X(F2_comps_10,:));
f1=age_scale00*sum(Aformer.*X(F1_comps_00,:))+age_scale10*sum(Apwid.*X(F1_comps_10,:));
f0=age_scale00*sum(Aformer.*X(F0_comps_00,:))+age_scale10*sum(Apwid.*X(F0_comps_10,:));


figure;area(t1,[lt' hcc' dc' f4' f3' f2' f1' f0'])
xticks(0:5:81)
set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
legend({'lt' 'hcc' 'dc' 'f4' 'f3' 'f2' 'f1' 'f0'},'location','best');
grid()
title('Number of PWID at each stage of chronic HCV')
set(gca,'xlim',[0 81])





f=sum(X(chronic_prev_current,:))./sum(X(total_current,:));
sum2=f(53); % target 0.6
sum3=sum(X(chronic_prev_current,66))/sum(X(total_current,66)); % target 0.5sum4= (f0(66)+f1(66))/(f0(66)+f1(66)+f2(66)+f3(66)+f4(66)+hcc(66)+dc(66)+lt(66))
sum4= (f0(66)+f1(66))/(f0(66)+f1(66)+f2(66)+f3(66)+f4(66)+hcc(66)+dc(66)+lt(66)); % target = 0.56

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

[t1,x1] = ode45(@odeeq_v5,0:1/12:81,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix);
X=x1';
% for 2000 to 2031 with 0:1/12:80 this is 601:973
% in this scale time2015 is 2015 whichis index 781 in above line
time2016 = 793; % Nite this is 2016   was 781
% the compartments ordered in compvec, Table B3 Scott supplement
% scaline is scale10=40000/sum(X([chronic_nums10],66)); % current
% get individual multipliers
% ********** This is for current

qualy_vals10=zeros(11,1); % This conatins PWID qualy for each chronic department
qualy_vals00=zeros(11,1); % This conatins former qualy for each chronic department

for i = 1 : 11
    Xtest=scale10*sum(X(compvec10(i,:),:));
    Ytest=scale00*sum(X(compvec00(i,:),:));
    %testy=testy+Xtest(time2015);
    qualy_vals10(i)=calc_qalyv2(Xtest(time2016:973),2016,2031,1/12,utility_vec(i),3); % 2000 to 2030 in steps of 1 month
    qualy_vals00(i)=calc_qalyv2(Ytest(time2016:973),2016,2031,1/12,utility_vec(i),3); % 2000 to 2030 in steps of 1 month
end
total_qual_current = sum(qualy_vals10)+sum(qualy_vals00);

% above with age weighting
% test thsi should be 184,000
%  scale10*sum(X(compvec10(3:end,:),time2015-12))+scale00*sum(X(compvec00(3:end,:),time2015-12))
% Apwidtest=repmat(repelem(age_weights,9),1,973);
% Aformertest=repmat(repelem(age_weights_former,9),1,973);
% age_scale10*sum(Apwidtest(:,1).*X(compvec10(3:end,:),time2015-12))+age_scale00*sum(Aformertest(:,1).*X(compvec00(3:end,:),time2015-12))

qualy_vals10_scaled=zeros(11,1); % This conatins PWID qualy for each chronic department
qualy_vals00_scaled=zeros(11,1); % This conatins former qualy for each chronic department
Apwid=repmat(repelem(age_weights,1),1,973);
Aformer=repmat(repelem(age_weights_former,1),1,973);
for i = 1 : 11
    Xtest_scaled=age_scale10*sum(Apwid.*X(compvec10(i,:),:));
    Ytest_scaled=age_scale00*sum(Aformer.*X(compvec00(i,:),:));
    %testy=testy+Xtest(time2015);
    qualy_vals10_scaled(i)=calc_qalyv2(Xtest_scaled(time2016:973),2016,2031,1/12,utility_vec(i),3); % 2000 to 2030 in steps of 1 month
    qualy_vals00_scaled(i)=calc_qalyv2(Ytest_scaled(time2016:973),2016,2031,1/12,utility_vec(i),3); % 2000 to 2030 in steps of 1 month
end
total_qual_current_scaled = sum(qualy_vals10_scaled)+sum(qualy_vals00_scaled);



% add in the fatality deaths
death_vec= [-log(1-0.138) -log(1-0.605) -log(1-0.169) -log(1-0.034)]; %DC,HCC,LT1,LT2
%death_vec= [0.138 0.605 0.169 0.034]; %DC,HCC,LT1,LT2 these are already rates
% total deaths
time2015=673;
% test for 12,13,14 and 15
testcompvec10= [DC_comps_10;HCC_comps_10;LT1_comps_10;LT2_comps_10];
testcompvec00= [DC_comps_00;HCC_comps_00;LT1_comps_00;LT2_comps_00];
death_vals10=zeros(4,1);
death_vals00=zeros(4,1);
Xdeaths = zeros(1,2030-2016+1);
m=size(X);
Deathmat=zeros(9,m(2));
for i = 1:4
    Deathmat=Deathmat + (scale10*X(testcompvec10(i,:),:)+scale00*X(testcompvec00(i,:),:))*death_vec(i);
    Xtest=scale10*sum(X(testcompvec10(i,:),:));
    Ytest=scale00*sum(X(testcompvec00(i,:),:));
    Xdeaths=Xdeaths+...
        just_integrate(Xtest(time2015:973),2016,2031,1/12)*death_vec(i)+...
        just_integrate(Ytest(time2015:973),2016,2031,1/12)*death_vec(i);
    death_vals10(i)=calc_qalyv2(Xtest(time2015:973),2016,2031,1/12,1,0); % 2000 to 2030 in steps of 1 month
    death_vals00(i)=calc_qalyv2(Ytest(time2015:973),2016,2031,1/12,1,0); % 2000 to 2030 in steps of 1 month
end
tot_liver_deaths=sum(death_vals00.*death_vec')+sum(death_vals10.*death_vec');

% tot liver deaths with age weighting
Apwid=repmat(repelem(age_weights,1),1,973);
Aformer=repmat(repelem(age_weights_former,1),1,973);
death_vals10_scaled=zeros(4,1);
death_vals00_scaled=zeros(4,1);
Xdeaths_scaled = zeros(1,2030-2016+1);
m=size(X);
Deathmat_scaled=zeros(9,m(2));
for i = 1:4
    Deathmat_scaled=Deathmat_scaled + (age_scale10*Apwid.*X(testcompvec10(i,:),:)+scale00*Aformer.*X(testcompvec00(i,:),:))*death_vec(i);
    Xtest=age_scale10*sum(Apwid.*X(testcompvec10(i,:),:));
    Ytest=age_scale00*sum(Aformer.*X(testcompvec00(i,:),:));
    Xdeaths_scaled=Xdeaths_scaled+...
        just_integrate(Xtest(time2015:973),2016,2031,1/12)*death_vec(i)+...
        just_integrate(Ytest(time2015:973),2016,2031,1/12)*death_vec(i);
    death_vals10_scaled(i)=calc_qalyv2(Xtest(time2015:973),2016,2031,1/12,1,0); % 2000 to 2030 in steps of 1 month
    death_vals00_scaled(i)=calc_qalyv2(Ytest(time2015:973),2016,2031,1/12,1,0); % 2000 to 2030 in steps of 1 month
end
tot_liver_deaths_scaled=sum(death_vals00_scaled.*death_vec')+sum(death_vals10_scaled.*death_vec');
t=1950:1/12:2031;
D=zeros(9,2031-2016+1);
for i = 1 : 9
    k=0;
    for j = 2016:2031
        jf=find(t==j):find(t==(j+1));
        x=t(jf);
        y=Deathmat(i,jf);
        k =k + 1;
        D(i,k) = trapz(x,y);
    end
end
Q_tot_liver_deaths=qalydeaths(D,3,2016,2030); % fatal qualys
% scaled of above
t=1950:1/12:2031;
D_scaled=zeros(9,2031-2016+1);
for i = 1 : 9
    k=0;
    for j = 2016:2031
        jf=find(t==j):find(t==(j+1));
        x=t(jf);
        y=Deathmat_scaled(i,jf);
        k =k + 1;
        D_scaled(i,k) = trapz(x,y);
    end
end
Q_tot_liver_deaths_scaled=qalydeaths(D_scaled,3,2016,2030); % fatal qualys

% incidence for compartments
% for acute duration is 12/52 years
% prev = inci*duration
% incidence = prev/duration
% incidence = prev * rate
rate=parm_current(2);
T00=rate*scale00*sum(X(A_comps_00,:));
T10=rate*scale10*sum(X(A_comps_10,:));
% integrate
inci_acute00=calc_qalyv2(T00(time2015:973),2016,2031,1/12,1,0); % 2000 to 2030 in steps of 1 month
inci_acute10=calc_qalyv2(T10(time2015:973),2016,2031,1/12,1,0); % 2000 to 2030 in steps of 1 month

inci_acute=inci_acute00+inci_acute10;
iout= just_integrate(T10(time2015:973),2016,2031,1/12)+just_integrate(T00(time2015:973),2016,2031,1/12);

% test to see if integration makes any differebce
% yes above is more accurate
%[~,x1] = ode45(@odeeq_v5,0:80,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix);
%X2=x1';
%T002=rate*scale00*sum(X2(A_comps_00,:));
%T102=rate*scale10*sum(X2(A_comps_10,:));
% integrate
%inci_acute002=calc_qalyv2(T002(65:81),2015,2030,1,1,0); % 2000 to 2030 in steps of 1 month
%inci_acute102=calc_qalyv2(T102(65:81),2015,2030,1,1,0); % 2000 to 2030 in steps of 1 month
%inci_acute2=inci_acute002+inci_acute102;

% look at number of liver transplants
% DC to LT the rate is parm_names(13)
r_DCLT=parm_current(13);
T00=r_DCLT*scale00*sum(X(DC_comps_00,:));
T10=r_DCLT*scale10*sum(X(DC_comps_10,:));
% integrate
DCLT00=calc_qalyv2(T00(time2015:973),2016,2031,1/12,1,0); % 2000 to 2030 in steps of 1 month
DCLT10=calc_qalyv2(T10(time2015:973),2016,2031,1/12,1,0); % 2000 to 2030 in steps of 1 month

% HCC to LT the rate is parm_names(15)
r_HCCLT=parm_current(15);
T00=r_HCCLT*scale00*sum(X(HCC_comps_00,:));
T10=r_HCCLT*scale10*sum(X(HCC_comps_10,:));
% integrate
HCCLT00=calc_qalyv2(T00(time2015:973),2016,2031,1/12,1,0); % 2000 to 2030 in steps of 1 month
HCCLT10=calc_qalyv2(T10(time2015:973),2016,2031,1/12,1,0); % 2000 to 2030 in steps of 1 month
total_liver_transplants = DCLT00 + DCLT10 + HCCLT00 + HCCLT10;

% compare with the prev of LT1
% rate is 1 as LT1 is the first year
T00=scale00*sum(X(LT1_comps_00,:));
T10=scale10*sum(X(LT1_comps_10,:));
% integrate
LT1T00=calc_qalyv2(T00(time2015:973),2016,2031,1/12,1,0); % 2000 to 2030 in steps of 1 month
LT1T10=calc_qalyv2(T10(time2015:973),2016,2031,1/12,1,0); % 2000 to 2030 in steps of 1 month
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
    cost_vals10(i)=calc_qalyv2(Xtest(time2015:973),2016,2031,1/12,costanual(i),3); % 2000 to 2030 in steps of 1 month
    cost_vals00(i)=calc_qalyv2(Ytest(time2015:973),2016,2031,1/12,costanual(i),3); % 2000 to 2030 in steps of 1 month
end
total_annual_costs = sum(cost_vals10)+sum(cost_vals00);
%As above for scaled
cost_vals10_scaled=zeros(7,1); % This conatins qualy for each chronic department
cost_vals00_scaled=zeros(7,1); % This conatins qualy for each chronic department
for i = 1 : 7
    Xtest=age_scale10*sum(Apwid.*X(compvec10(i,:),:));
    Ytest=age_scale00*sum(Aformer.*X(compvec00(i,:),:));
    %testy=testy+Xtest(time2015);
    cost_vals10_scaled(i)=calc_qalyv2(Xtest(time2015:973),2016,2031,1/12,costanual(i),3); % 2000 to 2030 in steps of 1 month
    cost_vals00_scaled(i)=calc_qalyv2(Ytest(time2015:973),2016,2031,1/12,costanual(i),3); % 2000 to 2030 in steps of 1 month
end
total_annual_costs_scaled = sum(cost_vals10_scaled)+sum(cost_vals00_scaled);



% total liver costs
% look at number of liver transplants
% DC to LT the rate is parm_names(13)
LTcost=costoneoff(3);
r_DCLT=parm_current(13);
T00=r_DCLT*scale00*sum(X(DC_comps_00,:));
T10=r_DCLT*scale10*sum(X(DC_comps_10,:));
% integrate
DCLT00=calc_qalyv2(T00(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
DCLT10=calc_qalyv2(T10(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month

% HCC to LT the rate is parm_names(15)
r_HCCLT=parm_current(15);
T00=r_HCCLT*scale00*sum(X(HCC_comps_00,:));
T10=r_HCCLT*scale10*sum(X(HCC_comps_10,:));
% integrate
HCCLT00=calc_qalyv2(T00(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
HCCLT10=calc_qalyv2(T10(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
total_cost_liver_transplants = DCLT00 + DCLT10 + HCCLT00 + HCCLT10;
% as above scaled
Apwid=repmat(repelem(age_weights,1),1,973);
Aformer=repmat(repelem(age_weights_former,1),1,973);
T00_scaled=r_DCLT*age_scale00*sum(Aformer.*X(DC_comps_00,:));
T10_scaled=r_DCLT*age_scale10*sum(Apwid.*X(DC_comps_10,:));
% integrate
DCLT00_scaled=calc_qalyv2(T00_scaled(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
DCLT10_scaled=calc_qalyv2(T10_scaled(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month

% HCC to LT the rate is parm_names(15)
r_HCCLT=parm_current(15);
T00_scaled=r_HCCLT*age_scale00*sum(Aformer.*X(HCC_comps_00,:));
T10_scaled=r_HCCLT*age_scale10*sum(Apwid.*X(HCC_comps_10,:));
% integrate
HCCLT00_scaled=calc_qalyv2(T00_scaled(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
HCCLT10_scaled=calc_qalyv2(T10_scaled(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
total_cost_liver_transplants_scaled = DCLT00_scaled + DCLT10_scaled + HCCLT00_scaled + HCCLT10_scaled;

% one off costs for diagnosis of incident chronic cases
rate=parm_current(2);
delta=parm_current(1);
T00=(1-delta)*rate*scale00*sum(X(A_comps_00,:));
T10=(1-delta)*rate*scale10*sum(X(A_comps_10,:));
% integrate
inci_acute00=calc_qalyv2(T00(time2015:973),2016,2031,1/12,costoneoff(4),3); % 2000 to 2030 in steps of 1 month
inci_acute10=calc_qalyv2(T10(time2015:973),2016,2031,1/12,costoneoff(4),3); % 2000 to 2030 in steps of 1 month
yearly_diagnosed_cases_in_chronic_stages=inci_acute00+inci_acute10;
% above for scaled
T00_scaled=(1-delta)*rate*age_scale00*sum(Aformer.*X(A_comps_00,:));
T10_scaled=(1-delta)*rate*age_scale10*sum(Apwid.*X(A_comps_10,:));
% integrate
inci_acute00_scaled=calc_qalyv2(T00_scaled(time2015:973),2016,2031,1/12,costoneoff(4),3); % 2000 to 2030 in steps of 1 month
inci_acute10_scaled=calc_qalyv2(T10_scaled(time2015:973),2016,2031,1/12,costoneoff(4),3); % 2000 to 2030 in steps of 1 month
yearly_diagnosed_cases_in_chronic_stages_scaled=inci_acute00_scaled+inci_acute10_scaled;

% one off costs for transition incidence
% 'r_F3F4','r_F4DC','r_F4HCC'
% F3 to F4
r_F3F4=parm_current(9);
T00=r_F3F4*scale00*sum(X(F3_comps_00,:));
T10=r_F3F4*scale10*sum(X(F3_comps_10,:));
% integrate
F3F400=calc_qalyv2(T00(time2015:973),2016,2031,1/12,costoneoff(1),3); % 2000 to 2030 in steps of 1 month
F3F410=calc_qalyv2(T10(time2015:973),2016,2031,1/12,costoneoff(1),3); % 2000 to 2030 in steps of 1 month
% F4 to HCC
r_F4HCC=parm_current(11);
T00=r_F4HCC*scale00*sum(X(F4_comps_00,:));
T10=r_F4HCC*scale10*sum(X(F4_comps_10,:));
% integrate
F4HCC00=calc_qalyv2(T00(time2015:973),2016,2031,1/12,costoneoff(2),3); % 2000 to 2030 in steps of 1 month
F4HCC10=calc_qalyv2(T10(time2015:973),2016,2031,1/12,costoneoff(2),3); % 2000 to 2030 in steps of 1 month
total_transistion_costs = F3F400 + F3F410 + F4HCC00 + F4HCC10;
% as above scaled
T00_scaled=r_F3F4*age_scale00*sum(Aformer.*X(F3_comps_00,:));
T10_scaled=r_F3F4*age_scale10*sum(Apwid.*X(F3_comps_10,:));
% integrate
F3F400_scaled=calc_qalyv2(T00_scaled(time2015:973),2016,2031,1/12,costoneoff(1),3); % 2000 to 2030 in steps of 1 month
F3F410_scaled=calc_qalyv2(T10_scaled(time2015:973),2016,2031,1/12,costoneoff(1),3); % 2000 to 2030 in steps of 1 month
% F4 to HCC
r_F4HCC=parm_current(11);
T00_scaled=r_F4HCC*age_scale00*sum(Aformer.*X(F4_comps_00,:));
T10_scaled=r_F4HCC*age_scale10*sum(Apwid.*X(F4_comps_10,:));
% integrate
F4HCC00_scaled=calc_qalyv2(T00_scaled(time2015:973),2016,2031,1/12,costoneoff(2),3); % 2000 to 2030 in steps of 1 month
F4HCC10_scaled=calc_qalyv2(T10_scaled(time2015:973),2016,2031,1/12,costoneoff(2),3); % 2000 to 2030 in steps of 1 month
total_transistion_costs_scaled = F3F400_scaled + F3F410_scaled + F4HCC00_scaled + F4HCC10_scaled;


total_costs = total_annual_costs+total_cost_liver_transplants+ ...
    yearly_diagnosed_cases_in_chronic_stages+total_transistion_costs;
total_costs_scaled = total_annual_costs_scaled+total_cost_liver_transplants_scaled+ ...
    yearly_diagnosed_cases_in_chronic_stages_scaled+total_transistion_costs_scaled;

total_qual_current =total_qual_current + Q_tot_liver_deaths;
total_qual_current_scaled =total_qual_current_scaled + Q_tot_liver_deaths_scaled;

import java.text.*
v=DecimalFormat;
outy{1,1}='No treatment';
outy{1,2}=['$' char(v.format(round(total_costs)))];
outy{1,3}=char(v.format(round(total_qual_current)));
outy{1,4}=char(v.format(round(tot_liver_deaths)));
outy{1,5}='';

outy_scaled{1,1}='No treatment';
outy_scaled{1,2}=['$' char(v.format(round(total_costs_scaled)))];
outy_scaled{1,3}=char(v.format(round(total_qual_current_scaled)));
outy_scaled{1,4}=char(v.format(round(tot_liver_deaths_scaled)));
outy_scaled{1,5}='';

cell2table(outy,'VariableNames',{'Intervention' 'TotalCost','TotalQALYS' 'LiverRelatedDeaths' 'ICER'})
cell2table(outy_scaled,'VariableNames',{'Intervention' 'TotalCost','TotalQALYS' 'LiverRelatedDeaths' 'ICER'})
%% run with treatemnt
% This assumes the treatemnt is applied to current and former
% at advanced stages

ntreat=5952; % Number of treatments per year, 2015 and onwards
pcom=0.892; % prob of current PWID completing treatment
alpha=0.95; % prob of acheiving SVR
%[t1,x1t] = ode45(@odeeq_v5T,0:1/12:80,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix,...
%                scale00,scale10,ntreat,pcom,alpha,65); % note 2015 marker is 65 years
%XT=x1t';

% check that with p = 1 it is the same
p=1;
[t1,x1tp] = ode45(@odeeq_v5T_p,0:1/12:81,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix,...
    scale00,scale10,ntreat,pcom,alpha,66,p); % note 2015 marker is 66 years
XT=x1tp';
%sum(abs(XT(:)-XTp(:)))
% ** same
% make some checks when ntreat = 0
figure;plot(t1,100*sum(XT(chronic_prev_current,:))./sum(XT(total_current,:)),'r--','LineWidth',2)
set(gca,'ylim',[0 70])
xticks(0:5:81)
set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
title('% of current PWID with HCV from 1950 to 2030')
set(gca,'xlim',[0 81])
grid()

rate=parm_current(2);
delta=parm_current(1);
T00=(1-delta)*rate*scale00*sum(XT(A_comps_00,:));
T01=(1-delta)*rate*scale00*sum(XT(A_comps_01,:));
T11=(1-delta)*rate*scale10*sum(XT(A_comps_11,:));
T10=(1-delta)*rate*scale10*sum(XT(A_comps_10,:));
% integrate
yinci_PWID=just_integrate(T10(time2015:973),2016,2031,1/12)+...
    just_integrate(T00(time2015:973),2016,2031,1/12)+...
    just_integrate(T01(time2015:973),2016,2031,1/12)+...
    just_integrate(T11(time2015:973),2016,2031,1/12);
figure;plot(2016:2030,100*((yinci_PWID/yinci_PWID(1)) - 1),'r-','LineWidth',2);
grid()
title('% incidence change for treating advanced disease')
% as above for age scaled
T00_scaled=(1-delta)*rate*age_scale00*sum(Aformer.*XT(A_comps_00,:));
T01_scaled=(1-delta)*rate*age_scale00*sum(Aformer.*XT(A_comps_01,:));
T11_scaled=(1-delta)*rate*age_scale10*sum(Apwid.*XT(A_comps_11,:));
T10_scaled=(1-delta)*rate*age_scale10*sum(Apwid.*XT(A_comps_10,:));
% integrate
yinci_PWID_scaled=just_integrate(T10_scaled(time2015:973),2016,2031,1/12)+...
    just_integrate(T00_scaled(time2015:973),2016,2031,1/12)+...
    just_integrate(T01_scaled(time2015:973),2016,2031,1/12)+...
    just_integrate(T11_scaled(time2015:973),2016,2031,1/12);
figure;plot(2016:2030,100*((yinci_PWID_scaled/yinci_PWID_scaled(1)) - 1),'r-','LineWidth',2);
grid()
title('% incidence change for treating advanced disease scaled')

% total subjects with advanced liver disease
%ad00=scale00*sum(XT([F3_comps_00 F4_comps_00 DC_comps_00 HCC_comps_00 LT1_comps_00 LT2_comps_00],:));
%ad01=scale00*sum(XT([F3_comps_01 F4_comps_01 DC_comps_01 HCC_comps_01 LT1_comps_01 LT2_comps_01],:));
%ad10=scale10*sum(XT([F3_comps_10 F4_comps_10 DC_comps_10 HCC_comps_10 LT1_comps_10 LT2_comps_10],:));
%ad11=scale10*sum(XT([F3_comps_11 F4_comps_11 DC_comps_11 HCC_comps_11 LT1_comps_11 LT2_comps_11],:));
ad00=scale00*sum(XT([F0_comps_00 F1_comps_00  F2_comps_00 F3_comps_00 F4_comps_00 DC_comps_00 HCC_comps_00 LT1_comps_00 LT2_comps_00],:));
ad01=scale00*sum(XT([F3_comps_01 F4_comps_01  LT1_comps_01 LT2_comps_01],:));
ad10=scale10*sum(XT([F0_comps_10 F1_comps_10  F2_comps_10 F3_comps_00 F4_comps_10 DC_comps_10 HCC_comps_10 LT1_comps_10 LT2_comps_10],:));
ad11=scale10*sum(XT([F3_comps_11 F4_comps_11  LT1_comps_11 LT2_comps_11],:));


advanced_tot=just_integrate(ad00(time2015:973),2016,2031,1/12)+...
    just_integrate(ad10(time2015:973),2016,2031,1/12);

ratio_of_treat_toadvanced = 15*ntreat/sum(advanced_tot);
% This is the proportion of IDU - aquired ifections that are treated


% deaths after treatment
% test for 12,13,14 and 15
death_vals10=zeros(4,1);
death_vals00=zeros(4,1);
death_vals11=zeros(4,1);
death_vals01=zeros(4,1);
XTdeaths=zeros(1,2030-2016+1);
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
        just_integrate(Xtest(time2015:973),2016,2031,1/12)*death_vec(i)+...
        just_integrate(Ytest(time2015:973),2016,2031,1/12)*death_vec(i)+...
        just_integrate(X1test(time2015:973),2016,2031,1/12)*death_vec(i)+...
        just_integrate(Y1test(time2015:973),2016,2031,1/12)*death_vec(i);
    death_vals10(i)=calc_qalyv2(Xtest(time2015:973),2016,2031,1/12,1,0); % 2000 to 2030 in steps of 1 month
    death_vals00(i)=calc_qalyv2(Ytest(time2015:973),2016,2031,1/12,1,0); % 2000 to 2030 in steps of 1 month
    death_vals11(i)=calc_qalyv2(X1test(time2015:973),2016,2031,1/12,1,0); % 2000 to 2030 in steps of 1 month
    death_vals01(i)=calc_qalyv2(Y1test(time2015:973),2016,2031,1/12,1,0); % 2000 to 2030 in steps of 1 month
end
tot_liver_deaths_trt=sum(death_vals00.*death_vec')+sum(death_vals10.*death_vec')+...
    sum(death_vals01.*death_vec')+sum(death_vals11.*death_vec');

%figure;plot(2015:2029,Xdeaths,'g-',2015:2029,XTdeaths,'r-','LineWidth',2);
%grid()
figure;plot(2016:2030,100*((XTdeaths/XTdeaths(1))),'r-','LineWidth',2);
grid()
title('% deaths change for treating advanced disease')

% as above for scaled age
% deaths after treatment
% test for 12,13,14 and 15
death_vals10_scaled=zeros(4,1);
death_vals00_scaled=zeros(4,1);
death_vals11_scaled=zeros(4,1);
death_vals01_scaled=zeros(4,1);
XTdeaths_scaled=zeros(1,2030-2016+1);
% Now can have deaths in the failed treatment compartments
testcompvec01_failed=[192:20:360;193:20:360;194:20:360;195:20:360];
testcompvec11_failed=[552:20:720;553:20:720;554:20:720;555:20:720];
% death_vec= order is %DC,HCC,LT1,LT2 these are the probablities per year
for i = 1:4
    Xtest_scaled=age_scale10*sum(Apwid.*XT(testcompvec10(i,:),:));
    Ytest_scaled=age_scale00*sum(Aformer.*XT(testcompvec00(i,:),:));
    X1test_scaled=age_scale10*sum(Apwid.*XT(testcompvec11_failed(i,:),:));
    Y1test_scaled=age_scale00*sum(Aformer.*XT(testcompvec01_failed(i,:),:));
    XTdeaths_scaled=XTdeaths_scaled+...
        just_integrate(Xtest_scaled(time2015:973),2016,2031,1/12)*death_vec(i)+...
        just_integrate(Ytest_scaled(time2015:973),2016,2031,1/12)*death_vec(i)+...
        just_integrate(X1test_scaled(time2015:973),2016,2031,1/12)*death_vec(i)+...
        just_integrate(Y1test_scaled(time2015:973),2016,2031,1/12)*death_vec(i);
    death_vals10_scaled(i)=calc_qalyv2(Xtest_scaled(time2015:973),2016,2031,1/12,1,0); % 2000 to 2030 in steps of 1 month
    death_vals00_scaled(i)=calc_qalyv2(Ytest_scaled(time2015:973),2016,2031,1/12,1,0); % 2000 to 2030 in steps of 1 month
    death_vals11_scaled(i)=calc_qalyv2(X1test_scaled(time2015:973),2016,2031,1/12,1,0); % 2000 to 2030 in steps of 1 month
    death_vals01_scaled(i)=calc_qalyv2(Y1test_scaled(time2015:973),2016,2031,1/12,1,0); % 2000 to 2030 in steps of 1 month
end
tot_liver_deaths_trt_scaled=sum(death_vals00_scaled.*death_vec')+sum(death_vals10_scaled.*death_vec')+...
    sum(death_vals01_scaled.*death_vec')+sum(death_vals11_scaled.*death_vec');

%figure;plot(2015:2029,Xdeaths,'g-',2015:2029,XTdeaths,'r-','LineWidth',2);
%grid()
figure;plot(2016:2030,100*((XTdeaths_scaled/XTdeaths_scaled(1))),'r-','LineWidth',2);
grid()
title('% deaths change for treating advanced disease scaled')


chronic_prev_current = [chronic_nums10]; % with no treatments chronic_nums11 all zero
total_current = [361:540];
figure;plot(t1,100*sum(XT(chronic_prev_current,:))./sum(XT(total_current,:)),'r--','LineWidth',2)
set(gca,'ylim',[0 70])
xticks(0:5:81)
set(gca,'xticklabel',1950:5:2031,'xticklabelrotation',90)
title('% of current PWID with HCV from 1950 to 2030')
set(gca,'xlim',[0 81])
grid()

%% now run with p = 0

ntreat=5525; % Number of treatments per year, 2015 and onwards
pcom=0.892; % prob of current PWID completing treatment
alpha=0.95; % prob of acheiving SVR
p=0;
[t1,x1t] = ode45(@odeeq_v5T_p,0:1/12:81,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix,...
    scale00,scale10,ntreat,pcom,alpha,66,p); % note 2015 marker is 65 years
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
yinci_PWID=just_integrate(T10(time2015:973),2016,2031,1/12)+...
    just_integrate(T00(time2015:973),2016,2031,1/12)+...
    just_integrate(T01(time2015:973),2016,2031,1/12)+...
    just_integrate(T11(time2015:973),2016,2031,1/12);
figure;plot(2016:2030,100*((yinci_PWID/yinci_PWID(1))),'r-','LineWidth',2);
grid()
title('% incidence change for treating current PWID')

% deaths
% deaths after treatment
% test for 12,13,14 and 15
death_vals10=zeros(4,1);
death_vals00=zeros(4,1);
death_vals11=zeros(4,1);
death_vals01=zeros(4,1);
XTdeaths=zeros(1,2030-2016+1);
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
        just_integrate(Xtest(time2015:973),2016,2031,1/12)*death_vec(i)+...
        just_integrate(Ytest(time2015:973),2016,2031,1/12)*death_vec(i)+...
        just_integrate(X1test(time2015:973),2016,2031,1/12)*death_vec(i)+...
        just_integrate(Y1test(time2015:973),2016,2031,1/12)*death_vec(i);
    %        XTdeaths=XTdeaths+...
    %        just_integrate(Xtest(time2015:961),2015,2030,1/12)*death_vec(i)+...
    %        just_integrate(Ytest(time2015:961),2015,2030,1/12)*death_vec(i);
    death_vals10(i)=calc_qalyv2(Xtest(time2015:961),2016,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
    death_vals00(i)=calc_qalyv2(Ytest(time2015:961),2016,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
    death_vals11(i)=calc_qalyv2(X1test(time2015:961),2016,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
    death_vals01(i)=calc_qalyv2(Y1test(time2015:961),2016,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
end
tot_liver_deaths_trt=sum(death_vals00.*death_vec')+sum(death_vals10.*death_vec')+...
    sum(death_vals01.*death_vec')+sum(death_vals11.*death_vec');
%figure;plot(2015:2029,Xdeaths,'g-',2015:2029,XTdeaths,'r-','LineWidth',2);
%grid()

figure;plot(2016:2030,100*((XTdeaths/XTdeaths(1))),'r-','LineWidth',2);
grid()
title('% deaths change for treating current PWID')


%% run with mixture
n0 = 6906; % was 4906
n1 = 6500; % was 6000
y1=71; % above is applied for 6 years after 2015 (year 66)

pcom=0.892; % prob of current PWID completing treatment
alpha=0.95; % prob of acheiving SVR
ploty = 1; % set to 1 to plot the solution
rate=parm_current(2);
delta=parm_current(1);
k = 0;
%global outyt
%global outyF4
%outyt=[];
%outyF4=[];
for ii = n0 % was 4906
    for j = n1 % was 6000
        n0=ii;
        n1=j;
        k=k+1;
         [t1,x1t] = ode45(@odeeq_v5T_2p,0:1/12:81,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix,...
             scale00,scale10,n0,n1,pcom,alpha,66,y1); % note 2015 marker is 65 years
%               [t1,x1t] = ode45(@odeeq_v5T_2p_scaled,0:1/12:81,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix,...
%             age_scale00,age_scale10,n0,n1,pcom,alpha,66,y1,age_weights',age_weights_former'); % note 2015 marker is 65 years
        XT=x1t';
        %save('testy','outyt','outyF4')
        T00=(1-delta)*rate*scale00*sum(XT(A_comps_00,:));
        T01=(1-delta)*rate*scale00*sum(XT(A_comps_01,:));
        T11=(1-delta)*rate*scale10*sum(XT(A_comps_11,:));
        T10=(1-delta)*rate*scale10*sum(XT(A_comps_10,:));
        % integrate
        yinci_PWID=just_integrate(T10(time2015:973),2016,2031,1/12)+...
            just_integrate(T00(time2015:973),2016,2031,1/12)+...
            just_integrate(T01(time2015:973),2016,2031,1/12)+...
            just_integrate(T11(time2015:973),2016,2031,1/12);
        if ploty==1
            figure;plot(2016:2030,100*((yinci_PWID/yinci_PWID(1))),'r-','LineWidth',2);
            grid()
            title('% incidence change for two treatment regimes')
        end
        
        % deaths
        % deaths after treatment
        % test for 12,13,14 and 15
        m=size(XT);
        Deathmat=zeros(9,m(2));
        death_vals10=zeros(4,1);
        death_vals00=zeros(4,1);
        death_vals11=zeros(4,1);
        death_vals01=zeros(4,1);
        XTdeaths=zeros(1,2030-2016+1);
        % Now can have deaths in the failed treatment compartments
        testcompvec01_failed=[192:20:360;193:20:360;194:20:360;195:20:360];
        testcompvec11_failed=[552:20:720;553:20:720;554:20:720;555:20:720];
        % death_vec= order is %DC,HCC,LT1,LT2 these are the probablities per year
        for i = 1:4
            Deathmat=Deathmat +(scale10*XT(testcompvec10(i,:),:)+scale00*XT(testcompvec00(i,:),:)+...
                scale10*XT(testcompvec11_failed(i,:),:)+scale00*XT(testcompvec01_failed(i,:),:))*death_vec(i);
            Xtest=scale10*sum(XT(testcompvec10(i,:),:));
            Ytest=scale00*sum(XT(testcompvec00(i,:),:));
            X1test=scale10*sum(XT(testcompvec11_failed(i,:),:));
            Y1test=scale00*sum(XT(testcompvec01_failed(i,:),:));
            XTdeaths=XTdeaths+...
                just_integrate(Xtest(time2015:973),2016,2031,1/12)*death_vec(i)+...
                just_integrate(Ytest(time2015:973),2016,2031,1/12)*death_vec(i)+...
                just_integrate(X1test(time2015:973),2016,2031,1/12)*death_vec(i)+...
                just_integrate(Y1test(time2015:973),2016,2031,1/12)*death_vec(i);
            %        XTdeaths=XTdeaths+...
            %        just_integrate(Xtest(time2015:961),2015,2030,1/12)*death_vec(i)+...
            %        just_integrate(Ytest(time2015:961),2015,2030,1/12)*death_vec(i);
            death_vals10(i)=calc_qalyv2(Xtest(time2015:961),2016,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
            death_vals00(i)=calc_qalyv2(Ytest(time2015:961),2016,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
            death_vals11(i)=calc_qalyv2(X1test(time2015:961),2016,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
            death_vals01(i)=calc_qalyv2(Y1test(time2015:961),2016,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
        end
        Ttot_liver_deaths=sum(death_vals00.*death_vec')+sum(death_vals10.*death_vec')+...
            sum(death_vals01.*death_vec')+sum(death_vals11.*death_vec');
        
        % death qalys
        t=1950:1/12:2031;
        D=zeros(9,2031-2016+1);
        for ii = 1 : 9
            k=0;
            for jj = 2016:2031
                jf=find(t==jj):find(t==(jj+1));
                x=t(jf);
                y=Deathmat(ii,jf);
                k =k + 1;
                D(ii,k) = trapz(x,y);
            end
        end
        Q_Ttot_liver_deaths=qalydeaths(D,3,2016,2030); % fatal qualys
        
        % scaled liver deaths and qalys
        Deathmat_scaled=zeros(9,m(2));
        death_vals10_scaled=zeros(4,1);
        death_vals00_scaled=zeros(4,1);
        death_vals11_scaled=zeros(4,1);
        death_vals01_scaled=zeros(4,1);
        XTdeaths_scaled=zeros(1,2030-2016+1);
        % Now can have deaths in the failed treatment compartments
        testcompvec01_failed=[192:20:360;193:20:360;194:20:360;195:20:360];
        testcompvec11_failed=[552:20:720;553:20:720;554:20:720;555:20:720];
        % death_vec= order is %DC,HCC,LT1,LT2 these are the probablities per year
        for i = 1:4
            Deathmat_scaled=Deathmat_scaled +(age_scale10*Apwid.*XT(testcompvec10(i,:),:)+age_scale00*Aformer.*XT(testcompvec00(i,:),:)+...
                age_scale10*Apwid.*XT(testcompvec11_failed(i,:),:)+age_scale00*Aformer.*XT(testcompvec01_failed(i,:),:))*death_vec(i);
            Xtest_scaled=age_scale10*sum(Apwid.*XT(testcompvec10(i,:),:));
            Ytest_scaled=age_scale00*sum(Aformer.*XT(testcompvec00(i,:),:));
            X1test_scaled=age_scale10*sum(Apwid.*XT(testcompvec11_failed(i,:),:));
            Y1test_scaled=age_scale00*sum(Aformer.*XT(testcompvec01_failed(i,:),:));
            XTdeaths_scaled=XTdeaths_scaled+...
                just_integrate(Xtest_scaled(time2015:973),2016,2031,1/12)*death_vec(i)+...
                just_integrate(Ytest_scaled(time2015:973),2016,2031,1/12)*death_vec(i)+...
                just_integrate(X1test_scaled(time2015:973),2016,2031,1/12)*death_vec(i)+...
                just_integrate(Y1test_scaled(time2015:973),2016,2031,1/12)*death_vec(i);
            %        XTdeaths=XTdeaths+...
            %        just_integrate(Xtest(time2015:961),2015,2030,1/12)*death_vec(i)+...
            %        just_integrate(Ytest(time2015:961),2015,2030,1/12)*death_vec(i);
            death_vals10_scaled(i)=calc_qalyv2(Xtest_scaled(time2015:961),2016,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
            death_vals00_scaled(i)=calc_qalyv2(Ytest_scaled(time2015:961),2016,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
            death_vals11_scaled(i)=calc_qalyv2(X1test_scaled(time2015:961),2016,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
            death_vals01_scaled(i)=calc_qalyv2(Y1test_scaled(time2015:961),2016,2030,1/12,1,0); % 2000 to 2030 in steps of 1 month
        end
        Ttot_liver_deaths_scaled=sum(death_vals00_scaled.*death_vec')+sum(death_vals10_scaled.*death_vec')+...
            sum(death_vals01_scaled.*death_vec')+sum(death_vals11_scaled.*death_vec');
        
        % death qalys
        t=1950:1/12:2031;
        D_scaled=zeros(9,2031-2016+1);
        for ii = 1 : 9
            k=0;
            for jj = 2016:2031
                jf=find(t==jj):find(t==(jj+1));
                x=t(jf);
                y=Deathmat_scaled(ii,jf);
                k =k + 1;
                D_scaled(ii,k) = trapz(x,y);
            end
        end
        Q_Ttot_liver_deaths_scaled=qalydeaths(D_scaled,3,2016,2030); % fatal qualys
        
        %figure;plot(2015:2029,Xdeaths,'g-',2015:2029,XTdeaths,'r-','LineWidth',2);
        %grid()
        if ploty==1
            figure;plot(2016:2030,100*((XTdeaths/XTdeaths(1))),'r-','LineWidth',2);
            grid()
            title('% deaths change for two treatment regimes')
        end
        errorsq=(100*((XTdeaths(end)/XTdeaths(1)))-35)^2 + (100*((yinci_PWID(end)/yinci_PWID(1)))-20)^2;
        Xout(k,1) = ii;
        Xout(k,2) = j;
        Xout(k,3) = errorsq;
    end
end


% calculate the qalys for the intervention
% T prefix is for treatment in vector XT
% Now subjects cab be in all compartments
Tcompvec10= [S0_comps_10;S1_comps_10;S2_comps_10;S3_comps_10;S4_comps_10;A_comps_10;...
    F0_comps_10;F1_comps_10;F2_comps_10;F3_comps_10;F4_comps_10;...
    HCC_comps_10;DC_comps_10;LT1_comps_10;LT2_comps_10;...
    T0_comps_10;T1_comps_10;T2_comps_10;T3_comps_10;T4_comps_10];
Tcompvec00= [S0_comps_00;S1_comps_00;S2_comps_00;S3_comps_00;S4_comps_00;A_comps_00;...
    F0_comps_00;F1_comps_00;F2_comps_00;F3_comps_00;F4_comps_00;...
    HCC_comps_00;DC_comps_00;LT1_comps_00;LT2_comps_00;...
    T0_comps_00;T1_comps_00;T2_comps_00;T3_comps_00;T4_comps_00];
Tcompvec11= [S0_comps_11;S1_comps_11;S2_comps_11;S3_comps_11;S4_comps_11;A_comps_11;...
    F0_comps_11;F1_comps_11;F2_comps_11;F3_comps_11;F4_comps_11;...
    HCC_comps_11;DC_comps_11;LT1_comps_11;LT2_comps_11;...
    T0_comps_11;T1_comps_11;T2_comps_11;T3_comps_11;T4_comps_11];
Tcompvec01= [S0_comps_01;S1_comps_01;S2_comps_01;S3_comps_01;S4_comps_01;A_comps_01;...
    F0_comps_01;F1_comps_01;F2_comps_01;F3_comps_01;F4_comps_01;...
    HCC_comps_01;DC_comps_01;LT1_comps_01;LT2_comps_01;...
    T0_comps_01;T1_comps_01;T2_comps_01;T3_comps_01;T4_comps_01];
% order for below is S0-S4,A,F0-F4,(DCC,HCC,LT1,LT2) and T0-T4
Tutility_vec=[0.93 0.93 0.93 0.93 0.93 ...
    0.77 ...
    0.77 0.77 0.77 0.66 0.55 ...
    0.45 0.45 0.45 0.67 ...
    0.77 0.77 0.77 0.66 0.55];

Tqualy_vals00=zeros(20,1); % This conatins qualy for all compartments 00 (former, not failed)
Tqualy_vals01=zeros(20,1); % This conatins qualy for all compartments 01 (former, failed)
Tqualy_vals10=zeros(20,1); % This conatins qualy for all compartments 10 (current, not failed)
Tqualy_vals11=zeros(20,1); % This conatins qualy for all compartments 11 (current, failed)

for i = 1 : 20
    Xtest00=scale00*sum(XT(Tcompvec00(i,:),:));
    Xtest01=scale00*sum(XT(Tcompvec01(i,:),:));
    Xtest10=scale10*sum(XT(Tcompvec10(i,:),:));
    Xtest11=scale10*sum(XT(Tcompvec11(i,:),:));
    Tqualy_vals00(i)=calc_qalyv2(Xtest00(time2015:973),2016,2031,1/12,Tutility_vec(i),3); % 2000 to 2030 in steps of 1 month
    Tqualy_vals01(i)=calc_qalyv2(Xtest01(time2015:973),2016,2031,1/12,Tutility_vec(i),3); % 2000 to 2030 in steps of 1 month
    Tqualy_vals10(i)=calc_qalyv2(Xtest10(time2015:973),2016,2031,1/12,Tutility_vec(i),3); % 2000 to 2030 in steps of 1 month
    Tqualy_vals11(i)=calc_qalyv2(Xtest11(time2015:973),2016,2031,1/12,Tutility_vec(i),3); % 2000 to 2030 in steps of 1 month
end
Ttotal_qual_current = sum(Tqualy_vals00)+sum(Tqualy_vals01)+sum(Tqualy_vals10)+sum(Tqualy_vals11);

% scaled
Tqualy_vals00_scaled=zeros(20,1); % This conatins qualy for all compartments 00 (former, not failed)
Tqualy_vals01_scaled=zeros(20,1); % This conatins qualy for all compartments 01 (former, failed)
Tqualy_vals10_scaled=zeros(20,1); % This conatins qualy for all compartments 10 (current, not failed)
Tqualy_vals11_scaled=zeros(20,1); % This conatins qualy for all compartments 11 (current, failed)

for i = 1 : 20
    Xtest00_scaled=age_scale00*sum(Aformer.*XT(Tcompvec00(i,:),:));
    Xtest01_scaled=age_scale00*sum(Aformer.*XT(Tcompvec01(i,:),:));
    Xtest10_scaled=age_scale10*sum(Apwid.*XT(Tcompvec10(i,:),:));
    Xtest11_scaled=age_scale10*sum(Apwid.*XT(Tcompvec11(i,:),:));
    Tqualy_vals00_scaled(i)=calc_qalyv2(Xtest00_scaled(time2015:973),2016,2031,1/12,Tutility_vec(i),3); % 2000 to 2030 in steps of 1 month
    Tqualy_vals01_scaled(i)=calc_qalyv2(Xtest01_scaled(time2015:973),2016,2031,1/12,Tutility_vec(i),3); % 2000 to 2030 in steps of 1 month
    Tqualy_vals10_scaled(i)=calc_qalyv2(Xtest10_scaled(time2015:973),2016,2031,1/12,Tutility_vec(i),3); % 2000 to 2030 in steps of 1 month
    Tqualy_vals11_scaled(i)=calc_qalyv2(Xtest11_scaled(time2015:973),2016,2031,1/12,Tutility_vec(i),3); % 2000 to 2030 in steps of 1 month
end
Ttotal_qual_current_scaled = sum(Tqualy_vals00_scaled)+sum(Tqualy_vals01_scaled)+sum(Tqualy_vals10_scaled)+sum(Tqualy_vals11_scaled);

% calculate the costs - this is the costs not including treatments

[costanual,costoneoff,namesof_anualcost,namesof_oneoff]=cost_params();
% anual costs
Tcompvec00= [F0_comps_00;F1_comps_00;F2_comps_00;F3_comps_00;F4_comps_00;HCC_comps_00;DC_comps_00];
Tcompvec01= [F0_comps_01;F1_comps_01;F2_comps_01;F3_comps_01;F4_comps_01;HCC_comps_01;DC_comps_01];
Tcompvec10= [F0_comps_10;F1_comps_10;F2_comps_10;F3_comps_10;F4_comps_10;HCC_comps_10;DC_comps_10];
Tcompvec11= [F0_comps_11;F1_comps_11;F2_comps_11;F3_comps_11;F4_comps_11;HCC_comps_11;DC_comps_11];

Tcost_vals00=zeros(7,1); % This conatins qualy for 00
Tcost_vals01=zeros(7,1); % This conatins qualy for 01
Tcost_vals10=zeros(7,1); % This conatins qualy for 10
Tcost_vals11=zeros(7,1); % This conatins qualy for 11

for i = 1 : 7
    Xtest00=scale00*sum(XT(Tcompvec00(i,:),:));
    Xtest01=scale00*sum(XT(Tcompvec01(i,:),:));
    Xtest10=scale10*sum(XT(Tcompvec10(i,:),:));
    Xtest11=scale10*sum(XT(Tcompvec11(i,:),:));
    Tcost_vals00(i)=calc_qalyv2(Xtest00(time2015:973),2016,2031,1/12,costanual(i),3); % 2000 to 2030 in steps of 1 month
    Tcost_vals01(i)=calc_qalyv2(Xtest01(time2015:973),2016,2031,1/12,costanual(i),3); % 2000 to 2030 in steps of 1 month
    Tcost_vals10(i)=calc_qalyv2(Xtest10(time2015:973),2016,2031,1/12,costanual(i),3); % 2000 to 2030 in steps of 1 month
    Tcost_vals11(i)=calc_qalyv2(Xtest11(time2015:973),2016,2031,1/12,costanual(i),3); % 2000 to 2030 in steps of 1 month
end
Ttotal_annual_costs = sum(Tcost_vals00)+sum(Tcost_vals01)+sum(Tcost_vals10)+sum(Tcost_vals11);

% as above for scaled
Tcost_vals00_scaled=zeros(7,1); % This conatins qualy for 00
Tcost_vals01_scaled=zeros(7,1); % This conatins qualy for 01
Tcost_vals10_scaled=zeros(7,1); % This conatins qualy for 10
Tcost_vals11_scaled=zeros(7,1); % This conatins qualy for 11

for i = 1 : 7
    Xtest00_scaled=age_scale00*sum(Aformer.*XT(Tcompvec00(i,:),:));
    Xtest01_scaled=age_scale00*sum(Aformer.*XT(Tcompvec01(i,:),:));
    Xtest10_scaled=age_scale10*sum(Apwid.*XT(Tcompvec10(i,:),:));
    Xtest11_scaled=age_scale10*sum(Apwid.*XT(Tcompvec11(i,:),:));
    Tcost_vals00_scaled(i)=calc_qalyv2(Xtest00_scaled(time2015:973),2016,2031,1/12,costanual(i),3); % 2000 to 2030 in steps of 1 month
    Tcost_vals01_scaled(i)=calc_qalyv2(Xtest01_scaled(time2015:973),2016,2031,1/12,costanual(i),3); % 2000 to 2030 in steps of 1 month
    Tcost_vals10_scaled(i)=calc_qalyv2(Xtest10_scaled(time2015:973),2016,2031,1/12,costanual(i),3); % 2000 to 2030 in steps of 1 month
    Tcost_vals11_scaled(i)=calc_qalyv2(Xtest11_scaled(time2015:973),2016,2031,1/12,costanual(i),3); % 2000 to 2030 in steps of 1 month
end
Ttotal_annual_costs_scaled = sum(Tcost_vals00_scaled)+sum(Tcost_vals01_scaled)+sum(Tcost_vals10_scaled)+sum(Tcost_vals11_scaled);

% total liver costs
% look at number of liver transplants
% DC to LT the rate is parm_names(13)
LTcost=costoneoff(3);
r_DCLT=parm_current(13);
T00=r_DCLT*scale00*sum(XT(DC_comps_00,:));
T01=r_DCLT*scale00*sum(XT(DC_comps_01,:));
T10=r_DCLT*scale10*sum(XT(DC_comps_10,:));
T11=r_DCLT*scale10*sum(XT(DC_comps_11,:));
% integrate
DCLT00=calc_qalyv2(T00(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
DCLT01=calc_qalyv2(T01(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
DCLT10=calc_qalyv2(T10(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
DCLT11=calc_qalyv2(T11(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month

% HCC to LT the rate is parm_names(15)
r_HCCLT=parm_current(15);
T00=r_HCCLT*scale00*sum(XT(HCC_comps_00,:));
T01=r_HCCLT*scale00*sum(XT(HCC_comps_01,:));
T10=r_HCCLT*scale10*sum(XT(HCC_comps_10,:));
T11=r_HCCLT*scale10*sum(XT(HCC_comps_11,:));
% integrate
HCCLT00=calc_qalyv2(T00(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
HCCLT01=calc_qalyv2(T01(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
HCCLT10=calc_qalyv2(T10(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
HCCLT11=calc_qalyv2(T11(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
Ttotal_cost_liver_transplants = DCLT00 + DCLT10 + HCCLT00 + HCCLT10+...
    DCLT01 + DCLT11 + HCCLT01 + HCCLT11;
% as above scaled
T00_scaled=r_DCLT*age_scale00*sum(Aformer.*XT(DC_comps_00,:));
T01_scaled=r_DCLT*age_scale00*sum(Aformer.*XT(DC_comps_01,:));
T10_scaled=r_DCLT*age_scale10*sum(Apwid.*XT(DC_comps_10,:));
T11_scaled=r_DCLT*age_scale10*sum(Apwid.*XT(DC_comps_11,:));
% integrate
DCLT00_scaled=calc_qalyv2(T00_scaled(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
DCLT01_scaled=calc_qalyv2(T01_scaled(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
DCLT10_scaled=calc_qalyv2(T10_scaled(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
DCLT11_scaled=calc_qalyv2(T11_scaled(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month

% HCC to LT the rate is parm_names(15)
T00_scaled=r_HCCLT*age_scale00*sum(Aformer.*XT(HCC_comps_00,:));
T01_scaled=r_HCCLT*age_scale00*sum(Aformer.*XT(HCC_comps_01,:));
T10_scaled=r_HCCLT*age_scale10*sum(Apwid.*XT(HCC_comps_10,:));
T11_scaled=r_HCCLT*age_scale10*sum(Apwid.*XT(HCC_comps_11,:));
% integrate
HCCLT00_scaled=calc_qalyv2(T00_scaled(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
HCCLT01_scaled=calc_qalyv2(T01_scaled(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
HCCLT10_scaled=calc_qalyv2(T10_scaled(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
HCCLT11_scaled=calc_qalyv2(T11_scaled(time2015:973),2016,2031,1/12,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
Ttotal_cost_liver_transplants_scaled = DCLT00_scaled + DCLT10_scaled + HCCLT00_scaled + HCCLT10_scaled+...
    DCLT01_scaled + DCLT11_scaled + HCCLT01_scaled + HCCLT11_scaled;

% one off costs for diagnosis of incident chronic cases
rate=parm_current(2);
delta=parm_current(1);
T00=(1-delta)*rate*scale00*sum(XT(A_comps_00,:));
T01=(1-delta)*rate*scale00*sum(XT(A_comps_01,:));
T10=(1-delta)*rate*scale10*sum(XT(A_comps_10,:));
T11=(1-delta)*rate*scale10*sum(XT(A_comps_11,:));
% integrate
inci_acute00=calc_qalyv2(T00(time2015:973),2016,2031,1/12,costoneoff(4),3); % 2000 to 2030 in steps of 1 month
inci_acute01=calc_qalyv2(T01(time2015:973),2016,2031,1/12,costoneoff(4),3); % 2000 to 2030 in steps of 1 month
inci_acute10=calc_qalyv2(T10(time2015:973),2016,2031,1/12,costoneoff(4),3); % 2000 to 2030 in steps of 1 month
inci_acute11=calc_qalyv2(T11(time2015:973),2016,2031,1/12,costoneoff(4),3); % 2000 to 2030 in steps of 1 month
Tyearly_diagnosed_cases_in_chronic_stages=inci_acute00+inci_acute01+inci_acute10+inci_acute11;

% as above scaled
T00_scaled=(1-delta)*rate*age_scale00*sum(Aformer.*XT(A_comps_00,:));
T01_scaled=(1-delta)*rate*age_scale00*sum(Aformer.*XT(A_comps_01,:));
T10_scaled=(1-delta)*rate*age_scale10*sum(Apwid.*XT(A_comps_10,:));
T11_scaled=(1-delta)*rate*age_scale10*sum(Apwid.*XT(A_comps_11,:));
% integrate
inci_acute00_scaled=calc_qalyv2(T00_scaled(time2015:973),2016,2031,1/12,costoneoff(4),3); % 2000 to 2030 in steps of 1 month
inci_acute01_scaled=calc_qalyv2(T01_scaled(time2015:973),2016,2031,1/12,costoneoff(4),3); % 2000 to 2030 in steps of 1 month
inci_acute10_scaled=calc_qalyv2(T10_scaled(time2015:973),2016,2031,1/12,costoneoff(4),3); % 2000 to 2030 in steps of 1 month
inci_acute11_scaled=calc_qalyv2(T11_scaled(time2015:973),2016,2031,1/12,costoneoff(4),3); % 2000 to 2030 in steps of 1 month
Tyearly_diagnosed_cases_in_chronic_stages_scaled=inci_acute00_scaled+inci_acute01_scaled+inci_acute10_scaled+inci_acute11_scaled;

% one off costs for transition incidence
% 'r_F3F4','r_F4DC','r_F4HCC'
% F3 to F4
r_F3F4=parm_current(9);
T00=r_F3F4*scale00*sum(XT(F3_comps_00,:));
T01=r_F3F4*scale00*sum(XT(F3_comps_01,:));
T10=r_F3F4*scale10*sum(XT(F3_comps_10,:));
T11=r_F3F4*scale10*sum(XT(F3_comps_11,:));
% integrate
F3F400=calc_qalyv2(T00(time2015:973),2016,2031,1/12,costoneoff(1),3); % 2000 to 2030 in steps of 1 month
F3F401=calc_qalyv2(T01(time2015:973),2016,2031,1/12,costoneoff(1),3); % 2000 to 2030 in steps of 1 month
F3F410=calc_qalyv2(T10(time2015:973),2016,2031,1/12,costoneoff(1),3); % 2000 to 2030 in steps of 1 month
F3F411=calc_qalyv2(T11(time2015:973),2016,2031,1/12,costoneoff(1),3); % 2000 to 2030 in steps of 1 month

% F4 to HCC
r_F4HCC=parm_current(11);
T00=r_F4HCC*scale00*sum(XT(F4_comps_00,:));
T01=r_F4HCC*scale00*sum(XT(F4_comps_01,:));
T10=r_F4HCC*scale10*sum(XT(F4_comps_10,:));
T11=r_F4HCC*scale10*sum(XT(F4_comps_11,:));
% integrate
F4HCC00=calc_qalyv2(T00(time2015:973),2016,2031,1/12,costoneoff(2),3); % 2000 to 2030 in steps of 1 month
F4HCC01=calc_qalyv2(T01(time2015:973),2016,2031,1/12,costoneoff(2),3); % 2000 to 2030 in steps of 1 month
F4HCC10=calc_qalyv2(T10(time2015:973),2016,2031,1/12,costoneoff(2),3); % 2000 to 2030 in steps of 1 month
F4HCC11=calc_qalyv2(T11(time2015:973),2016,2031,1/12,costoneoff(2),3); % 2000 to 2030 in steps of 1 month


Ttotal_transistion_costs = F3F400 + F3F401 + F3F410 + F3F411 +...
    F4HCC00 + F4HCC01 + F4HCC10 + F4HCC11 ;

% as above scaled
T00_scaled=r_F3F4*age_scale00*sum(Aformer.*XT(F3_comps_00,:));
T01_scaled=r_F3F4*age_scale00*sum(Aformer.*XT(F3_comps_01,:));
T10_scaled=r_F3F4*age_scale10*sum(Apwid.*XT(F3_comps_10,:));
T11_scaled=r_F3F4*age_scale10*sum(Apwid.*XT(F3_comps_11,:));
% integrate
F3F400_scaled=calc_qalyv2(T00_scaled(time2015:973),2016,2031,1/12,costoneoff(1),3); % 2000 to 2030 in steps of 1 month
F3F401_scaled=calc_qalyv2(T01_scaled(time2015:973),2016,2031,1/12,costoneoff(1),3); % 2000 to 2030 in steps of 1 month
F3F410_scaled=calc_qalyv2(T10_scaled(time2015:973),2016,2031,1/12,costoneoff(1),3); % 2000 to 2030 in steps of 1 month
F3F411_scaled=calc_qalyv2(T11_scaled(time2015:973),2016,2031,1/12,costoneoff(1),3); % 2000 to 2030 in steps of 1 month

% F4 to HCC
T00_scaled=r_F4HCC*age_scale00*sum(Aformer.*XT(F4_comps_00,:));
T01_scaled=r_F4HCC*age_scale00*sum(Aformer.*XT(F4_comps_01,:));
T10_scaled=r_F4HCC*age_scale10*sum(Apwid.*XT(F4_comps_10,:));
T11_scaled=r_F4HCC*age_scale10*sum(Apwid.*XT(F4_comps_11,:));
% integrate
F4HCC00_scaled=calc_qalyv2(T00_scaled(time2015:973),2016,2031,1/12,costoneoff(2),3); % 2000 to 2030 in steps of 1 month
F4HCC01_scaled=calc_qalyv2(T01_scaled(time2015:973),2016,2031,1/12,costoneoff(2),3); % 2000 to 2030 in steps of 1 month
F4HCC10_scaled=calc_qalyv2(T10_scaled(time2015:973),2016,2031,1/12,costoneoff(2),3); % 2000 to 2030 in steps of 1 month
F4HCC11_scaled=calc_qalyv2(T11_scaled(time2015:973),2016,2031,1/12,costoneoff(2),3); % 2000 to 2030 in steps of 1 month
Ttotal_transistion_costs_scaled = F3F400_scaled + F3F401_scaled + F3F410_scaled + F3F411_scaled +...
    F4HCC00_scaled + F4HCC01_scaled + F4HCC10_scaled + F4HCC11_scaled ;


% costs of follow-up on F4
%
time2015_inODEunits=66;
moved_to_T4_former=zeros(9,length(t1));
moved_to_T4_current=zeros(9,length(t1));
total_not_T4=zeros(9,length(t1));
total_T4=zeros(9,length(t1));
for i = 1 : length(t1)
    N=XT(:,i);
    [~,~,moved_to_T4_former(:,i),moved_to_T4_current(:,i),total_not_T4(:,i),total_T4(:,i)]=treat_comps_pv3(n0,n1,y1,N,scale00,scale10,t1(i),pcom,alpha,time2015_inODEunits);
end
%outyF4z=outyF4(:,3:4:4*length(outyt)); % former mort_current,mort_former
%outyF4z2=outyF4(:,4:4:4*length(outyt)); % former mort_current,mort_former

totaltestsF4=calc_ongoingF4treats(t1,total_T4,2016,2031);
totaltests_notF4=calc_ongoingF4treats(t1,total_not_T4,2016,2031);

figure;plot(2016:2030,totaltestsF4)

% outyF4 in fours
% 1. F4former inci, 2. F4current, 3. inci total not F4, 4. total F4
%outyF4z=outyF4(:,1:4:4*length(outyt)); % former mort_current,mort_former
pyformer=calc_ongoingF4costs(t1,moved_to_T4_former,2016,2031,mort_former,3);


pycurrent=calc_ongoingF4costs(t1,moved_to_T4_current,2016,2031,mort_current,3);

Tongoingf4costs = 557.62*(pyformer+pycurrent);

% total costs of treatments
% 6906 for 15 years
% 7500 for 5 years
% cost of 35,806.67 for non F4
% cost of 62,377.15 for F4
discount_rate=0.03;
cost_discount = 1./((1+discount_rate).^(0:(length(totaltests_notF4)-1)));
Ttreatcost = sum((totaltests_notF4*35806.67 + totaltestsF4 * 62377.15).*cost_discount);
% note slight over estimate due to trap rule and discontuity error

Ttotal_costs = Ttotal_annual_costs+Ttotal_cost_liver_transplants+ ...
    Tyearly_diagnosed_cases_in_chronic_stages+Ttotal_transistion_costs+...
    Tongoingf4costs+Ttreatcost;

Ttotal_costs_scaled = Ttotal_annual_costs_scaled+Ttotal_cost_liver_transplants_scaled+ ...
    Tyearly_diagnosed_cases_in_chronic_stages_scaled+Ttotal_transistion_costs_scaled+...
    Tongoingf4costs+Ttreatcost;

Ttotal_qual_current=Ttotal_qual_current+Q_Ttot_liver_deaths;
Ttotal_qual_current_scaled=Ttotal_qual_current_scaled+Q_Ttot_liver_deaths_scaled;
ICER = (Ttotal_costs-total_costs)/(Ttotal_qual_current-total_qual_current);
ICER_scaled = (Ttotal_costs_scaled-total_costs_scaled)/(Ttotal_qual_current_scaled-total_qual_current_scaled);
outy{2,1}='Composite treatment';
outy{2,2}=['$' char(v.format(round(Ttotal_costs)))];
outy{2,3}=char(v.format(round(Ttotal_qual_current)));
outy{2,4}=char(v.format(round(Ttot_liver_deaths)));
outy{2,5}=['$' char(v.format(round(ICER)))];

cell2table(outy,'VariableNames',{'Intervention' 'TotalCost','TotalQALYS' 'LiverRelatedDeaths' 'ICER'})








