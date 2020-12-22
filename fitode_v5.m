function fitode_v4
% 30 August 2020
% 31 August 2020 chronic HCV prevlance calibrated at 50% for current PWID in 2015.
% delete figures
%************************************
% **** note M matrices 15,15 should be minus - it is a LT2 death
% also set the svr4 to dc and hcc to zero, not correctly implemented
delete(get(0,'children'))
% get the data
[N0,M,Mdash,mort_current,mort_former,extra_parms_vals, phiminusvals, phiplussvals,extra_parms_nams,age_matrix,parm_names,parm_current,parm_former]= create_params_v4;
% solve the diff equation
t = 1 : 81; % 1950 to 2030
%[t1,x1] = ode15s(@odeeq_v4,0:80,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix);
[t1,x1] = ode45(@odeeq_v4,0:80,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix);
[compnames,chronic_nams,...
 chronic_nums00,chronic_nums01,chronic_nums10,chronic_nums11,...
 F0_comps,F1_comps,F2_comps,F3_comps,F4_comps,DC_comps,HCC_comps,LT1_comps,LT2_comps,...
 S0_comps,S1_comps,S2_comps,S3_comps,S4_comps,A_comps,...
 T0_comps,T1_comps,T2_comps,T3_comps,T4_comps,...
 LT1_comps_00,LT1_comps_10,...
 LT2_comps_00,LT2_comps_10,...
 DC_comps_00,DC_comps_10,...
 HCC_comps_00,HCC_comps_10,...
 F0_comps_00,F0_comps_10,...
 F1_comps_00,F1_comps_10,...
 F2_comps_00,F2_comps_10,...
 F3_comps_00,F3_comps_10,...
 F4_comps_00,F4_comps_10]=coeffsof_v4;
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

scale00=144000/sum(X([chronic_nums00],66)); % former
scale10=40000/sum(X([chronic_nums10],66)); % current
figure;plot(t1,scale00*sum(X([chronic_nums00],:))+scale10*sum(X([chronic_nums10],:)))
xticks(0:5:80)
set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
grid()
title('Absolute numbers of PWID with chronic HCV, scaled to the 2015 estimates')

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
% **************** run to here
figure;plot(t1,sum(X,1))
set(gca,'ylim',[900 1100])
title('Should be constant - tests for leaks')
% **************** This is not constant - should always be 1,000
% check for leakage - there is the sum is not constant
% D01 and D11 should always be zero as there are treatments to fail
% these should all be zero as they have no entries into compartments  
% These are OK
D01test = sum(sum(X(181:360,:))==0);
D11test = sum(sum(X(541:720,:))==0);
% all ages other than first age group should be zero
% These are OK too
D00agetest = sum(sum(X(21:180,:))==0);
D10agetest = sum(sum(X(381:540,:))==0);
% So there is a problem in the first 20 rows of D00 and D10
% comps 2,3,4,5,16,17,18,19,20 should be zero s1-s4 and t1 
% All OK
D00zerocomptest = sum(sum(X([2,3,4,5,16,17,18,19,20],:))==0);
D10zerocomptest = sum(sum(X((380+[2,3,4,5,16,17,18,19,20]),:))==0);
% so problem is in S0,A,F0,F1,F2,F3,F4,DC,HCC,LT1,LT2


figure;
plot(t1,sum(X),'r-')
total_chron00=sum(X(chronic_nums00,:));
total_chron01=sum(X(chronic_nums01,:));
total_chron10=sum(X(chronic_nums10,:));
total_chron11=sum(X(chronic_nums11,:));

totY = total_chron00+total_chron01+total_chron10+total_chron11 + sum(X(A_comps,:));
figure;
plot(t1,totY,'r-')

total_S0=sum(X(S0_comps,:));
total_S1=sum(X(S1_comps,:));
total_S2=sum(X(S2_comps,:));
total_S3=sum(X(S3_comps,:));
total_S4=sum(X(S4_comps,:));
totS = total_S0+total_S1+total_S2+total_S3+total_S4;
figure;
plot(t1,totS,'r-')