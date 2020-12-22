function fitode_v4_abc(Nsamp,filenum)
% 30 August 2020
% 31 August 2020 chronic HCV prevlance calibrated at 50% for current PWID in 2015.
% This creates the data for the abc to move to R
% summary statistics
% gradient at 2002 for % of current PWID - should be zero
% value at this point - should be 60%
% value at 2015 for % of current PWID - should be 50%
% current & former F0 and F1 should be 66% in 2015 of all chronic PWID


% delete figures
%************************************
% **** note M matrices 15,15 should be minus - it is a LT2 death
% also set the svr4 to dc and hcc to zero, not correctly implemented
delete(get(0,'children'))
% get the data
% solve the diff equation
t = 0:80; % 1950 to 2030
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
chronic_prev_current = [chronic_nums10]; % with no treatments chronic_nums11 all zero
total_current = [361:540];

% params to change
R=zeros(Nsamp,5);
S=zeros(Nsamp,3);
for i = 1 : Nsamp
    i
[R(i,:),N0,M,Mdash,mort_current,mort_former,extra_parms_vals, phiminusvals, phiplussvals,extra_parms_nams,age_matrix,parm_names,parm_current,parm_former]=...
    create_params_v4_with_change;
%[t1,x1] = ode15s(@odeeq_v4,t,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix);
[t1,x1] = ode45(@odeeq_v4,t,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix);
X=x1';    
% chronic HCV prevlance calibrated at 50% for current PWID in 2015.
% chronic_prev_current = [chronic_nums10]; % with no treatments chronic_nums11 all zero
% total_current = [361:540];
% figure;plot(t1,100*sum(X(chronic_prev_current,:))./sum(X(total_current,:)),'r--','LineWidth',2)
% set(gca,'ylim',[0 70])
% xticks(0:5:80)
% set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
% grid()

scale00=144000/sum(X([chronic_nums00],66)); % former
scale10=40000/sum(X([chronic_nums10],66)); % current
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

f=sum(X(chronic_prev_current,:))./sum(X(total_current,:));
sum2=f(53); % target 0.6
sum3=sum(X(chronic_prev_current,66))/sum(X(total_current,66)); % target 0.5
sum4= (f0(66)+f1(66))/(f0(66)+f1(66)+f2(66)+f3(66)+f4(66)+hcc(66)+dc(66)+lt(66)); % target = 0.66
S(i,:)=[sum2,sum3,sum4];
end
T=[0.6 0.5 0.51];
save(['F:\HCV\version4\Smatv2_abc',num2str(filenum),'.mat'],'S') % summary stats (4)
save(['F:\HCV\version4\Rmatv2_abc',num2str(filenum),'.mat'],'R') % parameters (6)
save(['F:\HCV\version4\Tmatv2_abc',num2str(filenum),'.mat'],'T') % target (4)

% figure;plot(t1,scale00*sum(X([chronic_nums00],:))+scale10*sum(X([chronic_nums10],:)))
% xticks(0:5:80)
% set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
% grid()




% figure;area(t1,[lt' hcc' dc' f4' f3' f2' f1' f0'])
% xticks(0:5:80)
% set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
% grid()

% figure;plot(t1,sum(X,1))
% set(gca,'ylim',[900 1100])
% **************** This is not constant - should always be 1,000
% check for leakage - there is the sum is not constant
% D01 and D11 should always be zero as there are treatments to fail
% these should all be zero as they have no entries into compartments  
% These are OK
