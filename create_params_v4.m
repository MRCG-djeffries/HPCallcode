function [N0,M,Mdash,mort_current,mort_former,extra_parms_vals, phiminusvals, phiplussvals,extra_parms_nams,age_matrix,parm_names,parm_current,parm_former]= create_params_v4

% 29 August 2020 - this is the test version
% 30 Ausgust 2020 - version 3

% sets death rates to zero
% sets age matrix to zero

% This creates two matrices for the static parameters
% M is for former i=0
% Mdash is for current i=1, i.e. PWID
parm_names = {'delta','r_AF0','w','r_SVR4DC','r_SVR4HCC','r_F0F1','r_F1F2','r_F2F3','r_F3F4','r_F4DC','r_F4HCC',...
               'r_DCHCC','r_DCLT','r_DCdeath','r_HCCLT','r_HCCdeath',...
               'r_LT1death','r_LT2death'};
 %prob = 1-exp(-rate)
 %rate = -log(1-prob)
 parm_current(1)=0.26; parm_former(1) = parm_current(1); % spontaneous clearance
 parm_current(2)=52/12; parm_former(2) = parm_current(2);% acute duration is 12 weeks
 %parm_current(3)=52/16.56; parm_former(3) = parm_current(3);% duration to return to susceptible after treated
 % for no treatment
 parm_current(3)=0; parm_former(3) = parm_current(3);% duration to return to susceptible after treated
%  parm_current(4)=-log(1-0.02); parm_former(4) = parm_current(4);% SVR4 to DC
%  parm_current(5)=-log(1-0.02); parm_former(5) = parm_current(5);% SVR4 to HCC
%  
 parm_current(4)=0; parm_former(4) = parm_current(4);% SVR4 to DC
 parm_current(5)=0; parm_former(5) = parm_current(5);% SVR4 to HCC
 
 
 parm_current(6) = -log(1-0.116); parm_former(6) = -log(1-0.106); % F0 to F1 current then former
 parm_current(7) = -log(1-0.085); parm_former(7) = -log(1-0.074); % F1 to F2 current then former
 parm_current(8) = -log(1-0.085); parm_former(8) = -log(1-0.106); % F2 to F3 current then former
 parm_current(9) = -log(1-0.130); parm_former(9) = -log(1-0.105); % F3 to F4 current then former
 parm_current(10) = -log(1-0.037); parm_former(10) = parm_current(10); % F4 to DC current then former
 parm_current(11) = -log(1-0.01); parm_former(11) = parm_current(11); % F4 to HCC current then former
 parm_current(12) = -log(1-0.068); parm_former(12) = parm_current(12); % DC to HCC
 parm_current(13) = -log(1-0.033); parm_former(13) = parm_current(13); % DC to LT
 parm_current(15) = -log(1-0.1); parm_former(15) = parm_current(15); % HCC to LT
 
 
%  parm_current(2) = 0.25*parm_current(2);
%  parm_former(2) = 0.25*parm_former(2);
%  parm_current(6) = 0.05*parm_current(6);
%  parm_former(6) = 0.05*parm_former(6);
%  parm_current(7) = 0.5*parm_current(7);
%  parm_former(7) = 0.5*parm_former(7);
 
 
 
 parm_current(14) = -log(1-0.138) ; parm_former(14) = parm_current(14); % DC to death 
 parm_current(16) = -log(1-0.605); parm_former(16) = parm_current(16); % HCC death
 parm_current(17) = -log(1-0.169); parm_former(17) = parm_current(17); % LT to death year 1
 parm_current(18) = -log(1-0.034); parm_former(18) = parm_current(18); % LT to death year 2
 % ***********************************************
%  parm_current([14 16 17 18]) = 0;
%  parm_current([14 16 17 18]) = 0;
%  parm_former([14 16 17 18]) = 0;
%  parm_former([14 16 17 18]) = 0;

 phi = 0;% place holder
 
 Mdash=zeros(20,20);
 phiminusvals=sub2ind(size(Mdash),[1 2 3 4 5],[1 2 3 4 5]);
 phiplussvals=sub2ind(size(Mdash),[6 8 9 10 11],[1 2 3 4 5]);
 Mdash(1,1)=-phi;Mdash(1,6)=parm_current(1)*parm_current(2);Mdash(1,16)=parm_current(3);
 Mdash(2,2)=-phi;Mdash(2,17)=parm_current(3);
 Mdash(3,3)=-phi;Mdash(3,18)=parm_current(3);
 Mdash(4,4)=-phi;Mdash(4,19)=parm_current(3);
 Mdash(5,5)=-phi-parm_current(4)-parm_current(5);Mdash(5,20)=parm_current(3);
 Mdash(6,1)=phi;Mdash(6,6)=-parm_current(2);
 Mdash(7,6)=(1-parm_current(1))*parm_current(2);Mdash(7,7)=-parm_current(6);
 Mdash(8,2)=phi;Mdash(8,7)=parm_current(6);Mdash(8,8)=-parm_current(7);
 Mdash(9,3)=phi;Mdash(9,8)=parm_current(7);Mdash(9,9)=-parm_current(8);
 Mdash(10,4)=phi;Mdash(10,9)=parm_current(8);Mdash(10,10)=-parm_current(9);
 Mdash(11,5)=phi;Mdash(11,10)=parm_current(9);Mdash(11,11)=-parm_current(10)-parm_current(11);
 Mdash(12,11)=parm_current(10)+parm_current(4);Mdash(12,12) = -parm_current(12)-parm_current(13)-parm_current(14);
 Mdash(13,11)=parm_current(11)+parm_current(5);Mdash(13,12)=parm_current(12);Mdash(13,13)=-parm_current(15)-parm_current(16);
 Mdash(14,12)=parm_current(13);Mdash(14,13)=parm_current(15);Mdash(14,14)=-1-parm_current(17);
 Mdash(15,14)=1;Mdash(15,15)=-parm_current(18);
 Mdash(16,16)=-parm_current(3);
 Mdash(17,17)=-parm_current(3);
 Mdash(18,18)=-parm_current(3);
 Mdash(19,19)=-parm_current(3);
 Mdash(20,20)=-parm_current(3);

 
 M=zeros(20,20);
 M(1,1)=-phi;M(1,6)=parm_former(1)*parm_former(2);M(1,16)=parm_former(3);
 M(2,2)=-phi;M(2,17)=parm_former(3);
 M(3,3)=-phi;M(3,18)=parm_former(3);
 M(4,4)=-phi;M(4,19)=parm_former(3);
 M(5,5)=-phi-parm_former(4)-parm_former(5);M(5,20)=parm_former(3);
 M(6,1)=phi;M(6,6)=-parm_former(2);
 M(7,6)=(1-parm_former(1))*parm_former(2);M(7,7)=-parm_former(6);
 M(8,2)=phi;M(8,7)=parm_former(6);M(8,8)=-parm_former(7);
 M(9,3)=phi;M(9,8)=parm_former(7);M(9,9)=-parm_former(8);
 M(10,4)=phi;M(10,9)=parm_former(8);M(10,10)=-parm_former(9);
 M(11,5)=phi;M(11,10)=parm_former(9);M(11,11)=-parm_former(10)-parm_former(11);
 M(12,11)=parm_former(10)+parm_former(4);M(12,12) = -parm_former(12)-parm_former(13)-parm_former(14);
 M(13,11)=parm_former(11)+parm_former(5);M(13,12)=parm_former(12);M(13,13)=-parm_former(15)-parm_former(16);
 M(14,12)=parm_former(13);M(14,13)=parm_former(15);M(14,14)=-1-parm_former(17);
 M(15,14)=1;M(15,15)=-parm_former(18);
 M(16,16)=-parm_former(3);
 M(17,17)=-parm_former(3);
 M(18,18)=-parm_former(3);
 M(19,19)=-parm_former(3);
 M(20,20)=-parm_former(3);
 
 age_matrix=zeros(9,9);
 age_matrix(1,1) = -1/5;
 age_matrix(2,1) = 1/5;age_matrix(2,2) = -1/5;
 age_matrix(3,2) = 1/5;age_matrix(3,3) = -1/5;
 age_matrix(4,3) = 1/5;age_matrix(4,4) = -1/10;
 age_matrix(5,4) = 1/10;age_matrix(5,5) = -1/10;
 age_matrix(6,5) = 1/10;age_matrix(6,6) = -1/10;
 age_matrix(7,6) = 1/10;age_matrix(7,7) = -1/10;
 age_matrix(8,7) = 1/10;age_matrix(8,8) = -1/10;
 age_matrix(9,8) = 1/10;
 % ***********************************************
%age_matrix=zeros(9,9);
 
 curr_mort_pwid=[0.96 0.96 1.12 0.18 0.22 0.53 1.38 4.28 14.96]/1000; % mortality per year for PWID
 curr_mort_former=[0.044 0.051 0.062 0.1 0.222 0.534 1.376 4.282 14.956 ]/1000; % mortality per year for PWID
 
  % ***********************************************
%  curr_mort_pwid = zeros(1,9);
%  curr_mort_former = zeros(1,9);
 mort_current = repmat(curr_mort_pwid,20,1);
 mort_former = repmat(curr_mort_former,20,1);
 

 
 extra_parms_nams={'piv','relapse','nu'}; % infection rate (piv instead of pi), relapse to IDU, 1/duration of injecting span
 extra_parms_vals=[1.58264*0.056,-log(1-0.027),1/17]; % was 5.6, nu was 1/17
 
 
 % starting values t = 1950:2030
 % S001 = 560 # formwer PWID compartment 161
 % S101 = 400 # current PWID compartment 361
 % F101 = 40 # acutely infected with HCV  366 (A is comp 6 in X list)
 N0=zeros(20*9*4,1);
 %N0(161)=560;
 N0(1)=560;
 N0(361)=396;
 N0(367)=44;%44;
 
 
 % compartment order is 00 (1-20, 21 -40, etc to 161-180)
 % compartment order is 01 180+(1-20, 21 -40, etc to 161-180)
 % compartment order is 10 360+(1-20, 21 -40, etc to 161-180)
 % compartment order is 11 540+(1-20, 21 -40, etc to 161-180)
 
 

%  M=zeros(20,20);
%  M(1,1)=-phi;M(1,6)=delta*r_AF0;M(1,16)=w;
%  M(2,2)=-phi;M(2,17)=w;
%  M(3,3)=-phi;M(3,18)=w;
%  M(4,4)=-phi;M(4,19)=w;
%  M(5,5)=-phi-r_SVR4DC-r_SVR4HCC;M(5,20)=w;
%  M(6,1)=phi;M(6,6)=-R_AF0;
%  M(7,6)=(1-delta)*r_AF0;M(7,8)=-r_F0F1;
%  M(8,2)=phi;M(8,7)=r_F0F1;M(8,8)=-r_F1F2;
%  M(9,3)=phi;M(9,8)=r_F1F2;M(9,9)=-r_F2F3;
%  M(10,4)=phi;M(10,9)=r_F2F3;M(10,10)=-r_F3F4;
%  M(11,5)=phi;M(11,10)=r_F3F4;M(11,11)=-r_F4DC-r_F4HCC;
%  M(12,11)=r_F4DC+r_SVR4DC;M(12,12) = -r_DCHCC-r_DCLT-r_DCdeath;
%  M(13,11)=r_F4HCC+r_SVR4HCC;M(13,12)=r_DCHCC;M(13,13)=-r_HCCLT-r_HCCdeath;
%  M(14,12)=r_DCLT;M(14,13)=r_HCCLT;M(14,14)=-1-r_LT1death;
%  M(15,14)=1;M(15,15)=r_LT2death;
%  M(16,16)=-w;
%  M(17,17)=-w;
%  M(18,18)=-w;
%  M(19,19)=-w;
%  M(20,20)=-w;