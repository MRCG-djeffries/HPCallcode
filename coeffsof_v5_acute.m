function [compnames,chronic_nams,...
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
          T4_comps_00,T4_comps_10,T4_comps_01,T4_comps_11,...
          TA_comps_00,TA_comps_10,TA_comps_01,TA_comps_11]=coeffsof_v5_acute
% compnames - root name sof the 20 compartments
% chronic_nams - names of ij ordering
% chronic_numsij - compartment index for all chronically infected with HCV


% Gives the index numbers for categories
% There are 20*9*4 compartments
% 20 represents the X notation
% element 1 to 5 S0,S1,S2,S3,S4
% infection naïve or previously achieving spontaneous clearance or SVR through treatment from liver fibrosis stage F0 to F4 respectively
% elemment 6 is A acute stage
% element 7 to 11 F0,F1,F2,F3,F4 - the fibrosis stages, chronic
% element 12 DC stage, chronic
% element 13 HCC stage, chronic
% element 14 and 15 LT1 and LT2 - liver transplant stages, chronic
% element 16 to 20 T0 to T4 
% chronically infected and in treatment achieving sustained viral response (SVR) (T0 to T4—treated from liver fibrosis stage F0 to F4 respectively)
compnames={'S0','S1','S2','S3','S4','A','F0','F1','F2','F3','F4','DC','HCC','LT1','LT2', 'T0','T1','T2','T3','T4','TA'};
% each elemet has three subscripts i,j,k
% i = 0 formwerPWID, i = 1 current PWID
% j = 0 never failed treatment, j = 1 had treatment failure
% k = 1 to 9 the 9 age groups

% for the 20*9*4 compartments
% ordering is 
% (i=0,j=0,k = 1);(i=0,j=0,k = 2);...;(i=0,j=0,k = 9) former, never failed
% (i=0,j=1,k = 1);(i=0,j=0,k = 2);...;(i=0,j=0,k = 9) former, failed
% (i=1,j=0,k = 1);(i=0,j=0,k = 2);...;(i=0,j=0,k = 9) current, never failed
% (i=1,j=1,k = 1);(i=0,j=0,k = 2);...;(i=0,j=0,k = 9) current, failed

chronic_numsA =7:15;
chronic_nums00 = chronic_numsA;
chronic_nums01 = 189+chronic_numsA;
chronic_nums10 = 378+chronic_numsA;
chronic_nums11 = 567+chronic_numsA;

for i = 2 : 9
    chronic_nums00 = [chronic_nums00 21*(i-1)+chronic_numsA];
    chronic_nums01 = [chronic_nums01 189+21*(i-1)+chronic_numsA];
    chronic_nums10 = [chronic_nums10 378+21*(i-1)+chronic_numsA];
    chronic_nums11 = [chronic_nums11 567+21*(i-1)+chronic_numsA];    
end

chronic_nams={'00=former,never','01=former,failed','10=current,never','11=current,failed'};

A_comps_00=(6:21:189);
A_comps_01=A_comps_00+189;
A_comps_10=A_comps_00+378;
A_comps_11=A_comps_00+567;
A_comps=[A_comps_00,A_comps_01,A_comps_10,A_comps_11];

F0_comps_00 = (7:21:189);
F0_comps_01 = F0_comps_00+189;
F0_comps_10 = F0_comps_00+378;
F0_comps_11 = F0_comps_00+567;
F0_comps=[F0_comps_00,F0_comps_01,F0_comps_10,F0_comps_11];

F1_comps_00 = (8:21:189);
F1_comps_01 = F1_comps_00+189;
F1_comps_10 = F1_comps_00+378;
F1_comps_11 = F1_comps_00+567;
F1_comps=[F1_comps_00,F1_comps_01,F1_comps_10,F1_comps_11];

F2_comps_00 = (9:21:189);
F2_comps_01 = F2_comps_00+189;
F2_comps_10 = F2_comps_00+378;
F2_comps_11 = F2_comps_00+567;
F2_comps=[F2_comps_00,F2_comps_01,F2_comps_10,F2_comps_11];

F3_comps_00 = (10:21:189);
F3_comps_01 = F3_comps_00+189;
F3_comps_10 = F3_comps_00+378;
F3_comps_11 = F3_comps_00+567;
F3_comps=[F3_comps_00,F3_comps_01,F3_comps_10,F3_comps_11];

F4_comps_00 = (11:21:189);
F4_comps_01 = F4_comps_00+189;
F4_comps_10 = F4_comps_00+378;
F4_comps_11 = F4_comps_00+567;
F4_comps=[F4_comps_00,F4_comps_01,F4_comps_10,F4_comps_11];

DC_comps_00 = (12:21:189);
DC_comps_01 = DC_comps_00+189;
DC_comps_10 = DC_comps_00+378;
DC_comps_11 = DC_comps_00+567;
DC_comps=[DC_comps_00,DC_comps_01,DC_comps_10,DC_comps_11];

HCC_comps_00 = (13:21:189);
HCC_comps_01 = HCC_comps_00+189;
HCC_comps_10 = HCC_comps_00+378;
HCC_comps_11 = HCC_comps_00+567;
HCC_comps=[HCC_comps_00,HCC_comps_01,HCC_comps_10,HCC_comps_11];

LT1_comps_00 = (14:21:189);
LT1_comps_01 = LT1_comps_00+189;
LT1_comps_10 = LT1_comps_00+378;
LT1_comps_11 = LT1_comps_00+567;
LT1_comps=[LT1_comps_00,LT1_comps_01,LT1_comps_10,LT1_comps_11];

LT2_comps_00 = (15:21:189);
LT2_comps_01 = LT2_comps_00+189;
LT2_comps_10 = LT2_comps_00+378;
LT2_comps_11 = LT2_comps_00+567;
LT2_comps=[LT2_comps_00,LT2_comps_01,LT2_comps_10,LT2_comps_11];


S0_comps_00 = (1:21:189);
S0_comps_01 = S0_comps_00+189;
S0_comps_10 = S0_comps_00+378;
S0_comps_11 = S0_comps_00+567;
S0_comps=[S0_comps_00,S0_comps_01,S0_comps_10,S0_comps_11];
 
S1_comps_00 = (2:21:189);
S1_comps_01 = S1_comps_00+189;
S1_comps_10 = S1_comps_00+378;
S1_comps_11 = S1_comps_00+567;
S1_comps=[S1_comps_00,S1_comps_01,S1_comps_10,S1_comps_11];
 
S2_comps_00 = (3:21:189);
S2_comps_01 = S2_comps_00+189;
S2_comps_10 = S2_comps_00+378;
S2_comps_11 = S2_comps_00+567;
S2_comps=[S2_comps_00,S2_comps_01,S2_comps_10,S2_comps_11];
 
S3_comps_00 = (4:21:189);
S3_comps_01 = S3_comps_00+189;
S3_comps_10 = S3_comps_00+378;
S3_comps_11 = S3_comps_00+567;
S3_comps=[S3_comps_00,S3_comps_01,S3_comps_10,S3_comps_11];
 
S4_comps_00 = (5:21:189);
S4_comps_01 = S4_comps_00+189;
S4_comps_10 = S4_comps_00+378;
S4_comps_11 = S4_comps_00+567;
S4_comps=[S4_comps_00,S4_comps_01,S4_comps_10,S4_comps_11];

T0_comps_00 = (16:21:189);
T0_comps_01 = T0_comps_00+189;
T0_comps_10 = T0_comps_00+378;
T0_comps_11 = T0_comps_00+567;
T0_comps=[T0_comps_00,T0_comps_01,T0_comps_10,T0_comps_11];
 
T1_comps_00 = (17:21:189);
T1_comps_01 = T1_comps_00+189;
T1_comps_10 = T1_comps_00+378;
T1_comps_11 = T1_comps_00+567;
T1_comps=[T1_comps_00,T1_comps_01,T1_comps_10,T1_comps_11];
 
T2_comps_00 = (18:21:189);
T2_comps_01 = T2_comps_00+189;
T2_comps_10 = T2_comps_00+378;
T2_comps_11 = T2_comps_00+567;
T2_comps=[T2_comps_00,T2_comps_01,T2_comps_10,T2_comps_11];
 
T3_comps_00 = (19:21:189);
T3_comps_01 = T3_comps_00+189;
T3_comps_10 = T3_comps_00+378;
T3_comps_11 = T3_comps_00+567;
T3_comps=[T3_comps_00,T3_comps_01,T3_comps_10,T3_comps_11];
 
T4_comps_00 = (20:21:189);
T4_comps_01 = T4_comps_00+189;
T4_comps_10 = T4_comps_00+378;
T4_comps_11 = T4_comps_00+567;
T4_comps=[T4_comps_00,T4_comps_01,T4_comps_10,T4_comps_11];

TA_comps_00 = (21:21:189);
TA_comps_01 = TA_comps_00+189;
TA_comps_10 = TA_comps_00+378;
TA_comps_11 = TA_comps_00+567;
TA_comps=[TA_comps_00,TA_comps_01,TA_comps_10,TA_comps_11];



