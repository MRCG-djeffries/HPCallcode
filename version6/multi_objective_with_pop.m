function multi_objective_with_pop
% use the output of multi_objective
% and looks at the starting population to get the actual populations
% 


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
load('F:\HCV\version6\multi.mat')
%save('F:\HCV\version6\multi.mat' ,'x_ga','fval_ga')
extra_parms_vals=x_ga(1,:);
N0(1)=560;N0(361)=435.6;N0(367)=4.4;%4.36;
N0=404.8911*N0;
[F,X,t1] = objval(extra_parms_vals, N0,M,Mdash,mort_current,mort_former,phiminusvals, phiplussvals,age_matrix,chronic_nums00,chronic_nums10);
% chronic are plot
lt2=sum(X(LT2_comps_00,:))+sum(X(LT2_comps_10,:));
lt1=sum(X(LT1_comps_00,:))+sum(X(LT1_comps_10,:));
lt=lt1+lt2;
dc=sum(X(DC_comps_00,:))+sum(X(DC_comps_10,:));
hcc=sum(X(HCC_comps_00,:))+sum(X(HCC_comps_10,:));
f4=sum(X(F4_comps_00,:))+sum(X(F4_comps_10,:));
f3=sum(X(F3_comps_00,:))+sum(X(F3_comps_10,:));
f2=sum(X(F2_comps_00,:))+sum(X(F2_comps_10,:));
f1=sum(X(F1_comps_00,:))+sum(X(F1_comps_10,:));
f0=sum(X(F0_comps_00,:))+sum(X(F0_comps_10,:));
a0=sum(X(A_comps_00,:))+sum(X(A_comps_10,:));

figure;area(t1',[lt' hcc' dc' f4' f3' f2' f1' f0'])
xticks(0:5:81)
set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
legend({'lt' 'hcc' 'dc' 'f4' 'f3' 'f2' 'f1' 'f0'},'location','best');
grid()
title('Number of PWID at each stage of chronic HCV')
set(gca,'xlim',[0 81])

chronic_prev_current = chronic_nums10; % with no treatments chronic_nums11 all zero
total_current = 361:540;
figure;plot(t1,100*sum(X(chronic_prev_current,:))./sum(X(total_current,:)),'r--','LineWidth',2)
set(gca,'ylim',[0 70])
xticks(0:5:80)
set(gca,'xticklabel',1950:5:2031,'xticklabelrotation',90)
title('% of current PWID with HCV from 1950 to 2031')
set(gca,'xlim',[0 81])
grid()

bot = sum(X(total_current,:));
%bot = sum(X(1:180,:));


% cols are f0,f1,f2,f3,f4,dc,hcc,lt1,lt2
% rows are the age
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
%agemat=agemat-360;
mn=size(X);
ageprop=zeros(9,mn(2));
cols  = lines(9);
H = zeros(9,1);  % Store handles
figure;
for i = 1 : 9
    ageprop(i,:) = sum(X(agemat(i,2:10),:));
    H(i)=plot(t1,100*ageprop(i,:)./bot,'-','LineWidth',2,'Color',cols(i,:));
    hold on
end
set(gca,'ylim',[0 90])
xticks(0:5:80)
set(gca,'xticklabel',1950:5:2031,'xticklabelrotation',90)
title('% of current PWID with HCV from 1950 to 2031 by age')
set(gca,'xlim',[0 81])
grid()
legend(H, {'20-24','25-29','30-34','35-44','45-54','55-64','65-74','75-84','85+'},'location','best');




npts = 100;
opts_ga = optimoptions('gamultiobj','Display','iter','PopulationSize',2*npts,'FunctionTolerance',0.001);%,'ParetoFraction',.1);
A = []; b = [];
Aeq = []; beq = [];
lb=[0.01,0.001,1/20];
ub=[0.5,0.1,1.2];
numberOfVariables = 3;

% pareto
optsp = optimoptions('paretosearch','Display','iter','UseVectorized',false,'PlotFcn',{'psplotparetof' 'psplotparetox'});
% optsp.PlotFcn = [];
% optsp.ParetoSetSize = 35;
[xpsh,fvalpsh,~,pshoutput] = paretosearch(fun,numberOfVariables,A,b,Aeq,beq,lb,ub,[],optsp);



% call gamultiobj function
[x_ga,fval_ga,exitflag,gaoutput1] = gamultiobj(fun,numberOfVariables,A,b,Aeq,beq,lb,ub,opts_ga);
%z=ga(fun,numberOfVariables,A,b,Aeq,beq,lb,ub,[],opts_ga);
save('F:\HCV\version6\multi3.mat' ,'x_ga','fval_ga')
df=45

function [F,X,t1] = objval(x, N0,M,Mdash,mort_current,mort_former,phiminusvals, phiplussvals,age_matrix,chronic00,chronic10)  
extra_parms_vals = x;
[t1,x1] = ode45(@odeeq_v5,0:81,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix);
X=x1'; 
% f1 = (sum(X(chronic10,66))/sum(X(chronic00,66))-0.5)^2;
% f2 = (sum(X(chronic10,66))/sum(X(361:540,66))-44/144)^2;
% 40,000 infected PWID / 144,000 infected former
%  f1 = (sum(X(chronic10,66))/sum(X(chronic00,66))-40/144)^2;
%  f2 = (sum(X(chronic10,66))/sum(X(361:540,66))-0.5)^2;
 f1 = sum(X(chronic10,66)); % current HCV-infected PWID 40k
 f2 = sum(X(361:540,66)); % all PWID 80k
 f3 = sum(X(chronic00,66)); % former HCV-infected PWID 144k
F = [f1, f2,f3];
%F=f1+f2;