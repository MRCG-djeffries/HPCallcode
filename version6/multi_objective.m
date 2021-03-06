function multi_objective

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
fun = @(x)objval(x, N0,M,Mdash,mort_current,mort_former,phiminusvals, phiplussvals,age_matrix,chronic_nums00,chronic_nums10);
npts = 100;
opts_ga = optimoptions('gamultiobj','Display','iter','PopulationSize',2*npts,'FunctionTolerance',0.001);%,'ParetoFraction',.1);
A = []; b = [];
Aeq = []; beq = [];
lb=[0.01,0.001,1/20];
ub=[0.5,0.1,1.2];
numberOfVariables = 3;

% % pareto
% optsp = optimoptions('paretosearch','Display','iter','UseVectorized',false);%'PlotFcn',{'psplotparetof' 'psplotparetox'});
% % optsp.PlotFcn = [];
% % optsp.ParetoSetSize = 35;
% [xpsh,fvalpsh,~,pshoutput] = paretosearch(fun,numberOfVariables,A,b,Aeq,beq,lb,ub,[],optsp);



% call gamultiobj function
[x_ga,fval_ga,exitflag,gaoutput1] = gamultiobj(fun,numberOfVariables,A,b,Aeq,beq,lb,ub,opts_ga);
%z=ga(fun,numberOfVariables,A,b,Aeq,beq,lb,ub,[],opts_ga);
save('F:\HCV\version6\multi.mat' ,'x_ga','fval_ga')
df=45

function F = objval(x, N0,M,Mdash,mort_current,mort_former,phiminusvals, phiplussvals,age_matrix,chronic00,chronic10)  
extra_parms_vals = x;
[~,x1] = ode45(@odeeq_v5,0:81,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix);
X=x1'; 
% f1 = (sum(X(chronic10,66))/sum(X(chronic00,66))-0.5)^2;
% f2 = (sum(X(chronic10,66))/sum(X(361:540,66))-44/144)^2;
% 40,000 infected PWID / 144,000 infected former
 f1 = (sum(X(chronic10,66))/sum(X(chronic00,66))-40/144)^2;
 f2 = (sum(X(chronic10,66))/sum(X(361:540,66))-0.5)^2;
F = [f1, f2];
%F=f1+f2;