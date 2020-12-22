function cost_effective_v2
% 24 Spetember 2020
% try different ode solvers
% total_treat function
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
%N0(1)=560;N0(361)=440;
%N0(367:20:540)=[12 4*ones(1,8)];
extra_parms_vals(1)=0.133;% was 113
[t1,x1] = ode45(@odeeq_v5,0:81,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix);
X=x1';
figure;plot(t1,(1-0.26)*2.1667*(sum(X(A_comps_00,:))+sum(X(A_comps_10,:))))
 
% chronic HCV prevlance calibrated at 50% for current PWID in 2015.
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
mn=size(X);
ageprop=zeros(9,mn(2));
cols  = lines(9);
H = zeros(9,1);  % Store handles
figure;
for i = 1 : 9
    ageprop(i,:) = sum(X(agemat(i,:),:));
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

H = zeros(9,1);  % Store handles
age_weights=zeros(9,1);

figure;
for i = 1 : 9
    if i == 1
        age_weights(i) = 0.25*bot(66)/ageprop(1,66);
    else
        age_weights(i) = (0.75/8)*bot(66)/ageprop(i,66);
    end
    ageprop(i,:) = sum(X(agemat(i,:),:));
    H(i)=plot(t1,100*age_weights(i)*ageprop(i,:)./bot,'-','LineWidth',2,'Color',cols(i,:));
    hold on
end

A=repmat(repelem(age_weights,9),1,82);
B=repmat(repelem(age_weights,20),1,82);
set(gca,'ylim',[0 30])
xticks(0:5:80)
set(gca,'xticklabel',1950:5:2031,'xticklabelrotation',90)
title('% of current PWID with HCV from 1950 to 2031 age weighted')
set(gca,'xlim',[0 81])
grid()
legend(H, {'20-24','25-29','30-34','35-44','45-54','55-64','65-74','75-84','85+'},'location','best');
%age_weights_matrix
figure;plot(t1,100*sum(A.*X(chronic_prev_current,:))./sum(B.*X(total_current,:)),'r-','LineWidth',2)
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
    H(i)=plot(t1,100*ageprop(i,:)./bot,'-','LineWidth',2,'Color',cols(i,:));
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
a0=scale00*sum(X(A_comps_00,:))+scale10*sum(X(A_comps_10,:));
s4=scale00*sum(X(S4_comps_00,:))+scale10*sum(X(S4_comps_10,:));
s3=scale00*sum(X(S3_comps_00,:))+scale10*sum(X(S3_comps_10,:));
s2=scale00*sum(X(S2_comps_00,:))+scale10*sum(X(S2_comps_10,:));
s1=scale00*sum(X(S1_comps_00,:))+scale10*sum(X(S1_comps_10,:));
s0=scale00*sum(X(S0_comps_00,:))+scale10*sum(X(S0_comps_10,:));
base_tot=s0+s1+s2+s3+s4+a0+f0+f1+f2+f3+f4+hcc+dc+lt1+lt2;
figure;area(t1,[lt' hcc' dc' f4' f3' f2' f1' f0'])
xticks(0:5:81)
set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
legend({'lt' 'hcc' 'dc' 'f4' 'f3' 'f2' 'f1' 'f0'},'location','best');
grid()
title('Number of PWID at each stage of chronic HCV')
set(gca,'xlim',[0 81])

% s compartments
s4=scale00*sum(X(S4_comps_00,:))+scale10*sum(X(S4_comps_10,:));
s3=scale00*sum(X(S3_comps_00,:))+scale10*sum(X(S3_comps_10,:));
s2=scale00*sum(X(S2_comps_00,:))+scale10*sum(X(S2_comps_10,:));
s1=scale00*sum(X(S1_comps_00,:))+scale10*sum(X(S1_comps_10,:));
s0=scale00*sum(X(S0_comps_00,:))+scale10*sum(X(S0_comps_10,:));

figure;area(t1,[s4' s3' s2' s1' s0'])
xticks(0:5:81)
set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
legend({'s4' 's3' 's2' 's1' 's0'},'location','best');
grid()
title('Number of PWID at S compartments')
set(gca,'xlim',[0 82])

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
set(gca,'xlim',[67 82])

% s compartments
s4=age_scale00*sum(Aformer.*X(S4_comps_00,:))+age_scale10*sum(Apwid.*X(S4_comps_10,:));
s3=age_scale00*sum(Aformer.*X(S3_comps_00,:))+age_scale10*sum(Apwid.*X(S3_comps_10,:));
s2=age_scale00*sum(Aformer.*X(S2_comps_00,:))+age_scale10*sum(Apwid.*X(S2_comps_10,:));
s1=age_scale00*sum(Aformer.*X(S1_comps_00,:))+age_scale10*sum(Apwid.*X(S1_comps_10,:));
s0=age_scale00*sum(Aformer.*X(S0_comps_00,:))+age_scale10*sum(Apwid.*X(S0_comps_10,:));
a0=age_scale00*sum(Aformer.*X(A_comps_00,:))+age_scale10*sum(Apwid.*X(A_comps_10,:));

CC=[sum(s0(67:end)) ,sum(s1(67:end)) ,sum(s2(67:end)) ,sum(s3(67:end)) ,sum(s4(67:end)),...
    sum(a0(67:end)),...
    sum(f0(67:end)) ,sum(f1(67:end)) ,sum(f2(67:end)) ,sum(f3(67:end)) ,sum(f4(67:end)),...
    sum(dc(67:end)) ,sum(hcc(67:end)), sum(lt(67:end)),...
    0 0 0 0 0];

CC2=[sum(s0(1:65)) ,sum(s1(1:65)) ,sum(s2(1:65)) ,sum(s3(1:65)) ,sum(s4(1:65)),...
    sum(a0(1:65)),...
    sum(f0(1:65)) ,sum(f1(1:65)) ,sum(f2(1:65)) ,sum(f3(1:65)) ,sum(f4(1:65)),...
    sum(dc(1:65)) ,sum(hcc(1:65)), sum(lt(1:65)),...
    0 0 0 0 0];
totinmodelbase=s0(67:end)+s1(67:end)+s2(67:end)+s3(67:end)+s4(67:end)+...
    a0(67:end)+...
    f0(67:end)+f1(67:end)+f2(67:end)+f3(67:end)+f4(67:end)+...
    dc(67:end)+hcc(67:end)+lt(67:end);

figure;area(t1,[s4' s3' s2' s1' s0'])
xticks(0:5:81)
set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
legend({'s4' 's3' 's2' 's1' 's0'},'location','best');
grid()
title('Number of PWID at S compartments')
set(gca,'xlim',[0 82])



%% vector defs
compvec10= [S0_comps_10;A_comps_10;F0_comps_10;F1_comps_10;F2_comps_10;F3_comps_10;F4_comps_10;HCC_comps_10;DC_comps_10;LT1_comps_10;LT2_comps_10];
compvec00= [S0_comps_00;A_comps_00;F0_comps_00;F1_comps_00;F2_comps_00;F3_comps_00;F4_comps_00;HCC_comps_00;DC_comps_00;LT1_comps_00;LT2_comps_00];
utility_vec=[0.93 0.77 0.77 0.77 0.77 0.66 0.55 0.45 0.45 0.45 0.67]; % These are the utiltiy multipliers for
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
Tutility_vec=[0.93 0.93 0.93 0.93 0.93 ...
    0.77 ...
    0.77 0.77 0.77 0.66 0.55 ...
    0.45 0.45 0.45 0.67 ...
    0.77 0.77 0.77 0.66 0.55];%Number of PWDI

%% calculate the QALY
qualy_vals10=zeros(11,1); % This conatins PWID qualy for each chronic department
qualy_vals00=zeros(11,1); % This conatins former qualy for each chronic department
Apwid=repmat(repelem(age_weights,1),1,82);
Aformer=repmat(repelem(age_weights_former,1),1,82);
time2016=67; % for scale 1950:2031;
for i = 1 : 11
    Xtest=age_scale10*sum(Apwid.*X(compvec10(i,:),:));
    Ytest=age_scale00*sum(Aformer.*X(compvec00(i,:),:));
    %testy=testy+Xtest(time2015);
    qualy_vals10(i)=calc_qalyv2(Xtest(time2016:end),2016,2031,1,utility_vec(i),3); % 2000 to 2030 in steps of 1 month
    qualy_vals00(i)=calc_qalyv2(Ytest(time2016:end),2016,2031,1,utility_vec(i),3); % 2000 to 2030 in steps of 1 month
end
total_qual_current = sum(qualy_vals10)+sum(qualy_vals00);



Tqualy_vals00=zeros(20,1); % This conatins qualy for all compartments 00 (former, not failed)
Tqualy_vals01=zeros(20,1); % This conatins qualy for all compartments 01 (former, failed)
Tqualy_vals10=zeros(20,1); % This conatins qualy for all compartments 10 (current, not failed)
Tqualy_vals11=zeros(20,1); % This conatins qualy for all compartments 11 (current, failed)
XX=0;
for i = 1 : 20
    Xtest00=age_scale00*sum(Aformer.*X(Tcompvec00(i,:),:));
    Xtest01=age_scale00*sum(Aformer.*X(Tcompvec01(i,:),:));
    Xtest10=age_scale10*sum(Apwid.*X(Tcompvec10(i,:),:));
    Xtest11=age_scale10*sum(Apwid.*X(Tcompvec11(i,:),:));
    %     Xtest00=scale00*sum(X(Tcompvec00(i,:),:));
    %     Xtest01=scale00*sum(X(Tcompvec01(i,:),:));
    %     Xtest10=scale10*sum(X(Tcompvec10(i,:),:));
    %     Xtest11=scale10*sum(X(Tcompvec11(i,:),:));
    XX=XX+Xtest00+Xtest01+Xtest10+Xtest11;
    Tqualy_vals00(i)=calc_qalyv2(Xtest00(time2016:end),2016,2031,1,Tutility_vec(i),3); % 2000 to 2030 in steps of 1 month
    Tqualy_vals01(i)=calc_qalyv2(Xtest01(time2016:end),2016,2031,1,Tutility_vec(i),3); % 2000 to 2030 in steps of 1 month
    Tqualy_vals10(i)=calc_qalyv2(Xtest10(time2016:end),2016,2031,1,Tutility_vec(i),3); % 2000 to 2030 in steps of 1 month
    Tqualy_vals11(i)=calc_qalyv2(Xtest11(time2016:end),2016,2031,1,Tutility_vec(i),3); % 2000 to 2030 in steps of 1 month
end
total_qual_current_scaled = sum(Tqualy_vals00)+sum(Tqualy_vals01)+sum(Tqualy_vals10)+sum(Tqualy_vals11);
figure;plot(t1,sum(X,1))
set(gca,'ylim',[900 1100])
title('Should be constant - tests for leaks')

% the number of incident cases
% above for scaled
rate=parm_current(2);
delta=parm_current(1);
T00_scaled=(1-delta)*rate*age_scale00*sum(Aformer.*X(A_comps_00,:));
T10_scaled=(1-delta)*rate*age_scale10*sum(Apwid.*X(A_comps_10,:));
% integrate
inci_acute00_scaled=calc_qalyv2(T00_scaled(time2016:82),2016,2031,1,1,0); % 2000 to 2030 in steps of 1 month
inci_acute10_scaled=calc_qalyv2(T10_scaled(time2016:82),2016,2031,1,1,0); % 2000 to 2030 in steps of 1 month
yearly_diagnosed_cases_in_chronic_stages_scaled=inci_acute00_scaled+inci_acute10_scaled;
testcompvec10= [DC_comps_10;HCC_comps_10;LT1_comps_10;LT2_comps_10];
testcompvec00= [DC_comps_00;HCC_comps_00;LT1_comps_00;LT2_comps_00];
death_vec= [-log(1-0.138) -log(1-0.605) -log(1-0.169) -log(1-0.034)]; %DC,HCC,LT1,LT2
fatal_QALY_scaled=fatalQALYS(age_scale00,age_scale10,age_weights,age_weights_former,X,time2016,testcompvec10,testcompvec00,death_vec);
[costanual,costoneoff,namesof_anualcost,namesof_oneoff]=cost_params();
total_costs_scaled=base_cost(age_weights,age_weights_former,age_scale10,age_scale00,X,...
    time2016,costanual,costoneoff,parm_current,...
    DC_comps_00,DC_comps_10,HCC_comps_00,HCC_comps_10,...
    A_comps_00,A_comps_10,...
    F0_comps_10,F1_comps_10,F2_comps_10,...
    F0_comps_00,F1_comps_00,F2_comps_00,...
    F3_comps_00,F3_comps_10,...
    F4_comps_00,F4_comps_10);

% same as above don't need to do this
% with interval of 1/12
% [t1,x1] = ode45(@odeeq_v5,0:1/12:81,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix);
% X=x1';
% for 2000 to 2031 with 0:1/12:80 this is 601:973
% in this scale time2015 is 2015 whichis index 781 in above line
% Nite this is 2016   was 781
% the compartments ordered in compvec, Table B3 Scott supplement
% scaline is scale10=40000/sum(X([chronic_nums10],66)); % current
% get individual multipliers
% ********** This is for current

% qualy_vals10=zeros(11,1); % This conatins PWID qualy for each chronic department
% qualy_vals00=zeros(11,1); % This conatins former qualy for each chronic department
% time2015 = 793;
% Apwid=repmat(repelem(age_weights,1),1,973);
% Aformer=repmat(repelem(age_weights_former,1),1,973);
% for i = 1 : 11
%     Xtest=age_scale10*sum(Apwid.*X(compvec10(i,:),:));
%     Ytest=age_scale00*sum(Aformer.*X(compvec00(i,:),:));
%     %testy=testy+Xtest(time2015);
%     qualy_vals10(i)=calc_qalyv2(Xtest(time2015:end),2016,2031,1/12,utility_vec(i),3); % 2000 to 2030 in steps of 1 month
%     qualy_vals00(i)=calc_qalyv2(Ytest(time2015:end),2016,2031,1/12,utility_vec(i),3); % 2000 to 2030 in steps of 1 month
% end
% total_qual_current = sum(qualy_vals10)+sum(qualy_vals00);
[~,tot_liver_deaths_scaled,~]=fatalQALYS(age_scale00,age_scale10,age_weights,age_weights_former,X,time2016,testcompvec10,testcompvec00,death_vec);
[ltb,dcb,hccb,f4b,f3b,f2b,f1b,f0b,s4b,s3b,s2b,s1b,s0b,t4b,t3b,t2b,t1b,t0b,a0b...
    ]=inallcomps(age_scale00,age_scale10,age_weights,age_weights_former,X,...
    LT1_comps_00,LT1_comps_10,LT1_comps_01,LT1_comps_11,...
    LT2_comps_00,LT2_comps_10,LT2_comps_01,LT2_comps_11,...
    DC_comps_00,DC_comps_10,DC_comps_01,DC_comps_11,...
    HCC_comps_00,HCC_comps_10,HCC_comps_01,HCC_comps_11,...
    F4_comps_00,F4_comps_10,F4_comps_01,F4_comps_11,...
    F3_comps_00,F3_comps_10,F3_comps_01,F3_comps_11,...
    F2_comps_00,F2_comps_10,F2_comps_01,F2_comps_11,...
    F1_comps_00,F1_comps_10,F1_comps_01,F1_comps_11,...
    F0_comps_00,F0_comps_10,F0_comps_01,F0_comps_11,...
    S4_comps_00,S4_comps_10,S4_comps_01,S4_comps_11,...
    S3_comps_00,S3_comps_10,S3_comps_01,S3_comps_11,...
    S2_comps_00,S2_comps_10,S2_comps_01,S2_comps_11,...
    S1_comps_00,S1_comps_10,S1_comps_01,S1_comps_11,...
    S0_comps_00,S0_comps_10,S0_comps_01,S0_comps_11,...
    T4_comps_00,T4_comps_10,T4_comps_01,T4_comps_11,...
    T3_comps_00,T3_comps_10,T3_comps_01,T3_comps_11,...
    T2_comps_00,T2_comps_10,T2_comps_01,T2_comps_11,...
    T1_comps_00,T1_comps_10,T1_comps_01,T1_comps_11,...
    T0_comps_00,T0_comps_10,T0_comps_01,T0_comps_11,...
    A_comps_00,A_comps_10,A_comps_01,A_comps_11);
%% treatment
pcom=0.892; % prob of current PWID completing treatment
alpha=0.95; % prob of acheiving SVR
opty = 0;
if opty==1
n0 = 5506; % was 4906
n1 = 7400; % was 6000
y1=73; % above is applied for 6 years after 2015 (year 66)
% optimization
%fun = @optimization_approach;
x0 = [5000,6000];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [1000,10000];
ub = [10000,10000];
nonlcon = [];
 options = optimoptions('fmincon','Display','iter-detailed');
% xvalout = fmincon(@(x)optimization_approach(x,N0,M,Mdash,mort_current,mort_former,extra_parms_vals,...
%     phiminusvals, phiplussvals,age_matrix,...
%     age_scale00,age_scale10,pcom,alpha,age_weights,age_weights_former,death_vec,...
%     delta,rate,A_comps_00,A_comps_01,A_comps_10,A_comps_11,...
%     testcompvec10,testcompvec00,y1),x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
% integer
options = optimoptions('ga','ConstraintTolerance',1e-6,'PlotFcn', @gaplotbestf);
xvalout = ga(@(x)optimization_approach(x,N0,M,Mdash,mort_current,mort_former,extra_parms_vals,...
    phiminusvals, phiplussvals,age_matrix,...
    age_scale00,age_scale10,pcom,alpha,age_weights,age_weights_former,death_vec,...
    delta,rate,A_comps_00,A_comps_01,A_comps_10,A_comps_11,...
    testcompvec10,testcompvec00,y1),2,A,b,Aeq,beq,lb,ub,nonlcon,[1 2],options);

% now adjust for death correction
% optimization
%fun = @optimization_approach;
x0 = [0.5];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0];
ub = [1.5];
nonlcon = [];
 options = optimoptions('fmincon','Display','iter-detailed');
 Xbase=X;
 totbase=ltb+dcb+hccb+f4b+f3b+f2b+f1b+f0b+s4b+s3b+s2b+s1b+s0b+t4b+t3b+t2b+t1b+t0b+a0b;
xvalout_death = fmincon(@(x)optimization_approach_death(x,N0,M,Mdash,mort_current,mort_former,extra_parms_vals,...
    phiminusvals, phiplussvals,age_matrix,...
    age_scale00,age_scale10,pcom,alpha,age_weights,age_weights_former,death_vec,...
    xvalout(1),xvalout(2),xvalout(3),...
    delta,rate,A_comps_00,A_comps_01,A_comps_10,A_comps_11,...
    LT1_comps_00,LT1_comps_10,LT1_comps_01,LT1_comps_11,...
    LT2_comps_00,LT2_comps_10,LT2_comps_01,LT2_comps_11,...
    DC_comps_00,DC_comps_10,DC_comps_01,DC_comps_11,...
    HCC_comps_00,HCC_comps_10,HCC_comps_01,HCC_comps_11,...
    F4_comps_00,F4_comps_10,F4_comps_01,F4_comps_11,...
    F3_comps_00,F3_comps_10,F3_comps_01,F3_comps_11,...
    F2_comps_00,F2_comps_10,F2_comps_01,F2_comps_11,...
    F1_comps_00,F1_comps_10,F1_comps_01,F1_comps_11,...
    F0_comps_00,F0_comps_10,F0_comps_01,F0_comps_11,...
    S4_comps_00,S4_comps_10,S4_comps_01,S4_comps_11,...
    S3_comps_00,S3_comps_10,S3_comps_01,S3_comps_11,...
    S2_comps_00,S2_comps_10,S2_comps_01,S2_comps_11,...
    S1_comps_00,S1_comps_10,S1_comps_01,S1_comps_11,...
    S0_comps_00,S0_comps_10,S0_comps_01,S0_comps_11,...
    T4_comps_00,T4_comps_10,T4_comps_01,T4_comps_11,...
    T3_comps_00,T3_comps_10,T3_comps_01,T3_comps_11,...
    T2_comps_00,T2_comps_10,T2_comps_01,T2_comps_11,...
    T1_comps_00,T1_comps_10,T1_comps_01,T1_comps_11,...
    T0_comps_00,T0_comps_10,T0_comps_01,T0_comps_11,...    
    testcompvec10,testcompvec00,totbase),x0,A,b,Aeq,beq,lb,ub,nonlcon,options);





%   [t1,x1t] = ode45(@odeeq_v5T_2p,0:81,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix,...
%                scale00,scale10,n0,n1,pcom,alpha,66,y1); % note 2015 marker is 65 years
% [t1,x1] = ode45(@odeeq_v5,0:81,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix);
end
gridy = 0;
if gridy==1
     n0=8000:100:9000;
     n1=3000:100:4000;
     Nn0=zeros(length(n0)*length(n1),1);
     Nn1=N0;
     Y_inc=N0;
     Y_dea=N0;
     Z_inc=zeros(length(n0),length(n1));
     Z_dea=Z_inc;
     XX=Z_inc;
     YY=Z_inc;
    y1=74;
death_val=0.4;%xvalout_death(1);
k = 0;
ik=0;
ij=0;
    for i = n0
        ik = ik + 1;
        ij=0;
        for j = n1
            ij=ij+1;
            k = k +1
            [t1,x1t] = ode45(@odeeq_v5T_2p_scaled_death,0:81,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix,...
             age_scale00,age_scale10,i,j,pcom,alpha,67,y1,age_weights',age_weights_former',death_val); % note 2015 marker is 65 years
XT=x1t';
Apwid=repmat(repelem(age_weights,1),1,82);
Aformer=repmat(repelem(age_weights_former,1),1,82);
T00=(1-delta)*rate*age_scale00*sum(Aformer.*XT(A_comps_00,:));
T01=(1-delta)*rate*age_scale00*sum(Aformer.*XT(A_comps_01,:));
T11=(1-delta)*rate*age_scale10*sum(Apwid.*XT(A_comps_11,:));
T10=(1-delta)*rate*age_scale10*sum(Apwid.*XT(A_comps_10,:));
% integrate
yinci_PWID=just_integrate(T10(time2016:82),2016,2031,1)+...
    just_integrate(T00(time2016:82),2016,2031,1)+...
    just_integrate(T01(time2016:82),2016,2031,1)+...
    just_integrate(T11(time2016:82),2016,2031,1);
[~,~,XTdeaths]=fatalQALYS(age_scale00,age_scale10,age_weights,age_weights_former,XT,time2016,testcompvec10,testcompvec00,death_vec);
% r2inci=(yinci_PWID(end)/yinci_PWID(1)-0.2)^2;
% d2inci=(XTdeaths(end)/XTdeaths(1)-0.35)^2;
         Nn0(k)=i;
         Nn1(k)=j;
         Y_inc(k)=yinci_PWID(end)/yinci_PWID(1);
         Y_dea(k)=XTdeaths(end)/XTdeaths(1);
         Z_inc(ik,ij)=Y_inc(k);
         Z_dea(ik,ij)=Y_dea(k);
         XX(ik,ij) = i;
         YY(ik,ij)=j;
        end
    end
end
n0=9000;%*2;%xvalout(1); % 17600
n1=5500;%*2;%xvalout(2);  %6600
%n0=0;n1=0;
y1=70;  % was 74 y1-67 + 1 is the number of years of intervention starting in 2016
death_val=0.35;%*1.75;%xvalout_death(1);

% % ICER was 14,301
% n0=100;%xvalout(1); % 17600
% n1=100;%xvalout(2);  %6600
% y1=70;  % was 74
% death_val=-0.024;%xvalout_death(1);
% 
% n0=1000;%xvalout(1); % 17600
% n1=1000;%xvalout(2);  %6600
% y1=70;  % was 74
% death_val=0.025;%xvalout_death(1);
% ***************************************************
%   note with very high treat levels - where every incidence case is
%   treated, and former compartmenst become empty, then becomes stiff
%   no change over some compartments
tic
[t1out,x1t] = ode15s(@odeeq_v5T_2p_scaled_death,0:81,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix,...
    age_scale00,age_scale10,n0,n1,pcom,alpha,67,y1,age_weights',age_weights_former',death_val); % note 2015 marker is 65 years
toc
XT=x1t';
%XT(XT<0)=0; % can be below tol i.e. very small 0s - messes up transformation
timey=t1out;
X=XT;
% scaled version
Apwid=repmat(age_weights,1,82);
Aformer=repmat(age_weights_former,1,82);
lt2=age_scale00*sum(Aformer.*X(LT2_comps_00,:))+age_scale10*sum(Apwid.*X(LT2_comps_10,:))+...
    age_scale00*sum(Aformer.*X(LT2_comps_01,:))+age_scale10*sum(Apwid.*X(LT2_comps_11,:));
lt1=age_scale00*sum(Aformer.*X(LT1_comps_00,:))+age_scale10*sum(Apwid.*X(LT1_comps_10,:))+...
    age_scale00*sum(Aformer.*X(LT1_comps_01,:))+age_scale10*sum(Apwid.*X(LT1_comps_11,:));
lt=lt1+lt2;
dc=age_scale00*sum(Aformer.*X(DC_comps_00,:))+age_scale10*sum(Apwid.*X(DC_comps_10,:))+...
    age_scale00*sum(Aformer.*X(DC_comps_01,:))+age_scale10*sum(Apwid.*X(DC_comps_11,:));
hcc=age_scale00*sum(Aformer.*X(HCC_comps_00,:))+age_scale10*sum(Apwid.*X(HCC_comps_10,:))+...
    age_scale00*sum(Aformer.*X(HCC_comps_01,:))+age_scale10*sum(Apwid.*X(HCC_comps_11,:));
f4=age_scale00*sum(Aformer.*X(F4_comps_00,:))+age_scale10*sum(Apwid.*X(F4_comps_10,:))+...
    age_scale00*sum(Aformer.*X(F4_comps_01,:))+age_scale10*sum(Apwid.*X(F4_comps_11,:));
f3=age_scale00*sum(Aformer.*X(F3_comps_00,:))+age_scale10*sum(Apwid.*X(F3_comps_10,:))+...
    age_scale00*sum(Aformer.*X(F3_comps_01,:))+age_scale10*sum(Apwid.*X(F3_comps_11,:));
f2=age_scale00*sum(Aformer.*X(F2_comps_00,:))+age_scale10*sum(Apwid.*X(F2_comps_10,:))+...
    age_scale00*sum(Aformer.*X(F2_comps_01,:))+age_scale10*sum(Apwid.*X(F2_comps_11,:));
f1=age_scale00*sum(Aformer.*X(F1_comps_00,:))+age_scale10*sum(Apwid.*X(F1_comps_10,:))+...
    age_scale00*sum(Aformer.*X(F1_comps_01,:))+age_scale10*sum(Apwid.*X(F1_comps_11,:));
f0=age_scale00*sum(Aformer.*X(F0_comps_00,:))+age_scale10*sum(Apwid.*X(F0_comps_10,:))+...
    age_scale00*sum(Aformer.*X(F0_comps_01,:))+age_scale10*sum(Apwid.*X(F0_comps_11,:));
s4=age_scale00*sum(Aformer.*X(S4_comps_00,:))+age_scale10*sum(Apwid.*X(S4_comps_10,:))+...
    age_scale00*sum(Aformer.*X(S4_comps_01,:))+age_scale10*sum(Apwid.*X(S4_comps_11,:));
s3=age_scale00*sum(Aformer.*X(S3_comps_00,:))+age_scale10*sum(Apwid.*X(S3_comps_10,:))+...
    age_scale00*sum(Aformer.*X(S3_comps_01,:))+age_scale10*sum(Apwid.*X(S3_comps_11,:));
s2=age_scale00*sum(Aformer.*X(S2_comps_00,:))+age_scale10*sum(Apwid.*X(S2_comps_10,:))+...
    age_scale00*sum(Aformer.*X(S2_comps_01,:))+age_scale10*sum(Apwid.*X(S2_comps_11,:));
s1=age_scale00*sum(Aformer.*X(S1_comps_00,:))+age_scale10*sum(Apwid.*X(S1_comps_10,:))+...
    age_scale00*sum(Aformer.*X(S1_comps_01,:))+age_scale10*sum(Apwid.*X(S1_comps_11,:));
s0=age_scale00*sum(Aformer.*X(S0_comps_00,:))+age_scale10*sum(Apwid.*X(S0_comps_10,:))+...
    age_scale00*sum(Aformer.*X(S0_comps_01,:))+age_scale10*sum(Apwid.*X(S0_comps_11,:));
t4=age_scale00*sum(Aformer.*X(T4_comps_00,:))+age_scale10*sum(Apwid.*X(T4_comps_10,:))+...
    age_scale00*sum(Aformer.*X(T4_comps_01,:))+age_scale10*sum(Apwid.*X(T4_comps_11,:));
t3=age_scale00*sum(Aformer.*X(T3_comps_00,:))+age_scale10*sum(Apwid.*X(T3_comps_10,:))+...
    age_scale00*sum(Aformer.*X(T3_comps_01,:))+age_scale10*sum(Apwid.*X(T3_comps_11,:));
t2=age_scale00*sum(Aformer.*X(T2_comps_00,:))+age_scale10*sum(Apwid.*X(T2_comps_10,:))+...
    age_scale00*sum(Aformer.*X(T2_comps_01,:))+age_scale10*sum(Apwid.*X(T2_comps_11,:));
t1=age_scale00*sum(Aformer.*X(T1_comps_00,:))+age_scale10*sum(Apwid.*X(T1_comps_10,:))+...
    age_scale00*sum(Aformer.*X(T1_comps_01,:))+age_scale10*sum(Apwid.*X(T1_comps_11,:));
t0=age_scale00*sum(Aformer.*X(T0_comps_00,:))+age_scale10*sum(Apwid.*X(T0_comps_10,:))+...
    age_scale00*sum(Aformer.*X(T0_comps_01,:))+age_scale10*sum(Apwid.*X(T0_comps_11,:));

treat_tot=s0+s1+s2+s3+s4+a0+f0+f1+f2+f3+f4+hcc+dc+lt1+lt2+t0+t1+t2+t3+t4;
figure;area(t1out,[lt' hcc' dc' f4' f3' f2' f1' f0'])
xticks(0:5:81)
set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
legend({'lt' 'hcc' 'dc' 'f4' 'f3' 'f2' 'f1' 'f0'},'location','best');
grid()
title('Number of PWID at each stage of chronic HCV')
set(gca,'xlim',[67 82])

Tqualy_vals00=zeros(20,1); % This conatins qualy for all compartments 00 (former, not failed)
Tqualy_vals01=zeros(20,1); % This conatins qualy for all compartments 01 (former, failed)
Tqualy_vals10=zeros(20,1); % This conatins qualy for all compartments 10 (current, not failed)
Tqualy_vals11=zeros(20,1); % This conatins qualy for all compartments 11 (current, failed)
XX=0;
for i = 1 : 20
    Xtest00=age_scale00*sum(Aformer.*X(Tcompvec00(i,:),:));
    Xtest01=age_scale00*sum(Aformer.*X(Tcompvec01(i,:),:));
    Xtest10=age_scale10*sum(Apwid.*X(Tcompvec10(i,:),:));
    Xtest11=age_scale10*sum(Apwid.*X(Tcompvec11(i,:),:));
    %     Xtest00=scale00*sum(X(Tcompvec00(i,:),:));
    %     Xtest01=scale00*sum(X(Tcompvec01(i,:),:));
    %     Xtest10=scale10*sum(X(Tcompvec10(i,:),:));
    %     Xtest11=scale10*sum(X(Tcompvec11(i,:),:));
    XX=XX+Xtest00+Xtest01+Xtest10+Xtest11;
    Tqualy_vals00(i)=calc_qalyv2(Xtest00(time2016:end),2016,2031,1,Tutility_vec(i),3); % 2000 to 2030 in steps of 1 month
    Tqualy_vals01(i)=calc_qalyv2(Xtest01(time2016:end),2016,2031,1,Tutility_vec(i),3); % 2000 to 2030 in steps of 1 month
    Tqualy_vals10(i)=calc_qalyv2(Xtest10(time2016:end),2016,2031,1,Tutility_vec(i),3); % 2000 to 2030 in steps of 1 month
    Tqualy_vals11(i)=calc_qalyv2(Xtest11(time2016:end),2016,2031,1,Tutility_vec(i),3); % 2000 to 2030 in steps of 1 month
end
Ttotal_qual_current_scaled = sum(Tqualy_vals00)+sum(Tqualy_vals01)+sum(Tqualy_vals10)+sum(Tqualy_vals11);

% s compartments
s4=age_scale00*sum(Aformer.*X(S4_comps_00,:))+age_scale10*sum(Apwid.*X(S4_comps_10,:))+...
    age_scale00*sum(Aformer.*X(S4_comps_01,:))+age_scale10*sum(Apwid.*X(S4_comps_11,:));
s3=age_scale00*sum(Aformer.*X(S3_comps_00,:))+age_scale10*sum(Apwid.*X(S3_comps_10,:))+...
    age_scale00*sum(Aformer.*X(S3_comps_01,:))+age_scale10*sum(Apwid.*X(S3_comps_11,:));
s2=age_scale00*sum(Aformer.*X(S2_comps_00,:))+age_scale10*sum(Apwid.*X(S2_comps_10,:))+...
    age_scale00*sum(Aformer.*X(S2_comps_01,:))+age_scale10*sum(Apwid.*X(S2_comps_11,:));
s1=age_scale00*sum(Aformer.*X(S1_comps_00,:))+age_scale10*sum(Apwid.*X(S1_comps_10,:))+...
    age_scale00*sum(Aformer.*X(S1_comps_01,:))+age_scale10*sum(Apwid.*X(S1_comps_11,:));
s0=age_scale00*sum(Aformer.*X(S0_comps_00,:))+age_scale10*sum(Apwid.*X(S0_comps_10,:))+...
    age_scale00*sum(Aformer.*X(S0_comps_01,:))+age_scale10*sum(Apwid.*X(S0_comps_11,:));
a0=age_scale00*sum(Aformer.*X(A_comps_00,:))+age_scale10*sum(Apwid.*X(A_comps_10,:))+...
    age_scale00*sum(Aformer.*X(A_comps_01,:))+age_scale10*sum(Apwid.*X(A_comps_11,:));

t4=age_scale00*sum(Aformer.*X(T4_comps_00,:))+age_scale10*sum(Apwid.*X(T4_comps_10,:))+...
    age_scale00*sum(Aformer.*X(T4_comps_01,:))+age_scale10*sum(Apwid.*X(T4_comps_11,:));
t3=age_scale00*sum(Aformer.*X(T3_comps_00,:))+age_scale10*sum(Apwid.*X(T3_comps_10,:))+...
    age_scale00*sum(Aformer.*X(T3_comps_01,:))+age_scale10*sum(Apwid.*X(T3_comps_11,:));
t2=age_scale00*sum(Aformer.*X(T2_comps_00,:))+age_scale10*sum(Apwid.*X(T2_comps_10,:))+...
    age_scale00*sum(Aformer.*X(T2_comps_01,:))+age_scale10*sum(Apwid.*X(T2_comps_11,:));
t1=age_scale00*sum(Aformer.*X(T1_comps_00,:))+age_scale10*sum(Apwid.*X(T1_comps_10,:))+...
    age_scale00*sum(Aformer.*X(T1_comps_01,:))+age_scale10*sum(Apwid.*X(T1_comps_11,:));
t0=age_scale00*sum(Aformer.*X(T0_comps_00,:))+age_scale10*sum(Apwid.*X(T0_comps_10,:))+...
    age_scale00*sum(Aformer.*X(T0_comps_01,:))+age_scale10*sum(Apwid.*X(T0_comps_11,:));
CCtreat=[sum(s0(67:end)) ,sum(s1(67:end)) ,sum(s2(67:end)) ,sum(s3(67:end)) ,sum(s4(67:end)),...
    sum(a0(67:end)),...
    sum(f0(67:end)) ,sum(f1(67:end)) ,sum(f2(67:end)) ,sum(f3(67:end)) ,sum(f4(67:end)),...
    sum(dc(67:end)) ,sum(hcc(67:end)), sum(lt(67:end)),...
    sum(t0(67:end)) ,sum(t1(67:end)) ,sum(t2(67:end)) ,sum(t3(67:end)) ,sum(t4(67:end))];
totinmodeltreat=s0(67:end)+s1(67:end)+s2(67:end)+s3(67:end)+s4(67:end)+...
    a0(67:end)+...
    f0(67:end)+f1(67:end)+f2(67:end)+f3(67:end)+f4(67:end)+...
    dc(67:end)+hcc(67:end)+lt(67:end)+...
    t0(67:end)+t1(67:end)+t2(67:end)+t3(67:end)+t4(67:end);


CCtreat2=[sum(s0(1:65)) ,sum(s1(1:65)) ,sum(s2(1:65)) ,sum(s3(1:65)) ,sum(s4(1:65)),...
    sum(a0(1:65)),...
    sum(f0(1:65)) ,sum(f1(1:65)) ,sum(f2(1:65)) ,sum(f3(1:65)) ,sum(f4(1:65)),...
    sum(dc(1:65)) ,sum(hcc(1:65)), sum(lt(1:65)),...
    0 0 0 0 0];
% incidence for treat
T00_scaled=(1-delta)*rate*age_scale00*sum(Aformer.*X(A_comps_00,:));
T10_scaled=(1-delta)*rate*age_scale10*sum(Apwid.*X(A_comps_10,:));
T01_scaled=(1-delta)*rate*age_scale00*sum(Aformer.*X(A_comps_01,:));
T11_scaled=(1-delta)*rate*age_scale10*sum(Apwid.*X(A_comps_11,:));
% integrate
inci_acute00_scaled=calc_qalyv2(T00_scaled(time2016:82),2016,2031,1,1,0); % 2000 to 2030 in steps of 1 month
inci_acute10_scaled=calc_qalyv2(T10_scaled(time2016:82),2016,2031,1,1,0); % 2000 to 2030 in steps of 1 month
inci_acute01_scaled=calc_qalyv2(T01_scaled(time2016:82),2016,2031,1,1,0); % 2000 to 2030 in steps of 1 month
inci_acute11_scaled=calc_qalyv2(T11_scaled(time2016:82),2016,2031,1,1,0); % 2000 to 2030 in steps of 1 month

yearly_diagnosed_cases_in_chronic_stages_scaled_treat=inci_acute00_scaled+inci_acute10_scaled+...
    inci_acute01_scaled+inci_acute11_scaled;
[~,Ttot_liver_deaths_scaled,XTdeaths]=fatalQALYS(age_scale00,age_scale10,age_weights,age_weights_former,XT,time2016,testcompvec10,testcompvec00,death_vec);

Qtreat = Ttotal_qual_current_scaled;%+Tfatal_QALY_scaled;
Qbase = total_qual_current_scaled;%+fatal_QALY_scaled;

% costs


%figure;plot(2016:2030,totaltestsF4)
% ********************************************************
%   treat costs need to check if all treats are used
%   Now fixed in 
%   odeeq_v5T_2p_scaled_death 
%   also in trt_cost
%   ndash=min(ndash, age_scale00*sum(Aformer.*X(jformer))); % can't treat more patinets than there
% ********************************************************
Ttotal_costs_scaled=trt_cost(age_weights,age_weights_former,age_scale10,age_scale00,XT,timey,n0,n1,y1,pcom,alpha,...
    mort_former,mort_current,...
    time2016,costanual,costoneoff,parm_current,...
    DC_comps_00,DC_comps_10,HCC_comps_00,HCC_comps_10,...
    DC_comps_01,DC_comps_11,HCC_comps_01,HCC_comps_11,...
    A_comps_00,A_comps_10,...
    A_comps_01,A_comps_11,...
    F0_comps_10,F1_comps_10,F2_comps_10,...
    F0_comps_11,F1_comps_11,F2_comps_11,...
    F0_comps_00,F1_comps_00,F2_comps_00,...
    F0_comps_01,F1_comps_01,F2_comps_01,...
    F3_comps_00,F3_comps_10,F3_comps_01,F3_comps_11,...
    F4_comps_00,F4_comps_10,F4_comps_01,F4_comps_11);

ICER_scaled = (Ttotal_costs_scaled-total_costs_scaled)/(Qtreat-Qbase);
import java.text.*
v=DecimalFormat;
outy{1,1}='No treatment';
outy{1,2}=['$' char(v.format(round(total_costs_scaled)))];
outy{1,3}=char(v.format(round(Qbase)));
outy{1,4}=char(v.format(round(tot_liver_deaths_scaled)));
outy{1,5}='';

outy{2,1}='Composite treatment';
outy{2,2}=['$' char(v.format(round(Ttotal_costs_scaled)))];
outy{2,3}=char(v.format(round(Qtreat)));
outy{2,4}=char(v.format(round(Ttot_liver_deaths_scaled)));
outy{2,5}=['$' char(v.format(round(ICER_scaled)))];

cell2table(outy,'VariableNames',{'Intervention' 'TotalCost','TotalQALYS' 'LiverRelatedDeaths' 'ICER'})

yinci_PWID=plotintervention(XT,delta,rate,age_weights,age_weights_former,age_scale00,age_scale10,time2016,1,...
    A_comps_00,A_comps_10,A_comps_01,A_comps_11,XTdeaths);

[lt,dc,hcc,f4,f3,f2,f1,f0,s4,s3,s2,s1,s0,t4,t3,t2,t1,t0,a0...
    ]=inallcomps(age_scale00,age_scale10,age_weights,age_weights_former,XT,...
    LT1_comps_00,LT1_comps_10,LT1_comps_01,LT1_comps_11,...
    LT2_comps_00,LT2_comps_10,LT2_comps_01,LT2_comps_11,...
    DC_comps_00,DC_comps_10,DC_comps_01,DC_comps_11,...
    HCC_comps_00,HCC_comps_10,HCC_comps_01,HCC_comps_11,...
    F4_comps_00,F4_comps_10,F4_comps_01,F4_comps_11,...
    F3_comps_00,F3_comps_10,F3_comps_01,F3_comps_11,...
    F2_comps_00,F2_comps_10,F2_comps_01,F2_comps_11,...
    F1_comps_00,F1_comps_10,F1_comps_01,F1_comps_11,...
    F0_comps_00,F0_comps_10,F0_comps_01,F0_comps_11,...
    S4_comps_00,S4_comps_10,S4_comps_01,S4_comps_11,...
    S3_comps_00,S3_comps_10,S3_comps_01,S3_comps_11,...
    S2_comps_00,S2_comps_10,S2_comps_01,S2_comps_11,...
    S1_comps_00,S1_comps_10,S1_comps_01,S1_comps_11,...
    S0_comps_00,S0_comps_10,S0_comps_01,S0_comps_11,...
    T4_comps_00,T4_comps_10,T4_comps_01,T4_comps_11,...
    T3_comps_00,T3_comps_10,T3_comps_01,T3_comps_11,...
    T2_comps_00,T2_comps_10,T2_comps_01,T2_comps_11,...
    T1_comps_00,T1_comps_10,T1_comps_01,T1_comps_11,...
    T0_comps_00,T0_comps_10,T0_comps_01,T0_comps_11,...
    A_comps_00,A_comps_10,A_comps_01,A_comps_11);

figure;area(timey,lt+dc+hcc+f4+f3+f2+f1+f0+s4+s3+s2+s1+s0+t4+t3+t2+t1+t0+a0)
xticks(0:5:81)
set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
grid()
title('Number of PWID in all compartments - treatment model')
set(gca,'xlim',[67 82])



figure;area(timey,ltb+dcb+hccb+f4b+f3b+f2b+f1b+f0b+s4b+s3b+s2b+s1b+s0b+t4b+t3b+t2b+t1b+t0b+a0b)
xticks(0:5:81)
set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
grid()
title('Number of PWID in all compartments - base model')
set(gca,'xlim',[67 82])
figure;plot(timey,(ltb+dcb+hccb+f4b+f3b+f2b+f1b+f0b+s4b+s3b+s2b+s1b+s0b+t4b+t3b+t2b+t1b+t0b+a0b))
hold on;
plot(timey,(lt+dc+hcc+f4+f3+f2+f1+f0+s4+s3+s2+s1+s0+t4+t3+t2+t1+t0+a0),'r-')
xticks(0:5:81)
set(gca,'xticklabel',1950:5:2030,'xticklabelrotation',90)
grid()
title('Number of PWID in all compartments - both model')
set(gca,'xlim',[67 82])
df=45

function [Q_tot_liver_deaths_scaled,tot_liver_deaths_scaled,Xdeaths_scaled]=fatalQALYS(age_scale00,age_scale10,age_weights,age_weights_former,X,time2016,testcompvec10,testcompvec00,death_vec)
% tot liver deaths with age weighting
testcompvec01_failed=[192:20:360;193:20:360;194:20:360;195:20:360];
testcompvec11_failed=[552:20:720;553:20:720;554:20:720;555:20:720];
Apwid=repmat(repelem(age_weights,1),1,82);
Aformer=repmat(repelem(age_weights_former,1),1,82);
death_vals10_scaled=zeros(4,1);
death_vals00_scaled=zeros(4,1);
death_vals11_scaled=zeros(4,1);
death_vals01_scaled=zeros(4,1);
Xdeaths_scaled = zeros(1,2030-2016+1);
m=size(X);
Deathmat_scaled=zeros(9,m(2));
for i = 1:4
    Deathmat_scaled=Deathmat_scaled + (age_scale10*Apwid.*X(testcompvec10(i,:),:)+age_scale00*Aformer.*X(testcompvec00(i,:),:)*death_vec(i)+...
        age_scale10*Apwid.*X(testcompvec11_failed(i,:),:)+age_scale00*Aformer.*X(testcompvec01_failed(i,:),:))*death_vec(i);
    Xtest=age_scale10*sum(Apwid.*X(testcompvec10(i,:),:));
    Ytest=age_scale00*sum(Aformer.*X(testcompvec00(i,:),:));
    X1test=age_scale10*sum(Apwid.*X(testcompvec11_failed(i,:),:));
    Y1test=age_scale00*sum(Aformer.*X(testcompvec01_failed(i,:),:));
    Xdeaths_scaled=Xdeaths_scaled+...
        just_integrate(Xtest(time2016:82),2016,2031,1)*death_vec(i)+...
        just_integrate(Ytest(time2016:82),2016,2031,1)*death_vec(i)+...
        just_integrate(X1test(time2016:82),2016,2031,1)*death_vec(i)+...
        just_integrate(Y1test(time2016:82),2016,2031,1)*death_vec(i);
    death_vals10_scaled(i)=calc_qalyv2(Xtest(time2016:82),2016,2031,1,1,0); % 2000 to 2030 in steps of 1 month
    death_vals00_scaled(i)=calc_qalyv2(Ytest(time2016:82),2016,2031,1,1,0); % 2000 to 2030 in steps of 1 month
    death_vals11_scaled(i)=calc_qalyv2(X1test(time2016:82),2016,2030,1,1,0); % 2000 to 2030 in steps of 1 month
    death_vals01_scaled(i)=calc_qalyv2(Y1test(time2016:82),2016,2030,1,1,0); % 2000 to 2030 in steps of 1 month
    
end
tot_liver_deaths_scaled=sum(death_vals00_scaled.*death_vec')+sum(death_vals10_scaled.*death_vec')+...
    sum(death_vals01_scaled.*death_vec')+sum(death_vals11_scaled.*death_vec');
t=1950:2031;
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


function total_costs_scaled=base_cost(age_weights,age_weights_former,age_scale10,age_scale00,X,...
    time2016,costanual,costoneoff,parm_current,...
    DC_comps_00,DC_comps_10,HCC_comps_00,HCC_comps_10,...
    A_comps_00,A_comps_10,...
    F0_comps_10,F1_comps_10,F2_comps_10,...
    F0_comps_00,F1_comps_00,F2_comps_00,...
    F3_comps_00,F3_comps_10,...
    F4_comps_00,F4_comps_10)
compvec10= [F0_comps_10;F1_comps_10;F2_comps_10;F3_comps_10;F4_comps_10;HCC_comps_10;DC_comps_10];
compvec00= [F0_comps_00;F1_comps_00;F2_comps_00;F3_comps_00;F4_comps_00;HCC_comps_00;DC_comps_00];
Apwid=repmat(repelem(age_weights,1),1,82);
Aformer=repmat(repelem(age_weights_former,1),1,82);
cost_vals10_scaled=zeros(7,1); % This conatins qualy for each chronic department
cost_vals00_scaled=zeros(7,1); % This conatins qualy for each chronic department
for i = 1 : 7
    Xtest=age_scale10*sum(Apwid.*X(compvec10(i,:),:));
    Ytest=age_scale00*sum(Aformer.*X(compvec00(i,:),:));
    %testy=testy+Xtest(time2015);
    cost_vals10_scaled(i)=calc_qalyv2(Xtest(time2016:82),2016,2031,1,costanual(i),3); % 2000 to 2030 in steps of 1 month
    cost_vals00_scaled(i)=calc_qalyv2(Ytest(time2016:82),2016,2031,1,costanual(i),3); % 2000 to 2030 in steps of 1 month
end
total_annual_costs_scaled = sum(cost_vals10_scaled)+sum(cost_vals00_scaled);

% total liver costs
% look at number of liver transplants
% DC to LT the rate is parm_names(13)

r_DCLT=parm_current(13);
T00_scaled=r_DCLT*age_scale00*sum(Aformer.*X(DC_comps_00,:));
T10_scaled=r_DCLT*age_scale10*sum(Apwid.*X(DC_comps_10,:));
% integrate
DCLT00_scaled=calc_qalyv2(T00_scaled(time2016:82),2016,2031,1,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
DCLT10_scaled=calc_qalyv2(T10_scaled(time2016:82),2016,2031,1,costoneoff(3),3); % 2000 to 2030 in steps of 1 month

% HCC to LT the rate is parm_names(15)
r_HCCLT=parm_current(15);
T00_scaled=r_HCCLT*age_scale00*sum(Aformer.*X(HCC_comps_00,:));
T10_scaled=r_HCCLT*age_scale10*sum(Apwid.*X(HCC_comps_10,:));
% integrate
HCCLT00_scaled=calc_qalyv2(T00_scaled(time2016:82),2016,2031,1,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
HCCLT10_scaled=calc_qalyv2(T10_scaled(time2016:82),2016,2031,1,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
total_cost_liver_transplants_scaled = DCLT00_scaled + DCLT10_scaled + HCCLT00_scaled + HCCLT10_scaled;

% one off costs for diagnosis of incident chronic cases
rate=parm_current(2);
delta=parm_current(1);
T00_scaled=(1-delta)*rate*age_scale00*sum(Aformer.*X(A_comps_00,:));
T10_scaled=(1-delta)*rate*age_scale10*sum(Apwid.*X(A_comps_10,:));
% integrate
inci_acute00_scaled=calc_qalyv2(T00_scaled(time2016:82),2016,2031,1,costoneoff(4),3); % 2000 to 2030 in steps of 1 month
inci_acute10_scaled=calc_qalyv2(T10_scaled(time2016:82),2016,2031,1,costoneoff(4),3); % 2000 to 2030 in steps of 1 month
yearly_diagnosed_cases_in_chronic_stages_scaled=inci_acute00_scaled+inci_acute10_scaled;

% one off costs for transition incidence
% 'r_F3F4','r_F4DC','r_F4HCC'
% F3 to F4
r_F3F4=parm_current(9);
% F4 to HCC
% as above scaled
T00_scaled=r_F3F4*age_scale00*sum(Aformer.*X(F3_comps_00,:));
T10_scaled=r_F3F4*age_scale10*sum(Apwid.*X(F3_comps_10,:));
% integrate
F3F400_scaled=calc_qalyv2(T00_scaled(time2016:82),2016,2031,1,costoneoff(1),3); % 2000 to 2030 in steps of 1 month
F3F410_scaled=calc_qalyv2(T10_scaled(time2016:82),2016,2031,1,costoneoff(1),3); % 2000 to 2030 in steps of 1 month
% F4 to HCC
r_F4HCC=parm_current(11);
T00_scaled=r_F4HCC*age_scale00*sum(Aformer.*X(F4_comps_00,:));
T10_scaled=r_F4HCC*age_scale10*sum(Apwid.*X(F4_comps_10,:));
% integrate
F4HCC00_scaled=calc_qalyv2(T00_scaled(time2016:82),2016,2031,1,costoneoff(2),3); % 2000 to 2030 in steps of 1 month
F4HCC10_scaled=calc_qalyv2(T10_scaled(time2016:82),2016,2031,1,costoneoff(2),3); % 2000 to 2030 in steps of 1 month
total_transistion_costs_scaled = F3F400_scaled + F3F410_scaled + F4HCC00_scaled + F4HCC10_scaled;

total_costs_scaled = total_annual_costs_scaled+total_cost_liver_transplants_scaled+ ...
    yearly_diagnosed_cases_in_chronic_stages_scaled+total_transistion_costs_scaled;

function Ttotal_costs_scaled=trt_cost(age_weights,age_weights_former,age_scale10,age_scale00,XT,t1,n0,n1,y1,pcom,alpha,...
    mort_former,mort_current,...
    time2016,costanual,costoneoff,parm_current,...
    DC_comps_00,DC_comps_10,HCC_comps_00,HCC_comps_10,...
    DC_comps_01,DC_comps_11,HCC_comps_01,HCC_comps_11,...
    A_comps_00,A_comps_10,...
    A_comps_01,A_comps_11,...
    F0_comps_10,F1_comps_10,F2_comps_10,...
    F0_comps_11,F1_comps_11,F2_comps_11,...
    F0_comps_00,F1_comps_00,F2_comps_00,...
    F0_comps_01,F1_comps_01,F2_comps_01,...
    F3_comps_00,F3_comps_10,F3_comps_01,F3_comps_11,...
    F4_comps_00,F4_comps_10,F4_comps_01,F4_comps_11)
compvec10= [F0_comps_10;F1_comps_10;F2_comps_10;F3_comps_10;F4_comps_10;HCC_comps_10;DC_comps_10];
compvec00= [F0_comps_00;F1_comps_00;F2_comps_00;F3_comps_00;F4_comps_00;HCC_comps_00;DC_comps_00];
compvec11= [F0_comps_11;F1_comps_11;F2_comps_11;F3_comps_11;F4_comps_11;HCC_comps_11;DC_comps_11];
compvec01= [F0_comps_01;F1_comps_01;F2_comps_01;F3_comps_01;F4_comps_01;HCC_comps_01;DC_comps_01];
Apwid=repmat(repelem(age_weights,1),1,82);
Aformer=repmat(repelem(age_weights_former,1),1,82);

Tcost_vals00_scaled=zeros(7,1); % This conatins qualy for 00
Tcost_vals01_scaled=zeros(7,1); % This conatins qualy for 01
Tcost_vals10_scaled=zeros(7,1); % This conatins qualy for 10
Tcost_vals11_scaled=zeros(7,1); % This conatins qualy for 11

for i = 1 : 7
    Xtest00_scaled=age_scale00*sum(Aformer.*XT(compvec00(i,:),:));
    Xtest01_scaled=age_scale00*sum(Aformer.*XT(compvec01(i,:),:));
    Xtest10_scaled=age_scale10*sum(Apwid.*XT(compvec10(i,:),:));
    Xtest11_scaled=age_scale10*sum(Apwid.*XT(compvec11(i,:),:));
    Tcost_vals00_scaled(i)=calc_qalyv2(Xtest00_scaled(time2016:82),2016,2031,1,costanual(i),3); % 2000 to 2030 in steps of 1 month
    Tcost_vals01_scaled(i)=calc_qalyv2(Xtest01_scaled(time2016:82),2016,2031,1,costanual(i),3); % 2000 to 2030 in steps of 1 month
    Tcost_vals10_scaled(i)=calc_qalyv2(Xtest10_scaled(time2016:82),2016,2031,1,costanual(i),3); % 2000 to 2030 in steps of 1 month
    Tcost_vals11_scaled(i)=calc_qalyv2(Xtest11_scaled(time2016:82),2016,2031,1,costanual(i),3); % 2000 to 2030 in steps of 1 month
end
Ttotal_annual_costs_scaled = sum(Tcost_vals00_scaled)+sum(Tcost_vals01_scaled)+sum(Tcost_vals10_scaled)+sum(Tcost_vals11_scaled);

% total liver costs
% look at number of liver transplants
% DC to LT the rate is parm_names(13)

r_DCLT=parm_current(13);
r_HCCLT=parm_current(15);
T00_scaled=r_DCLT*age_scale00*sum(Aformer.*XT(DC_comps_00,:));
T01_scaled=r_DCLT*age_scale00*sum(Aformer.*XT(DC_comps_01,:));
T10_scaled=r_DCLT*age_scale10*sum(Apwid.*XT(DC_comps_10,:));
T11_scaled=r_DCLT*age_scale10*sum(Apwid.*XT(DC_comps_11,:));
% integrate
DCLT00_scaled=calc_qalyv2(T00_scaled(time2016:82),2016,2031,1,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
DCLT01_scaled=calc_qalyv2(T01_scaled(time2016:82),2016,2031,1,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
DCLT10_scaled=calc_qalyv2(T10_scaled(time2016:82),2016,2031,1,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
DCLT11_scaled=calc_qalyv2(T11_scaled(time2016:82),2016,2031,1,costoneoff(3),3); % 2000 to 2030 in steps of 1 month

% HCC to LT the rate is parm_names(15)
T00_scaled=r_HCCLT*age_scale00*sum(Aformer.*XT(HCC_comps_00,:));
T01_scaled=r_HCCLT*age_scale00*sum(Aformer.*XT(HCC_comps_01,:));
T10_scaled=r_HCCLT*age_scale10*sum(Apwid.*XT(HCC_comps_10,:));
T11_scaled=r_HCCLT*age_scale10*sum(Apwid.*XT(HCC_comps_11,:));
% integrate
HCCLT00_scaled=calc_qalyv2(T00_scaled(time2016:82),2016,2031,1,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
HCCLT01_scaled=calc_qalyv2(T01_scaled(time2016:82),2016,2031,1,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
HCCLT10_scaled=calc_qalyv2(T10_scaled(time2016:82),2016,2031,1,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
HCCLT11_scaled=calc_qalyv2(T11_scaled(time2016:82),2016,2031,1,costoneoff(3),3); % 2000 to 2030 in steps of 1 month
Ttotal_cost_liver_transplants_scaled = DCLT00_scaled + DCLT10_scaled + HCCLT00_scaled + HCCLT10_scaled+...
    DCLT01_scaled + DCLT11_scaled + HCCLT01_scaled + HCCLT11_scaled;

% one off costs for diagnosis of incident chronic cases
rate=parm_current(2);
delta=parm_current(1);
T00_scaled=(1-delta)*rate*age_scale00*sum(Aformer.*XT(A_comps_00,:));
T01_scaled=(1-delta)*rate*age_scale00*sum(Aformer.*XT(A_comps_01,:));
T10_scaled=(1-delta)*rate*age_scale10*sum(Apwid.*XT(A_comps_10,:));
T11_scaled=(1-delta)*rate*age_scale10*sum(Apwid.*XT(A_comps_11,:));
% integrate
inci_acute00_scaled=calc_qalyv2(T00_scaled(time2016:82),2016,2031,1,costoneoff(4),3); % 2000 to 2030 in steps of 1 month
inci_acute01_scaled=calc_qalyv2(T01_scaled(time2016:82),2016,2031,1,costoneoff(4),3); % 2000 to 2030 in steps of 1 month
inci_acute10_scaled=calc_qalyv2(T10_scaled(time2016:82),2016,2031,1,costoneoff(4),3); % 2000 to 2030 in steps of 1 month
inci_acute11_scaled=calc_qalyv2(T11_scaled(time2016:82),2016,2031,1,costoneoff(4),3); % 2000 to 2030 in steps of 1 month
Tyearly_diagnosed_cases_in_chronic_stages_scaled=inci_acute00_scaled+inci_acute01_scaled+inci_acute10_scaled+inci_acute11_scaled;

% one off costs for transition incidence
% 'r_F3F4','r_F4DC','r_F4HCC'
% F3 to F4
r_F3F4=parm_current(9);
T00_scaled=r_F3F4*age_scale00*sum(Aformer.*XT(F3_comps_00,:));
T01_scaled=r_F3F4*age_scale00*sum(Aformer.*XT(F3_comps_01,:));
T10_scaled=r_F3F4*age_scale10*sum(Apwid.*XT(F3_comps_10,:));
T11_scaled=r_F3F4*age_scale10*sum(Apwid.*XT(F3_comps_11,:));
% integrate
F3F400_scaled=calc_qalyv2(T00_scaled(time2016:82),2016,2031,1,costoneoff(1),3); % 2000 to 2030 in steps of 1 month
F3F401_scaled=calc_qalyv2(T01_scaled(time2016:82),2016,2031,1,costoneoff(1),3); % 2000 to 2030 in steps of 1 month
F3F410_scaled=calc_qalyv2(T10_scaled(time2016:82),2016,2031,1,costoneoff(1),3); % 2000 to 2030 in steps of 1 month
F3F411_scaled=calc_qalyv2(T11_scaled(time2016:82),2016,2031,1,costoneoff(1),3); % 2000 to 2030 in steps of 1 month

% F4 to HCC
r_F4HCC=parm_current(11);
T00_scaled=r_F4HCC*age_scale00*sum(Aformer.*XT(F4_comps_00,:));
T01_scaled=r_F4HCC*age_scale00*sum(Aformer.*XT(F4_comps_01,:));
T10_scaled=r_F4HCC*age_scale10*sum(Apwid.*XT(F4_comps_10,:));
T11_scaled=r_F4HCC*age_scale10*sum(Apwid.*XT(F4_comps_11,:));
% integrate
F4HCC00_scaled=calc_qalyv2(T00_scaled(time2016:82),2016,2031,1,costoneoff(2),3); % 2000 to 2030 in steps of 1 month
F4HCC01_scaled=calc_qalyv2(T01_scaled(time2016:82),2016,2031,1,costoneoff(2),3); % 2000 to 2030 in steps of 1 month
F4HCC10_scaled=calc_qalyv2(T10_scaled(time2016:82),2016,2031,1,costoneoff(2),3); % 2000 to 2030 in steps of 1 month
F4HCC11_scaled=calc_qalyv2(T11_scaled(time2016:82),2016,2031,1,costoneoff(2),3); % 2000 to 2030 in steps of 1 month
Ttotal_transistion_costs_scaled = F3F400_scaled + F3F401_scaled + F3F410_scaled + F3F411_scaled +...
    F4HCC00_scaled + F4HCC01_scaled + F4HCC10_scaled + F4HCC11_scaled ;

% costs of follow-up on F4
%
time2015_inODEunits=time2016; % This is the time when treatment starts (67 is 2016)
moved_to_T4_former=zeros(9,length(t1));
moved_to_T4_current=zeros(9,length(t1));
total_not_T4=zeros(9,length(t1));
total_T4=zeros(9,length(t1));
tot_treated=zeros(length(t1),5);
for i = 1 : length(t1)
    if i==70%time2016
        df=45;
    end
    N=XT(:,i);
    % note t1(i) needs 1 added to it as t1 is 0:81 so t1(67) is actually
    % 66
    [~,~,moved_to_T4_former(:,i),moved_to_T4_current(:,i),total_not_T4(:,i),total_T4(:,i),tot_treated(i,:)]=countTreats(n0,n1,y1,N,age_scale00,age_scale10,t1(i)+1,pcom,alpha,time2015_inODEunits,age_weights',age_weights_former');
end
%outyF4z=outyF4(:,3:4:4*length(outyt)); % former mort_current,mort_former
%outyF4z2=outyF4(:,4:4:4*length(outyt)); % former mort_current,mort_former

% already integrated - just sum see below
% totaltestsF4=calc_ongoingF4treats(t1,total_T4,2016,2031);
% totaltests_notF4=calc_ongoingF4treats(t1,total_not_T4,2016,2031);

totaltestsF4=sum(total_T4(:,67:end));
totaltests_notF4=sum(total_not_T4(:,67:end));
figure;plot(2016:2031,totaltestsF4+totaltests_notF4)

% should be the same
yy=sum(tot_treated');
hold on ;plot(2016:2031,yy(67:end),'r*')
title('test the total treats plot')
figure;pareto(yy(67:end))
figure;bar(2016:2031,cumsum(yy(67:end)))
% outyF4 in fours
% 1. F4former inci, 2. F4current, 3. inci total not F4, 4. total F4
%outyF4z=outyF4(:,1:4:4*length(outyt)); % former mort_current,mort_former
pyformer=calc_ongoingF4costs(t1,moved_to_T4_former,2016,2031,mort_former,3);
pycurrent=calc_ongoingF4costs(t1,moved_to_T4_current,2016,2031,mort_current,3);
Tongoingf4costs = 557.62*(pyformer+pycurrent);
discount_rate=0.03;
cost_discount = 1./((1+discount_rate).^(0:(length(totaltests_notF4)-1)));
Ttreatcost = sum((totaltests_notF4*35806.67 + totaltestsF4 * 62377.15).*cost_discount);

Ttotal_costs_scaled = Ttotal_annual_costs_scaled+Ttotal_cost_liver_transplants_scaled+ ...
    Tyearly_diagnosed_cases_in_chronic_stages_scaled+Ttotal_transistion_costs_scaled+...
    Tongoingf4costs+Ttreatcost;

function yinci_PWID=plotintervention(XT,delta,rate,age_weights,age_weights_former,age_scale00,age_scale10,time2016,ploty,...
    A_comps_00,A_comps_10,A_comps_01,A_comps_11,XTdeaths)
Apwid=repmat(repelem(age_weights,1),1,82);
Aformer=repmat(repelem(age_weights_former,1),1,82);
T00=(1-delta)*rate*age_scale00*sum(Aformer.*XT(A_comps_00,:));
T01=(1-delta)*rate*age_scale00*sum(Aformer.*XT(A_comps_01,:));
T11=(1-delta)*rate*age_scale10*sum(Apwid.*XT(A_comps_11,:));
T10=(1-delta)*rate*age_scale10*sum(Apwid.*XT(A_comps_10,:));
% integrate
yinci_PWID=just_integrate(T10(time2016:82),2016,2031,1)+...
    just_integrate(T00(time2016:82),2016,2031,1)+...
    just_integrate(T01(time2016:82),2016,2031,1)+...
    just_integrate(T11(time2016:82),2016,2031,1);
if ploty==1
    figure;plot(2016:2030,100*((yinci_PWID/yinci_PWID(1))),'r-','LineWidth',2);
    grid()
    title('% incidence change for two treatment regimes')
    
    if ploty==1
        figure;plot(2016:2030,100*((XTdeaths/XTdeaths(1))),'r-','LineWidth',2);
        grid()
        title('% deaths change for two treatment regimes')
    end
end
function [lt,dc,hcc,f4,f3,f2,f1,f0,s4,s3,s2,s1,s0,t4,t3,t2,t1,t0,a0...
    ]=inallcomps(age_scale00,age_scale10,age_weights,age_weights_former,X,...
    LT1_comps_00,LT1_comps_10,LT1_comps_01,LT1_comps_11,...
    LT2_comps_00,LT2_comps_10,LT2_comps_01,LT2_comps_11,...
    DC_comps_00,DC_comps_10,DC_comps_01,DC_comps_11,...
    HCC_comps_00,HCC_comps_10,HCC_comps_01,HCC_comps_11,...
    F4_comps_00,F4_comps_10,F4_comps_01,F4_comps_11,...
    F3_comps_00,F3_comps_10,F3_comps_01,F3_comps_11,...
    F2_comps_00,F2_comps_10,F2_comps_01,F2_comps_11,...
    F1_comps_00,F1_comps_10,F1_comps_01,F1_comps_11,...
    F0_comps_00,F0_comps_10,F0_comps_01,F0_comps_11,...
    S4_comps_00,S4_comps_10,S4_comps_01,S4_comps_11,...
    S3_comps_00,S3_comps_10,S3_comps_01,S3_comps_11,...
    S2_comps_00,S2_comps_10,S2_comps_01,S2_comps_11,...
    S1_comps_00,S1_comps_10,S1_comps_01,S1_comps_11,...
    S0_comps_00,S0_comps_10,S0_comps_01,S0_comps_11,...
    T4_comps_00,T4_comps_10,T4_comps_01,T4_comps_11,...
    T3_comps_00,T3_comps_10,T3_comps_01,T3_comps_11,...
    T2_comps_00,T2_comps_10,T2_comps_01,T2_comps_11,...
    T1_comps_00,T1_comps_10,T1_comps_01,T1_comps_11,...
    T0_comps_00,T0_comps_10,T0_comps_01,T0_comps_11,...
    A_comps_00,A_comps_10,A_comps_01,A_comps_11)
Apwid=repmat(age_weights,1,82);
Aformer=repmat(age_weights_former,1,82);
lt2=age_scale00*sum(Aformer.*X(LT2_comps_00,:))+age_scale10*sum(Apwid.*X(LT2_comps_10,:))+...
    age_scale00*sum(Aformer.*X(LT2_comps_01,:))+age_scale10*sum(Apwid.*X(LT2_comps_11,:));
lt1=age_scale00*sum(Aformer.*X(LT1_comps_00,:))+age_scale10*sum(Apwid.*X(LT1_comps_10,:))+...
    age_scale00*sum(Aformer.*X(LT1_comps_01,:))+age_scale10*sum(Apwid.*X(LT1_comps_11,:));
lt=lt1+lt2;
dc=age_scale00*sum(Aformer.*X(DC_comps_00,:))+age_scale10*sum(Apwid.*X(DC_comps_10,:))+...
    age_scale00*sum(Aformer.*X(DC_comps_01,:))+age_scale10*sum(Apwid.*X(DC_comps_11,:));
hcc=age_scale00*sum(Aformer.*X(HCC_comps_00,:))+age_scale10*sum(Apwid.*X(HCC_comps_10,:))+...
    age_scale00*sum(Aformer.*X(HCC_comps_01,:))+age_scale10*sum(Apwid.*X(HCC_comps_11,:));
f4=age_scale00*sum(Aformer.*X(F4_comps_00,:))+age_scale10*sum(Apwid.*X(F4_comps_10,:))+...
    age_scale00*sum(Aformer.*X(F4_comps_01,:))+age_scale10*sum(Apwid.*X(F4_comps_11,:));
f3=age_scale00*sum(Aformer.*X(F3_comps_00,:))+age_scale10*sum(Apwid.*X(F3_comps_10,:))+...
    age_scale00*sum(Aformer.*X(F3_comps_01,:))+age_scale10*sum(Apwid.*X(F3_comps_11,:));
f2=age_scale00*sum(Aformer.*X(F2_comps_00,:))+age_scale10*sum(Apwid.*X(F2_comps_10,:))+...
    age_scale00*sum(Aformer.*X(F2_comps_01,:))+age_scale10*sum(Apwid.*X(F2_comps_11,:));
f1=age_scale00*sum(Aformer.*X(F1_comps_00,:))+age_scale10*sum(Apwid.*X(F1_comps_10,:))+...
    age_scale00*sum(Aformer.*X(F1_comps_01,:))+age_scale10*sum(Apwid.*X(F1_comps_11,:));
f0=age_scale00*sum(Aformer.*X(F0_comps_00,:))+age_scale10*sum(Apwid.*X(F0_comps_10,:))+...
    age_scale00*sum(Aformer.*X(F0_comps_01,:))+age_scale10*sum(Apwid.*X(F0_comps_11,:));

% s compartments
s4=age_scale00*sum(Aformer.*X(S4_comps_00,:))+age_scale10*sum(Apwid.*X(S4_comps_10,:))+...
    age_scale00*sum(Aformer.*X(S4_comps_01,:))+age_scale10*sum(Apwid.*X(S4_comps_11,:));
s3=age_scale00*sum(Aformer.*X(S3_comps_00,:))+age_scale10*sum(Apwid.*X(S3_comps_10,:))+...
    age_scale00*sum(Aformer.*X(S3_comps_01,:))+age_scale10*sum(Apwid.*X(S3_comps_11,:));
s2=age_scale00*sum(Aformer.*X(S2_comps_00,:))+age_scale10*sum(Apwid.*X(S2_comps_10,:))+...
    age_scale00*sum(Aformer.*X(S2_comps_01,:))+age_scale10*sum(Apwid.*X(S2_comps_11,:));
s1=age_scale00*sum(Aformer.*X(S1_comps_00,:))+age_scale10*sum(Apwid.*X(S1_comps_10,:))+...
    age_scale00*sum(Aformer.*X(S1_comps_01,:))+age_scale10*sum(Apwid.*X(S1_comps_11,:));
s0=age_scale00*sum(Aformer.*X(S0_comps_00,:))+age_scale10*sum(Apwid.*X(S0_comps_10,:))+...
    age_scale00*sum(Aformer.*X(S0_comps_01,:))+age_scale10*sum(Apwid.*X(S0_comps_11,:));
a0=age_scale00*sum(Aformer.*X(A_comps_00,:))+age_scale10*sum(Apwid.*X(A_comps_10,:))+...
    age_scale00*sum(Aformer.*X(A_comps_01,:))+age_scale10*sum(Apwid.*X(A_comps_11,:));

t4=age_scale00*sum(Aformer.*X(T4_comps_00,:))+age_scale10*sum(Apwid.*X(T4_comps_10,:))+...
    age_scale00*sum(Aformer.*X(T4_comps_01,:))+age_scale10*sum(Apwid.*X(T4_comps_11,:));
t3=age_scale00*sum(Aformer.*X(T3_comps_00,:))+age_scale10*sum(Apwid.*X(T3_comps_10,:))+...
    age_scale00*sum(Aformer.*X(T3_comps_01,:))+age_scale10*sum(Apwid.*X(T3_comps_11,:));
t2=age_scale00*sum(Aformer.*X(T2_comps_00,:))+age_scale10*sum(Apwid.*X(T2_comps_10,:))+...
    age_scale00*sum(Aformer.*X(T2_comps_01,:))+age_scale10*sum(Apwid.*X(T2_comps_11,:));
t1=age_scale00*sum(Aformer.*X(T1_comps_00,:))+age_scale10*sum(Apwid.*X(T1_comps_10,:))+...
    age_scale00*sum(Aformer.*X(T1_comps_01,:))+age_scale10*sum(Apwid.*X(T1_comps_11,:));
t0=age_scale00*sum(Aformer.*X(T0_comps_00,:))+age_scale10*sum(Apwid.*X(T0_comps_10,:))+...
    age_scale00*sum(Aformer.*X(T0_comps_01,:))+age_scale10*sum(Apwid.*X(T0_comps_11,:));


function f=optimization_approach(x,N0,M,Mdash,mort_current,mort_former,extra_parms_vals,...
    phiminusvals, phiplussvals,age_matrix,...
    age_scale00,age_scale10,pcom,alpha,age_weights,age_weights_former,death_vec,...
    delta,rate,A_comps_00,A_comps_01,A_comps_10,A_comps_11,...
    testcompvec10,testcompvec00,y1)
% n0,n1,y1
n0=x(1);
n1=x(2);
%y1=x(3)
% starting values are 5000,5000,6
time2016=67;
[t1,x1t] = ode45(@odeeq_v5T_2p_scaled,0:81,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix,...
    age_scale00,age_scale10,n0,n1,pcom,alpha,67,y1,age_weights',age_weights_former'); % note 2015 marker is 65 years
XT=x1t';
Apwid=repmat(repelem(age_weights,1),1,82);
Aformer=repmat(repelem(age_weights_former,1),1,82);
T00=(1-delta)*rate*age_scale00*sum(Aformer.*XT(A_comps_00,:));
T01=(1-delta)*rate*age_scale00*sum(Aformer.*XT(A_comps_01,:));
T11=(1-delta)*rate*age_scale10*sum(Apwid.*XT(A_comps_11,:));
T10=(1-delta)*rate*age_scale10*sum(Apwid.*XT(A_comps_10,:));
% integrate
yinci_PWID=just_integrate(T10(time2016:82),2016,2031,1)+...
    just_integrate(T00(time2016:82),2016,2031,1)+...
    just_integrate(T01(time2016:82),2016,2031,1)+...
    just_integrate(T11(time2016:82),2016,2031,1);
[~,~,XTdeaths]=fatalQALYS(age_scale00,age_scale10,age_weights,age_weights_former,XT,time2016,testcompvec10,testcompvec00,death_vec);
r2inci=(yinci_PWID(end)/yinci_PWID(1)-0.2)^2;
d2inci=(XTdeaths(end)/XTdeaths(1)-0.35)^2;
f=r2inci+d2inci;

function f=optimization_approach_death(x,N0,M,Mdash,mort_current,mort_former,extra_parms_vals,...
    phiminusvals, phiplussvals,age_matrix,...
    age_scale00,age_scale10,pcom,alpha,age_weights,age_weights_former,death_vec,...
    n0,n1,y1,...
    delta,rate,A_comps_00,A_comps_01,A_comps_10,A_comps_11,...
    LT1_comps_00,LT1_comps_10,LT1_comps_01,LT1_comps_11,...
    LT2_comps_00,LT2_comps_10,LT2_comps_01,LT2_comps_11,...
    DC_comps_00,DC_comps_10,DC_comps_01,DC_comps_11,...
    HCC_comps_00,HCC_comps_10,HCC_comps_01,HCC_comps_11,...
    F4_comps_00,F4_comps_10,F4_comps_01,F4_comps_11,...
    F3_comps_00,F3_comps_10,F3_comps_01,F3_comps_11,...
    F2_comps_00,F2_comps_10,F2_comps_01,F2_comps_11,...
    F1_comps_00,F1_comps_10,F1_comps_01,F1_comps_11,...
    F0_comps_00,F0_comps_10,F0_comps_01,F0_comps_11,...
    S4_comps_00,S4_comps_10,S4_comps_01,S4_comps_11,...
    S3_comps_00,S3_comps_10,S3_comps_01,S3_comps_11,...
    S2_comps_00,S2_comps_10,S2_comps_01,S2_comps_11,...
    S1_comps_00,S1_comps_10,S1_comps_01,S1_comps_11,...
    S0_comps_00,S0_comps_10,S0_comps_01,S0_comps_11,...
    T4_comps_00,T4_comps_10,T4_comps_01,T4_comps_11,...
    T3_comps_00,T3_comps_10,T3_comps_01,T3_comps_11,...
    T2_comps_00,T2_comps_10,T2_comps_01,T2_comps_11,...
    T1_comps_00,T1_comps_10,T1_comps_01,T1_comps_11,...
    T0_comps_00,T0_comps_10,T0_comps_01,T0_comps_11,...
    testcompvec10,testcompvec00,totbase)
% n0,n1,y1
% n0=x(1);
% n1=x(2);
% y1=x(3);
death_val=x(1);
% starting values are 5000,5000,6
time2016=67;
[t1,x1t] = ode45(@odeeq_v5T_2p_scaled_death,0:81,N0,[],M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix,...
    age_scale00,age_scale10,n0,n1,pcom,alpha,67,y1,age_weights',age_weights_former',death_val); % note 2015 marker is 65 years
XT=x1t';
Apwid=repmat(repelem(age_weights,1),1,82);
Aformer=repmat(repelem(age_weights_former,1),1,82);

[lt,dc,hcc,f4,f3,f2,f1,f0,s4,s3,s2,s1,s0,t4,t3,t2,t1,t0,a0...
    ]=inallcomps(age_scale00,age_scale10,age_weights,age_weights_former,XT,...
    LT1_comps_00,LT1_comps_10,LT1_comps_01,LT1_comps_11,...
    LT2_comps_00,LT2_comps_10,LT2_comps_01,LT2_comps_11,...
    DC_comps_00,DC_comps_10,DC_comps_01,DC_comps_11,...
    HCC_comps_00,HCC_comps_10,HCC_comps_01,HCC_comps_11,...
    F4_comps_00,F4_comps_10,F4_comps_01,F4_comps_11,...
    F3_comps_00,F3_comps_10,F3_comps_01,F3_comps_11,...
    F2_comps_00,F2_comps_10,F2_comps_01,F2_comps_11,...
    F1_comps_00,F1_comps_10,F1_comps_01,F1_comps_11,...
    F0_comps_00,F0_comps_10,F0_comps_01,F0_comps_11,...
    S4_comps_00,S4_comps_10,S4_comps_01,S4_comps_11,...
    S3_comps_00,S3_comps_10,S3_comps_01,S3_comps_11,...
    S2_comps_00,S2_comps_10,S2_comps_01,S2_comps_11,...
    S1_comps_00,S1_comps_10,S1_comps_01,S1_comps_11,...
    S0_comps_00,S0_comps_10,S0_comps_01,S0_comps_11,...
    T4_comps_00,T4_comps_10,T4_comps_01,T4_comps_11,...
    T3_comps_00,T3_comps_10,T3_comps_01,T3_comps_11,...
    T2_comps_00,T2_comps_10,T2_comps_01,T2_comps_11,...
    T1_comps_00,T1_comps_10,T1_comps_01,T1_comps_11,...
    T0_comps_00,T0_comps_10,T0_comps_01,T0_comps_11,...
    A_comps_00,A_comps_10,A_comps_01,A_comps_11);
treatval=lt+dc+hcc+f4+f3+f2+f1+f0+s4+s3+s2+s1+s0+t4+t3+t2+t1+t0+a0;
% 76 is 2025
death_adjust=abs(trapz(totbase(67:end))-trapz(treatval(67:end)));
%death_adjust=sum((treatval(67:end)-totbase(67:end)).^2);%sum((treatval(end)-totbase(end)).^2)+sum((treatval(76)-totbase(76)).^2);
f=death_adjust;



function [phi,phidash,moved_to_T4_former,moved_to_T4_current,total_not_T4,total_T4,tot_treated]=countTreats(n0,n1,y1,X,age_scale00,age_scale10,t,pcom,alpha,time2015_inODEunits,pwid_weights,former_weights)
% This allows simultaneous applications of p=0 and p=1
% p=0 is PWID 
% p=1 is adanced for PWID and former
% y0 and y1 are the number of years p=0 and p = 1 are applied for
% y0 not used so dropped from argument list
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

% also outputs the number of subjects going to T4
% from line 1057
% extra outputs to check treat counts
n0=round(n0); % ofr the integer optimization
n1=round(n1);
% needed to count treats
sumy_split=zeros(4,1);
sumy=0;
ndash=0;

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

Apwid=[pwid_weights pwid_weights pwid_weights pwid_weights]';
Apwid_long=[pwid_weights pwid_weights pwid_weights pwid_weights pwid_weights pwid_weights pwid_weights]';
Aformer=[former_weights former_weights former_weights former_weights]';

Topwider=age_scale10*sum(Apwid_long.*X(c_10));
Adformer = age_scale00*sum(Aformer.*X(r_00)); % advanced stage totals
Adpwider = age_scale10*sum(Apwid.*X(r_10));
Adtotal=Adformer+Adpwider;

% p=0

fprop_p0 =  n0/Topwider;
gprop_p0 =  0;
hprop_p0 = n0/Topwider;

% p=1
if t<=y1
    fprop_p1=n1/Adtotal;
    gprop_p1=n1/Adtotal;
    hprop_p1 =  0;
    n=n0+n1;
else
    fprop_p1=0;
    gprop_p1=0;
    hprop_p1=0;
    n=n0;
end

phi=zeros(720,1);
phi_advanced_from_PWID_allocation=zeros(720,1);
phi_advanced_from_ADVA_allocation=zeros(720,1);
jAdpwider=[F3_pwider F4_pwider LT1_pwider LT2_pwider ]; % Advanced pwid
jMipwider=[F0_pwider F1_pwider F2_pwider]; % Mild pwid
jAdformer=[F3_former F4_former LT1_former LT2_former ]; % Advanced former

if t>=time2015_inODEunits
    if t<=y1
    df=45;
end
    % how any treated in F3,F4,LT1,LT2
    % f function
     phi(jMipwider)=min(hprop_p0*X(jMipwider),X(jMipwider));
     phi(jAdpwider)=min(fprop_p0*X(jAdpwider),X(jAdpwider))+min(fprop_p1*X(jAdpwider),X(jAdpwider));%min((fprop_p0+fprop_p1)*X(jAdpwider),X(jAdpwider));
     phi_advanced_from_PWID_allocation(jAdpwider)=min(fprop_p0*X(jAdpwider),X(jAdpwider));
     phi_advanced_from_ADVA_allocation(jAdpwider)=min(fprop_p1*X(jAdpwider),X(jAdpwider));
     phi(jAdformer)=min(gprop_p1*X(jAdformer),X(jAdformer));
   
    % if there are more treats than subjects then there will be surplus
    [sumy,sumy_split]=scale_sumphi_p_scaled_for_treats(phi_advanced_from_PWID_allocation,phi_advanced_from_ADVA_allocation,phi,age_scale00,age_scale10,jAdformer,jAdpwider,jMipwider,pwid_weights,former_weights);

    ndash=n-sumy;                        
    if ndash > 1e-5
        if t<=y1
            % allocation is over
            % advanced current
            % advanced former
            % mild current
            % distribute over mild former
            jformer=[F0_former F1_former F2_former ];
            Aformer=[former_weights former_weights former_weights]';
            scaledbot=max(1,age_scale00*sum(Aformer.*X(jformer))); % just stops NAN for empty compartments
            % extra will still be zero below as X(jformer) all zero
            ndash=min(ndash, age_scale00*sum(Aformer.*X(jformer))); % can't treat more patinets than there
            if scaledbot>0
            phi(jformer)=phi(jformer)+ndash*X(jformer)/scaledbot; % extra allocation
            end
        else % p=0
            %Treatments were allocated proportionally across liver disease stages, 
            %and if the number of treatments available was greater than the number 
            %of current PWID eligible for treatment, 
            %remaining treatments were allocated to former PWID, 
            %proportionally across disease stages.
            jformer=[F0_former F1_former F2_former F3_former F4_former LT1_former LT2_former];
            Aformer=[former_weights former_weights former_weights former_weights former_weights former_weights former_weights ]';
            scaledbot=max(1,age_scale00*sum(Aformer.*X(jformer)));
            ndash=min(ndash, age_scale00*sum(Aformer.*X(jformer))); % can't treat more patinets than there
             if scaledbot>0
             phi(jformer)=phi(jformer)+ndash*X(jformer)/scaledbot; % extra allocation
             end
        end
    end
end
tot_treated=sumy_split;
tot_treated(end+1)=ndash;
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
% test
%moved_to_T4=scale00*alpha*(phi(F3_former)+phi(F4_former)+phi(LT1_former)+phi(LT2_former))+...
%    scale10*alpha*pcom*(phi(F3_pwider)+phi(F4_pwider)+phi(LT1_pwider)+phi(LT2_pwider));
total_not_T4 = age_scale10*pwid_weights'.*(phi(F0_pwider)+phi(F1_pwider)+phi(F2_pwider)+phi(F3_pwider)+...
               phi(LT1_pwider)+phi(LT2_pwider))+...
               age_scale00*former_weights'.*(phi(F0_former)+phi(F1_former)+phi(F2_former)+phi(F3_former)+...
               phi(LT1_former)+phi(LT2_former));
total_T4=age_scale00*former_weights'.*phi(F4_former)+ age_scale10*pwid_weights'.*phi(F4_pwider);       

moved_to_T4_former=age_scale00*alpha*former_weights'.*phi(F4_former);
moved_to_T4_current=age_scale10*alpha*pcom*pwid_weights'.*phi(F4_pwider); % 9 by 1 for ages
function sumy=scale_sumphi_p_scaled(phi,age_scale00,age_scale10,jAdformer,jAdpwider,jMipwider,pwid_weights,former_weights)
% phi conatins the number of subjects treated
% it needs tobe scaled to compare to total treatments
% This is for any p
% jAdformer is the former advanced compartments
% jAdpwider is the current advanced compartments
% jMipwider is the current mild compartments
% sumy below is the f comps + h comps + g compartments
Apwid_A=[pwid_weights pwid_weights pwid_weights pwid_weights]';
Apwid_M=[pwid_weights pwid_weights pwid_weights]';
Aformer=[former_weights former_weights former_weights former_weights]';
sumy = age_scale10*sum(Apwid_A.*phi(jAdpwider))+...
       age_scale10*sum(Apwid_M.*phi(jMipwider))+...
       age_scale00*sum(Aformer.*phi(jAdformer));

   function [sumy,sumy_split]=scale_sumphi_p_scaled_for_treats(phi_advanced_from_PWID_allocation,phi_advanced_from_ADVA_allocation,phi,age_scale00,age_scale10,jAdformer,jAdpwider,jMipwider,pwid_weights,former_weights)
% phi conatins the number of subjects treated
% it needs tobe scaled to compare to total treatments
% This is for any p
% jAdformer is the former advanced compartments
% jAdpwider is the current advanced compartments
% jMipwider is the current mild compartments
% sumy below is the f comps + h comps + g compartments
Apwid_A=[pwid_weights pwid_weights pwid_weights pwid_weights]';
Apwid_M=[pwid_weights pwid_weights pwid_weights]';
Aformer=[former_weights former_weights former_weights former_weights]';
sumy = age_scale10*sum(Apwid_A.*phi(jAdpwider))+...
       age_scale10*sum(Apwid_M.*phi(jMipwider))+...
       age_scale00*sum(Aformer.*phi(jAdformer));
% below in order
% advanced PWID allocation from n0
% catch-up PWID allocation to advanced from n1
% mild PWID allocation from n0
% cath-up former allocation to F3,F4,LT1,LT2
sumy_split=zeros(4,1);
sumy_split(1)=age_scale10*sum(Apwid_A.*phi_advanced_from_PWID_allocation(jAdpwider));
sumy_split(2)=age_scale10*sum(Apwid_A.*phi_advanced_from_ADVA_allocation(jAdpwider));
sumy_split(3)=age_scale10*sum(Apwid_M.*phi(jMipwider));
sumy_split(4)=age_scale00*sum(Aformer.*phi(jAdformer));







