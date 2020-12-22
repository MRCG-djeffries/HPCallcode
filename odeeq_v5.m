function hcv=odeeq_v5(t,N,M,Mdash,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix)
 % N are the compartments - there are
 % 80 (0,0 0,1 1,0 1,1) by 9 ages 
 % 0,0 are rows 1 : 20
 % 0,1 are rows 21 : 40
 % 1,0 are rows 41 : 60
 % 1,1 are rows 61 : 80
 % i,j notation: i = 0 is former PWID and i = 1 is current PWID
 %             : j = 0 is never failed trt and j = 1 is failed trt
N(N<0)=0;
 hcv = zeros(20*9*4,1);
 
 Irows = [6,7,8,9,10,11,12,13,14,15]; % rows of infectious PWID
 death_vec= [-log(1-0.138) -log(1-0.605) -log(1-0.169) -log(1-0.034)]; %DC,HCC,LT1,LT2
 
 
 % no now changed
%  top = sum(N(Irows+40,:),'all') + sum(N(Irows+60,:),'all'); % 40 is j = 0, 60 is j = 1 and i = 1 for both
%  bot = sum(N(((1:20)+40),:),'all') + sum(N((1:20)+60),'all'); % 40 is j = 0, 60 is j = 1 and i = 1 for both
 top10_index = reshape((repmat(Irows,9,1)+(360:20:520)')',length(Irows)*9,1);
 top11_index = reshape((repmat(Irows,9,1)+(540:20:700)')',length(Irows)*9,1);
 top = sum(N(top10_index))+sum(N(top11_index)); % 40 is j = 0, 60 is j = 1 and i = 1 for both
 bot = sum(N(361:end)); % 40 is j = 0, 60 is j = 1 and i = 1 for both
 I = top/bot;
%  if t>=0 && t<50
%      scaleI = 1+1.5*t/50;
%  elseif t>=50 && t<=55
%       scaleI = 2.5-1.5*(t-50)/5;
%  else
%      scaleI=1;
%  end
%  

 if t>=0 && t<51
     scaleI = 1+1.5*t/51;
 elseif t>51 && t<=56
      scaleI = 2.5-1.5*(t-51)/5;
%  elseif t>56 && t<=66
%      scaleI=1;
  else
     scaleI=1;
 end
 
 phi = scaleI*extra_parms_vals(1)*I;
 dum = Mdash(5,5);
 Mdash(phiminusvals(1:end-1))=-phi;
 Mdash(phiminusvals(end))=-phi+dum; %dum = -r_svr4DC-rsvr4HCC
 Mdash(phiplussvals)=phi;
 
 %**********************************************************
 % Note phi is assumed to be zero for former
 %**********************************************************
 
 % former, failed, ages grps 1 to 9, i=0, j = 0
 % current, failed, ages grps 1 to 9, i=1, j = 0
 % former, never failed, ages grps 1 to 9, i=0, j = 1
 % current, never failed, ages grps 1 to 9, i=1, j = 1
 
 X00 = reshape(N(1:20*9),20,9);
 X01 = reshape(N(181:(180+20*9)),20,9);
 X10 = reshape(N(361:(360+20*9)),20,9);
 X11 = reshape(N(541:(540+20*9)),20,9);
 
 % with deaths set as zero - this should be closed
 % age movement also set as zero
 % there are treatment so no failed treatment
 % so only D10 and D00
 % only for the youngest age group since no tranbsistsion
 d00=M*X00-extra_parms_vals(2)*X00+extra_parms_vals(3)*X10-mort_former.*X00+(age_matrix*X00')';
 d01=M*X01-extra_parms_vals(2)*X01+extra_parms_vals(3)*X11-mort_former.*X01+(age_matrix*X01')';
 d10=Mdash*X10+extra_parms_vals(2)*X00-extra_parms_vals(3)*X10-mort_current.*X10+(age_matrix*X10')';
 d11=Mdash*X11+extra_parms_vals(2)*X01-extra_parms_vals(3)*X11-mort_current.*X11+(age_matrix*X11')';
 
%  d00=M*X00-extra_parms_vals(2)*X00+extra_parms_vals(3)*X10-mort_former.*X00+(age_matrix*X00')';
%  d01=M*X01-extra_parms_vals(2)*X01+extra_parms_vals(3)*X11-mort_former.*X01+(age_matrix*X01')';
%  d10=Mdash*X10+extra_parms_vals(2)*X00-extra_parms_vals(3)*X10-mort_current.*X10+(age_matrix*X10')';
%  d11=Mdash*X11+extra_parms_vals(2)*X01-extra_parms_vals(3)*X11-mort_current.*X11+(age_matrix*X11')';
d10(1,1) =d10(1,1) + sum(mort_former.*X00,'all') + sum(mort_former.*X01,'all') + sum(mort_current.*X10,'all') + sum(mort_current.*X11,'all')...
                   + death_vec(1)*sum(X00(12,:),'all')  + death_vec(2)*sum(X00(13,:),'all')  + death_vec(3)*sum(X00(14,:),'all')  + death_vec(4)*sum(X00(15,:),'all')...
                   + death_vec(1)*sum(X01(12,:),'all')  + death_vec(2)*sum(X01(13,:),'all')  + death_vec(3)*sum(X01(14,:),'all')  + death_vec(4)*sum(X01(15,:),'all')...
                   + death_vec(1)*sum(X10(12,:),'all')  + death_vec(2)*sum(X10(13,:),'all')  + death_vec(3)*sum(X10(14,:),'all')  + death_vec(4)*sum(X10(15,:),'all')...
                   + death_vec(1)*sum(X11(12,:),'all')  + death_vec(2)*sum(X11(13,:),'all')  + death_vec(3)*sum(X11(14,:),'all')  + death_vec(4)*sum(X11(15,:),'all');

hcv=[reshape(d00,180,1);reshape(d01,180,1);reshape(d10,180,1);reshape(d11,180,1)];

 
 
 
 
                          
 
 
 
 
 
 
 
 
 
 
 