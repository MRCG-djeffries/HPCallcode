function hcv=odeeq_v5T_2p_scaled_death_acute(t,N,Macute,Mdash_acute,mort_current,mort_former,extra_parms_vals,phiminusvals, phiplussvals,age_matrix,...
                      age_scale00,age_scale10,treatprob,pcom,alpha,time2015_inODEunits,y1,pwid_weights,former_weights,...
                      death_val,percent_treated)
 % N are the compartments - there are
 % 80 (0,0 0,1 1,0 1,1) by 9 ages 
 % 0,0 are rows 1 : 20
 % 0,1 are rows 21 : 40
 % 1,0 are rows 41 : 60
 % 1,1 are rows 61 : 80
 % i,j notation: i = 0 is former PWID and i = 1 is current PWID
 %             : j = 0 is never failed trt and j = 1 is failed trt
 % T version - this means the treatment version
 
 % This version inputs two treatments
 % n0 and n1 are the number of treatments
 % no is for p=0, i.e. current PWID which reduces incidence
 % n1 is for p=1, i.e. advanced PWID and former which reduces mortality
 % y1 is year number that this is applied for after 2015 (66) which for 
 % 5 years is 71
%global outyt
%global outyF4
% if sum(isnan(N))>0
%     df = 56
% end
 N(N<0)=0;  
 hcv = zeros(21*9*4,1);
 
 Irows = [6,7,8,9,10,11,12,13,14,15]; % rows of infectious PWID
 death_vec= [-log(1-0.138) -log(1-0.605) -log(1-0.169) -log(1-0.034)]; %DC,HCC,LT1,LT2
 
%  (360:20:520)')',length(Irows)*9,1);
%  top11_index = reshape((repmat(Irows,9,1)+(540:20:700)'
 % no now changed
%  top = sum(N(Irows+40,:),'all') + sum(N(Irows+60,:),'all'); % 40 is j = 0, 60 is j = 1 and i = 1 for both
%  bot = sum(N(((1:20)+40),:),'all') + sum(N((1:20)+60),'all'); % 40 is j = 0, 60 is j = 1 and i = 1 for both
 top10_index = reshape((repmat(Irows,9,1)+(378:21:546)')',length(Irows)*9,1);
 top11_index = reshape((repmat(Irows,9,1)+(567:21:735)')',length(Irows)*9,1);
 top = sum(N(top10_index))+sum(N(top11_index)); % 40 is j = 0, 60 is j = 1 and i = 1 for both
 bot = sum(N(379:end)); % 40 is j = 0, 60 is j = 1 and i = 1 for both
 I = top/bot;
%  if t>=0 && t<50
%      scaleI = 1+1.5*t/50;
%  elseif t>=50 && t<=55
%       scaleI = 2.5-1.5*(t-50)/5;
%  else
%      scaleI=1;
%  end

   if t>=0 && t<51
     scaleI = 1+1.5*t/51;
 elseif t>51 && t<=56
      scaleI = 2.5-1.5*(t-51)/5;
%  elseif t>56 && t<=66
%      scaleI=1;
  else
     scaleI=1;
   end
  phiminusvals=phiminusvals+[0 1 2 3 4]; % allows for 21 by 21 matrix column order
  phiplussvals=phiplussvals+[0 1 2 3 4];
 phi = scaleI*extra_parms_vals(1)*I;
 dum = Mdash_acute(5,5);
 Mdash_acute(phiminusvals(1:end-1))=-phi;
 Mdash_acute(phiminusvals(end))=-phi+dum; %dum = -r_svr4DC-rsvr4HCC
 Mdash_acute(phiplussvals)=phi;

 % former, failed, ages grps 1 to 9, i=0, j = 0
 % current, failed, ages grps 1 to 9, i=1, j = 0
 % former, never failed, ages grps 1 to 9, i=0, j = 1
 % current, never failed, ages grps 1 to 9, i=1, j = 1
 
 X00 = reshape(N(1:21*9),21,9);
 X01 = reshape(N(190:(189+21*9)),21,9);
 X10 = reshape(N(379:(378+21*9)),21,9);
 X11 = reshape(N(568:(567+21*9)),21,9);
 
 % with deaths set as zero - this should be closed
 % age movement also set as zero
 % there are treatment so no failed treatment
 % so only D10 and D00
 % only for the youngest age group since no tranbsistsion

%   if t>67
%       Mdash_acute(6,1) = 0;
%   end
 
 d00=Macute*X00-extra_parms_vals(2)*X00+extra_parms_vals(3)*X10-mort_former.*X00+(age_matrix*X00')';
 d01=Macute*X01-extra_parms_vals(2)*X01+extra_parms_vals(3)*X11-mort_former.*X01+(age_matrix*X01')';
 d10=Mdash_acute*X10+extra_parms_vals(2)*X00-extra_parms_vals(3)*X10-mort_current.*X10+(age_matrix*X10')';
 d11=Mdash_acute*X11+extra_parms_vals(2)*X01-extra_parms_vals(3)*X11-mort_current.*X11+(age_matrix*X11')';

%  d00=M*X00-extra_parms_vals(2)*X00+extra_parms_vals(3)*X10-mort_former.*X00+(age_matrix*X00')';
%  d01=M*X01-extra_parms_vals(2)*X01+extra_parms_vals(3)*X11-mort_former.*X01+(age_matrix*X01')';
%  d10=Mdash*X10+extra_parms_vals(2)*X00-extra_parms_vals(3)*X10-mort_current.*X10+(age_matrix*X10')';
%  d11=Mdash*X11+extra_parms_vals(2)*X01-extra_parms_vals(3)*X11-mort_current.*X11+(age_matrix*X11')';
if t>time2015_inODEunits 
    death_mult=1.01+death_val*(t-time2015_inODEunits)/15;
else
    death_mult=1;
end

death_mult=1;
d10(1,1) =d10(1,1) + death_mult*(sum(mort_former.*X00,'all') + sum(mort_former.*X01,'all') + sum(mort_current.*X10,'all') + sum(mort_current.*X11,'all')...
                   + death_vec(1)*sum(X00(12,:),'all')  + death_vec(2)*sum(X00(13,:),'all')  + death_vec(3)*sum(X00(14,:),'all')  + death_vec(4)*sum(X00(15,:),'all')...
                   + death_vec(1)*sum(X01(12,:),'all')  + death_vec(2)*sum(X01(13,:),'all')  + death_vec(3)*sum(X01(14,:),'all')  + death_vec(4)*sum(X01(15,:),'all')...
                   + death_vec(1)*sum(X10(12,:),'all')  + death_vec(2)*sum(X10(13,:),'all')  + death_vec(3)*sum(X10(14,:),'all')  + death_vec(4)*sum(X10(15,:),'all')...
                   + death_vec(1)*sum(X11(12,:),'all')  + death_vec(2)*sum(X11(13,:),'all')  + death_vec(3)*sum(X11(14,:),'all')  + death_vec(4)*sum(X11(15,:),'all'));
incibit=0;
[phi,phidash,~,~,~,~]=treat_comps_pv3_scaled_acute_v2(treatprob,y1,N,age_scale00,age_scale10,t,pcom,alpha,time2015_inODEunits,pwid_weights,former_weights,incibit);   

if t>68
    df=45
end

%  if sum(isnan(phi))>0
%       df = 56
%  end
%  if sum(isnan(phidash))>0
%      df = 56
%  end
%outyt=[outyt t];
%outyF4=[outyF4 moved_to_T4_former moved_to_T4_current,total_not_T4,total_T4];
%save(['/Users/celeste/HCV/version5/F4dat/f',replace(num2str(round(t,4)),'.','_')],'moved_to_T4')
% d00 is -phi00 + phidash00

 if t>67
     alpha=0.95;
     pcom=1;
     dum = percent_treated * Mdash_acute(6,1)*X10(1,:);
     %dum1 = percent_treated * Mdash_acute(6,1)*X11(1,:);
     phi(384:21:567)=dum;
     %phi1(384:21:567)=dum1;
     A_pwid_treat = 399:21:567;
     phidash(A_pwid_treat)=alpha*pcom*dum;

 end


trtmodel00=-phi(1:189)+phidash(1:189);
trtmodel01=(1-alpha)*phi(1:189); % treatment failures,alpha is in phidash above
trtmodel10=-phi(379:567)+phidash(379:567);
trtmodel11=(1-alpha*pcom)*phi(379:567); %treatment failures aplha*pcom is in phidash above

hcv=[reshape(d00,189,1)+trtmodel00;reshape(d01,189,1)+trtmodel01;...
     reshape(d10,189,1)+trtmodel10;reshape(d11,189,1)+trtmodel11];
 



 
 
 
 
                          
 
 
 
 
 
 
 
 
 
 
 