function [inci,chron]=incidence(N,t,extra_parms_vals,scale_former,scale_current,param_current)
 % N are the compartments - there are
 % 80 (0,0 0,1 1,0 1,1) by 9 ages 
 % 0,0 are rows 1 : 20
 % 0,1 are rows 21 : 40
 % 1,0 are rows 41 : 60
 % 1,1 are rows 61 : 80
 % i,j notation: i = 0 is former PWID and i = 1 is current PWID
 %             : j = 0 is never failed trt and j = 1 is failed trt
 
 
 
 Irows = [6,7,8,9,10,11,12,13,14,15]; % rows of infectious PWID

 delta=param_current(1);
 %pa_to_mild = 1 - exp(-param_current(2));
 pa_to_mild = param_current(2);
 pmult=(1-delta)*pa_to_mild;
 
 % no now changed
%  top = sum(N(Irows+40,:),'all') + sum(N(Irows+60,:),'all'); % 40 is j = 0, 60 is j = 1 and i = 1 for both
%  bot = sum(N(((1:20)+40),:),'all') + sum(N((1:20)+60),'all'); % 40 is j = 0, 60 is j = 1 and i = 1 for both
 top10_index = reshape((repmat(Irows,9,1)+(360:20:520)')',length(Irows)*9,1);
 top11_index = reshape((repmat(Irows,9,1)+(540:20:700)')',length(Irows)*9,1);
 top = sum(N(top10_index))+sum(N(top11_index)); % 40 is j = 0, 60 is j = 1 and i = 1 for both
 bot = sum(N(361:end)); % 40 is j = 0, 60 is j = 1 and i = 1 for both
 I = top/bot;
 if t>=0 && t<50
     scaleI = 1+1.5*t/50;
 elseif t>=50 && t<=55
      scaleI = 2.5-1.5*(t-50)/5;
 else
     scaleI=1;
 end
 phi = scaleI*extra_parms_vals(1)*I;
 %phiprob = 1 -exp(-phi);
 inci=phi*(scale_former*N(1:20:180)+scale_former*N(181:20:360)+scale_current*N(361:20:540)+scale_current*N(541:20:720));
 inci=sum(inci);
 inci=(scale_former*N(6:20:180)+scale_former*N(186:20:360)+scale_current*N(366:20:540)+scale_current*N(546:20:720));
 inci=sum(inci);
 chron=pmult*(scale_former*N(6:20:180)+scale_former*N(186:20:360)+scale_current*N(366:20:540)+scale_current*N(546:20:720));
 chron=sum(chron);
