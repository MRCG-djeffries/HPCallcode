function [total_liver_transplants]=costdetail(X,LT1_comps_00,LT1_comps_10,scale00,scale10,time2015,param_current)

% Get the number of liver transplants, by integrating over the LT
% compartment

compvec= [LT1_comps_00;LT1_comps_10];
liver_vals=zeros(2,1); % This conatins qualy for each chronic department
scalevals=[scale00,scale10];
for i = 1 : 2
    Xtest=scalevals(i)*sum(X(compvec(i,:),:));
    liver_vals(i)=calc_qalyv2(Xtest(time2015:961),2015,2030,1/12,1,0);
end
total_liver_transplants= sum(liver_vals);

% Calculate the same using the incidence
probdclt = param_current(13);
probhcclt = param_current(15);


function [inci,chron]=incidenceLT(N,t,rdclt,rhcclt,scale_former,scale_current)
 % N are the compartments - there are
 % 80 (0,0 0,1 1,0 1,1) by 9 ages 
 % 0,0 are rows 1 : 20
 % 0,1 are rows 21 : 40
 % 1,0 are rows 41 : 60
 % 1,1 are rows 61 : 80
 % i,j notation: i = 0 is former PWID and i = 1 is current PWID
 %             : j = 0 is never failed trt and j = 1 is failed trt
 
 % This calculates the incidence into the Liver transplant LT! compartment
 
 Irows = [6,7,8,9,10,11,12,13,14,15]; % rows of infectious PWID

 delta=param_current(1);
 pa_to_mild = 1 - exp(-param_current(2));
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
 phiprob = 1 -exp(-phi);
 inci=phiprob*(scale_former*N(1:20:180)+scale_former*N(181:20:360)+scale_current*N(361:20:540)+scale_current*N(541:20:720));
 inci=sum(inci);
 chron=pmult*inci;