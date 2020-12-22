function plotscale

T=0:82;
for i = 1 : length(T)
    t=T(i);
if t>=0 && t<=51
     scaleI(i) = 1+1.5*t/51;
 elseif t>51 && t<=56
      scaleI(i) = 2.5-1.5*(t-51)/5;
%  elseif t>56 && t<=66
%      scaleI=1;
  else
     scaleI(i)=1;
end
end
figure;plot(T,scaleI,'linewidth',2)
xticks(0:5:80)
set(gca,'xticklabel',1950:5:2031,'xticklabelrotation',90)
xlabel('Year')
ylabel('Scale')
set(gca,'ylim',[0 2.5])
xlim([0 81])
grid()
