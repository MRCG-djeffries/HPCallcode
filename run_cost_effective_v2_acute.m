function run_cost_effective_v2_acute()
runy=0;
if runy==1
for i= 50:100
    p=i/100;
    cost_effective_v2_acute(p)
end
else
    load(['foutaplha1100'])
    inci1=incy;
    load(['foutaplha99100'])
    inci2=incy;    
    load(['foutaplha95100'])
    inci3=incy;
    figure;plot(2016:2030,inci1,'r-','LineWidth',2);
            grid()
            title('% incidence change in current PWID becoming chronic')
            set(gca,'xlim',[2016 2030])
            hold on;
plot(2016:2030,inci2,'g-','LineWidth',2);
plot(2016:2030,inci3,'b-','LineWidth',2);
legend({'100%','99%','95%'})
set(gca,'ylim',[0 100])
ylabel('%')            
    
%     for i= 50:100
%         load(['fout',num2str(i)])
%         if i==50
%             figure;plot(2016:2030,incy,'r-','LineWidth',2);
%             grid()
%             title('% incidence change in current PWID becoming chronic')
%             set(gca,'xlim',[2016 2030])
%         else
%             hold on
%             plot(2016:2030,incy,'r-','LineWidth',2);
%         end
%     end
end

