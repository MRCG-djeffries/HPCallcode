function Qout=qalydeaths(D,discount_rate,ystart,yend)
% Fatal loss qualy
% from
% Supporting Information
% Optimising age coverage of seasonal influenza
% vaccination in England: A mathematical and health
% economic evaluation
% Edward M. Hill1*, Stavros Petrou2,3, Henry Forster4, Simon de Lusignan3,5, Ivelina Yonova3,5, Matt
% J. Keeling1

% D is deaths per age by years
% age are the 9 age bands
% years are the years of model output, example 2016 to 2030
if discount_rate >1
    discount_rate=discount_rate/100;
end
group_ages = [22.5 27.5 32.5 40 50 60 70 80 90];
% calculate E(a) the age specific quality of life weight at age a
E=zeros(length(group_ages),1);
for a = 1 : length(group_ages)
    E(a) = 0;
    for i = 1 : life_expectancy(group_ages(a))
        E(a) = E(a) + quolw(group_ages(a)+i)/(1+discount_rate)^i;
    end
end
Q=zeros(length(group_ages),1);
for a = 1 : length(group_ages)
    Q(a)=0;
    for y=ystart:yend
        i=y-ystart+1;
        Q(a) = Q(a) + D(a,i)*E(a)*(1/(1+discount_rate))^(y-ystart);
    end
end
Qout=sum(Q);


function Qw=quolw(age)

Qw = interp1([22.5 27.5 32.5 40 50 60 70 80 90],...
             [0.929 0.919 0.919 0.893 0.855 0.810 0.773 0.703 0.703],...
             age,'linear',-10);
if sum(Qw==-10)>0
    j=find(Qw==-10);
    for i = 1 : length(j)
        if age(j(i))<22.5
            Qw(j(i)) = 0.929;
        else
            Qw(j(i)) = 0.703;
        end
    end
end