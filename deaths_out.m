function [D00_non_liver,D01_non_liver,D10_non_liver,D11_non_liver,D00_liver,D01_liver,D10_liver,D11_liver]=deaths_out(Nin,tbeg,tint,tend,mort_former,mort_current,death_vec)

% mort_former are the non liver related deaths
% mort_current as above
% death_vec are the DC,HCC,LT1,LT2 deaths

% N is the 720 by tstart:tinit:tend matrix
% cols represent the times
t=tbeg:tint:tend;
icol=0;
D00_non_liver = zeros(180,tend-tbeg);
for i = tbeg:(tend-1)
    j=find(t==i):find(t==(i+1));
    icol=icol+1;
    for k = 1 : length(j)
        N=Nin(:,j(k));
        X00{k} = reshape(N(1:20*9),20,9);
        X01{k} = reshape(N(181:(180+20*9)),20,9);
        X10{k} = reshape(N(361:(360+20*9)),20,9);
        X11{k} = reshape(N(541:(540+20*9)),20,9);
    end
     D00_non_liver(:,icol)=integrate_deaths(mort_former,X00,t(j));
     D01_non_liver(:,icol)=integrate_deaths(mort_former,X01,t(j));
     D10_non_liver(:,icol)=integrate_deaths(mort_current,X10,t(j));
     D11_non_liver(:,icol)=integrate_deaths(mort_current,X11,t(j));
     
     D00_liver(:,icol)=liver_integrate_deaths(death_vec,X00,t(j));
     D01_liver(:,icol)=liver_integrate_deaths(death_vec,X01,t(j));
     D10_liver(:,icol)=liver_integrate_deaths(death_vec,X10,t(j));
     D11_liver(:,icol)=liver_integrate_deaths(death_vec,X11,t(j));
end

%D00_non_liver{i}=sum(mort_former.*X00,'all');
%D01_non_liver{i}=sum(mort_former.*X01,'all');
%D10_non_liver{i}=sum(mort_current.*X10,'all');
%D11_non_liver{i}=sum(mort_current.*X11,'all');
%D00_liver{i}=death_vec(1)*sum(X00(12,:),'all')  + death_vec(2)*sum(X00(13,:),'all')  + death_vec(3)*sum(X00(14,:),'all')  + death_vec(4)*sum(X00(15,:),'all');
%D01_liver{i}=death_vec(1)*sum(X01(12,:),'all')  + death_vec(2)*sum(X01(13,:),'all')  + death_vec(3)*sum(X01(14,:),'all')  + death_vec(4)*sum(X01(15,:),'all');
%D10_liver{i}=death_vec(1)*sum(X10(12,:),'all')  + death_vec(2)*sum(X10(13,:),'all')  + death_vec(3)*sum(X10(14,:),'all')  + death_vec(4)*sum(X10(15,:),'all');
%D11_liver{i}=death_vec(1)*sum(X11(12,:),'all')  + death_vec(2)*sum(X11(13,:),'all')  + death_vec(3)*sum(X11(14,:),'all')  + death_vec(4)*sum(X11(15,:),'all');
end

function Ivec=integrate_deaths(mult,X,t)
% integrates over the deaths in 1 year
% writes as a column ength 180, 20 lots of age group 1, 2 etc
% t is the current interval i.e 1 to 2 in steps of 1/12
Iout=zeros(20,9);
n=length(t);
for i=12:15
    for j = 1 : 9
        y=zeros(n,1);
        for k=1:n
            dum=X{k}; % the 20 by 9 matrix for time k
            y(k) = dum(i,j);
        end
        Iout(i,j)=trapz(t,y);
    end
end
Ivec=reshape(mult.*Iout,20*9,1);
end

function Ivec=liver_integrate_deaths(death_vec,X,t)
% integrates over the deaths in 1 year
% writes as a column ength 180, 20 lots of age group 1, 2 etc
% t is the current interval i.e 1 to 2 in steps of 1/12
mult=[death_vec(1)*ones(1,9);death_vec(2)*ones(1,9);death_vec(3)*ones(1,9);death_vec(4)*ones(1,9)];
Iout=zeros(4,9);
n=length(t);
for i=12:15
    for j = 1 : 9
        y=zeros(n,1);
        for k=1:n
            dum=X{k}; % the 20 by 9 matrix for time k
            y(k) = dum(i,j);
        end
        Iout((i-11),j)=trapz(t,y);
    end
end
Ivec=reshape(mult.*Iout,4*9,1);
end


