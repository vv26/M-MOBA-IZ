% Modified edition, e.g., according to the table and the idzone area is modified according to the discussion in Mar,2017 
clear;%% use same random number and switch to equal allocation when pj=0
n0=5;
%Mu=[0.5 5.5;1.9 4.2;2.8 3.3;3 3;3.9 2.1;4.3 1.8;4.6 1.5;3.8 6.3;4.8 5.5;5.2 5;5.9 4.1;6.3 3.8;6.7 7.2;7 7;7.9 6.1;9 9];
%Mu=[1 2;3 1;5 5];
%Mu=[1 8;2 5;2.3 5.3;4 2;4.1 1.2;3 7;3.8 3;1.5 6;3 4;3.5 8;];
Mu=[1 8;2 5;3.5 5.01;3 2;2.5 8;3 7;3.05 2.2;1.5 6;2.1 5.2;2.5 4;2.6 3.9;2 7;2.5 6];
[r,cl]=size(Mu);tsig=1.5*ones(r,2);
% tsig=1+(3-1).*rand(r,2);%% random number[1,3]
matlabzero=10^(-6);
f0=paretofront(r,Mu);

d1=0.2;
d2=0.2;
A=sort(Mu(f0==1,1));
B=sort(Mu(f0==1,2),'descend');
A=[A;inf];B=[inf;B];
%p=sum(f0)+1;

% for i=1:r
%     if f0(i)==0
% %         for j=1:p
% %             if (Mu(i,1)>=A(j)&&Mu(i,1)<=(A(j)+d1)&&Mu(i,2)>=B(j+1)&&Mu(i,2)<=(B(j)+d2))
% %                 f00(i)=(1==1);%%%%to make sure f00(i) is logic 1 instead of simple integer 1
% %                 break;
% %             end
% %             if (Mu(i,1)>=(A(j)+d1)&&Mu(i,1)<=A(j+1)&&Mu(i,2)>=B(j+1)&&Mu(i,2)<=(B(j+1)+d2))
% %                 f00(i)=(1==1);
% %                 break;
% %             end
% %         end
%     end
% end
 f00=borderline(r,f0,Mu,d1,d2);



jn=2;   %%how many different values of budget will be tested
ct=zeros(3,jn);%%the times of the right choice;ct(1,) equal allocation; ct(2,) allocation accoding to probability
cti=zeros(3,jn);
T=zeros(4,jn);
budgets=0;%%the smallest budget
budgeti=1;%%budget step
mre=1; %% the maximum repitation
record1=zeros(r,mre,jn);%% record in each iteration which alternative gets the sampling and the pj accordingly
record2=zeros(r,mre,jn);%% record in each iteration which alternative gets the sampling and the pj accordingly
recordmobaequal=zeros(mre,jn);
recordidequal=zeros(mre,jn);
for re=1:mre
    for i=1:r%++++++++++++++++++++++++++++++++++++++
        for j=1:cl
            sps(i,j,1:jn)=normrnd(Mu(i,j),tsig(i,j),jn,1);
        end;
    end;

%%%%%%%%%%%%%%%the initial sample mean value and variance%%%%%%%%%%%%%%%%%%
    xb0=zeros(r,2);sig0=zeros(r,2);
    for i=1:n0
        X(1:r,2*i-1:2*i)= (normrnd(Mu,tsig));
        xb0= (xb0+X(1:r,2*i-1:2*i));
    end
    xb0 = (xb0/n0);
    for i=1:n0
        sig0 = (sig0+(X(1:r,2*i-1:2*i)-xb0).^2);
    end
    sig0 = (sqrt(sig0/(n0-1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%euqal allocation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xb=xb0;sig=sig0;n=n0*ones(r,1);budget=0;
    
    for j=1:jn
        tbudget=budgets+(j-1)*budgeti;
        while(budget<tbudget)
            mn=mod(budget,r)+1;
            budget=budget+1;
            X=sps(mn,:,n(mn)-n0+1);%++++++++++++++++++
            xb(mn,1:2)= ((n(mn)*xb(mn,1:2)+X)/(n(mn)+1));
            n(mn)=n(mn)+1;
        end;
        f1=paretofront(r,xb);
        
        f10=borderline(r,f1,xb,d1,d2);
        ct(1,j)=ct(1,j)+(sum(f1==f0)==r);
        cti(1,j)=cti(1,j)+((sum(f00(f1~=f0)&f10(f1~=f0))+sum(f1==f0))==r); %???sum?????????borderline????sum????
        
        %%cti(1,j)=cti(1,j)+(sum(f1(~f00)==f0(~f00))==(r-sum(f00)));
    end;

%%%%%%%%%%%%%%%%allocation by indifference zone%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xb=xb0; sig=sig0;n=n0*ones(r,1);budget=0;nonalgorithm=0;
    
    for k=1:r
        a=xb(k,1)<xb([1:k-1,k+1:r],1);
        b=xb(k,2)<xb([1:k-1,k+1:r],2);
        c=a|b;
        cn=sum(c);
        f1(k)=(cn==(r-1));%%if point k is on the pareto front based on five observations, f1(k)=1,0 otherwise
    end;
    
    for j=1:jn
        tbudget=budgets+(j-1)*budgeti;
        while(budget<tbudget)
            pj=(paretotidz(xb,sig,f1,f10,r,n,1,d1,d2));
            if sum(pj)>matlabzero
                [m,mn]=max(pj);%%m is the largest number of pj and mn is the correspoding alternative's sequence number 
                %[m,mn]= min(pj); 
                T(1,j)=T(1,j)+1;%%how many times use algorithm
                recordidequal(re,j)=mn+100;
            else
                pj=(paretotidz(xb,sig,f1,f10,r,n,10,d1,d2));%% switch twice to tau=10
                if sum(pj)>matlabzero
                    [m,mn]=max(pj);
                    T(1,j)=T(1,j)+1;
                    recordidequal(re,j)=mn+200;
                else
                    T(2,j)=T(2,j)+1; 
                    mn=mod(nonalgorithm,r)+1;
                    recordidequal(re,j)=mn+300;
                    nonalgorithm=nonalgorithm+1;
                end
            end
            
            X=sps(mn,:,n(mn)-n0+1);%++++++++++++++++++++++++++++++++++++++
            sig(mn,1:2)=(sqrt((n(mn)-1)/n(mn)*sig(mn,1:2).^2+1/(n(mn)+1)*(X-xb(mn,1:2)).^2));
            xb(mn,1:2)=((n(mn)*xb(mn,1:2)+X)/(n(mn)+1));
            n(mn)=n(mn)+1;
            budget=budget+1;
            
            for k=1:r 
                a=xb(k,1)<xb([1:k-1,k+1:r],1);
                b=xb(k,2)<xb([1:k-1,k+1:r],2);
                c=a|b;
                cn=sum(c);
                f1(k)=(cn==(r-1));
            end;
            record1(:,re,j)=pj;%%record in each iteration which alternative gets the sampling and the pj accordingly
        end;
        
        f10=borderline(r,f1,xb,d1,d2);%% dont need be upper since while 
        ct(2,j)=ct(2,j)+(sum(f1==f0)==r);
        cti(2,j)=cti(2,j)+((sum(f00(f1~=f0)&f10(f1~=f0))+sum(f1==f0))==r);
    end;
     %%%%%%%%%%%%%%%%allocation by mmoba%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xb=xb0; sig=sig0;n=n0*ones(r,1);budget=0;nonalgorithm=0;
    
    for k=1:r
        a=xb(k,1)<xb([1:k-1,k+1:r],1);
        b=xb(k,2)<xb([1:k-1,k+1:r],2);
        c=a|b;
        cn=sum(c);
        f1(k)=(cn==(r-1));%%if point k is on the pareto front based on five observations, f1(k)=1,0 otherwise
    end
    for j=1:jn
        tbudget=budgets+(j-1)*budgeti;
        while(budget<tbudget)
            pj= (paretot(xb,sig,f1,r,n,1));
            if sum(pj)>matlabzero
                [m,mn]=max(pj);%%m is the largest number of pj and mn is the correspoding alternative's sequence number 
                %[m,mn]= min(pj); 
                T(3,j)=T(3,j)+1;%%how many times use algorithm
                recordmobaequal(re,j)=mn+100;%%record in each iteration which alternative gets the sampling and the pj accordingly
            else
                pj=(paretot(xb,sig,f1,r,n,10));%% switch twice to tau=10
                if sum(pj)>matlabzero
                    [m,mn]=max(pj);
                    T(3,j)=T(3,j)+1;%%how many times use algorithm
                    recordmobaequal(re,j)=mn+200;%%record in each iteration which alternative gets the sampling and the pj accordingly
                else
                    T(4,j)=T(4,j)+1; 
                    mn=mod(nonalgorithm,r)+1;
                    nonalgorithm=nonalgorithm+1;
                    recordmobaequal(re,j)=mn+300;
                end
            end
            
            X=sps(mn,:,n(mn)-n0+1);%++++++++++++++++++++++++++++++++++++++
            sig(mn,1:2)=(sqrt((n(mn)-1)/n(mn)*sig(mn,1:2).^2+1/(n(mn)+1)*(X-xb(mn,1:2)).^2));
            xb(mn,1:2)=((n(mn)*xb(mn,1:2)+X)/(n(mn)+1));
            n(mn)=n(mn)+1;
            budget=budget+1;
            
            for k=1:r 
                a=xb(k,1)<xb([1:k-1,k+1:r],1);
                b=xb(k,2)<xb([1:k-1,k+1:r],2);
                c=a|b;
                cn=sum(c);
                f1(k)=(cn==(r-1));
            end
            record2(:,re,j)=pj;%%record in each iteration which alternative gets the sampling and the pj accordingly
        end
        f10=borderline(r,f1,xb,d1,d2);
        ct(3,j)=ct(3,j)+(sum(f1==f0)==r);
        cti(3,j)=cti(3,j)+((sum(f00(f1~=f0)&f10(f1~=f0))+sum(f1==f0))==r);
    end
    fprintf('re=%d\n',re);
end
x=budgets+[0:(jn-1)]*budgeti;
ct1=ct(1,1:jn);
ct2=ct(2,1:jn);
ct3=ct(3,1:jn);
cp1=ct1/mre;
cp2=ct2/mre;
cp3=ct3/mre;
figure
plot(x,cp1,'-r',x,cp2,'-.b',x,cp3,'-.g')
legend('Equal','MMOBA-ID','MMOBA PCS')
xlabel('budget')
ylabel('P{CS}')
title('Probability of Correct Selection');

ct1=cti(1,1:jn);
ct2=cti(2,1:jn);
ct3=cti(3,1:jn);
cp1=ct1/mre;
cp2=ct2/mre;
cp3=ct3/mre;
figure
plot(x,cp1,'-r',x,cp2,'-.b',x,cp3,'-.g')
legend('Equal','MMOBA-ID','MMOBA PCS')
xlabel('budget')
ylabel('P{GS}')
title('Probability of Good Selection with Indifference Zone');