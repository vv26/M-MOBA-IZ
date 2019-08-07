function pj=paretotidz(xb,sig,f1,f10,r,n,tau,d1,d2)

for k=1:r
    pj(k)=0;
    if f1(k)==1
        f2=zeros(1,r);
        for i=[1:k-1,k+1:r]
            a=xb(i,1)<=xb([1:k-1,k+1:r],1);
            b=xb(i,2)<=xb([1:k-1,k+1:r],2);
            c=a|b;
            cn=sum(c);
            f2(i)=(cn==(r-1));
        end;
    else
        f2=f1;
    end;
    
    tmp=xb(:,1);
    tmp(f2==0)=-inf;
    [A,num]=sort(tmp);
    num=num(A>-inf);
    A=A(A>-inf);
    
    A=sort(xb(f2==1,1));
    B=sort(xb(f2==1,2));
    %     Af1=sort(xb(f1==1,1));%% PF alternatives on f1
    %     Bf1=sort(xb(f1==1,2));
    l=sum(f2);
    
    A1=[-inf;A;inf];
    B1=[-inf;B;inf];
    %     A1f1=[Af1;inf];
    %     B1f1=[Bf1;inf];
    
    if f1(k)==1
        if l==sum(f1)-1 %% if no new solutions on PF when solution k is removed
            % do the following no matter if f10(k)==0 or f10(k)==1
            % on the left of f2 pareto
            for i=1:l+1
                if i<=l
                    j=num(i);
                end
                if i>l || f10(j)==0   %% for a PF is ND
                    xab=[A1(i),A1(i+1)];yab=[B1(l-i+2),B1(l-i+3)];
                    pj(k)=pj(k)+td(xab,yab,xb,sig,k,n,tau);
                    if td(xab,yab,xb,sig,k,n,tau)<0
                        error('negative value 1')
                    end
                else
                    xab=[A1(i),A1(i+1)];yab=[B1(l-i+2)-d2,B1(l-i+3)];
                    pj(k)=pj(k)+td(xab,yab,xb,sig,k,n,tau);
                    if td(xab,yab,xb,sig,k,n,tau)<0
                        error('negative value 2')
                    end
                    %                     if i<=l
                    xab=[A1(i+1)-d1,A1(i+1)];yab=[B1(l-i+1),B1(l-i+2)-d2];
                    pj(k)=pj(k)+td(xab,yab,xb,sig,k,n,tau);%% need to add the stripe area for each possible alternative
                    %                         if td(xab,yab,xb,sig,k,n,tau)<0
                    %                             error('negative value 3, k is %d',k)
                    %                         end
                    if  (i>1) && (f10(num(i-1))==1)
                        xab = [A1(i)-d1,A1(i)]; yab = [B1(l-i+2)-d2, B1(l-i+2)];
                        pj(k)=pj(k)+td(xab,yab,xb,sig,k,n,tau);
                        if td(xab,yab,xb,sig,k,n,tau)<0
                            error('negative value 4')
                        end
                    end
                    %                     end
                end
            end;
            
            if f10(k)==1 %% k is BND, more area on the right of f2
                A2 = [A;inf];
                for i=1:l
                    xab1= [A2(i),A2(i+1)];yab1=[B(l+1-i),inf];
                    xab2= xab1+d1; yab2=yab1+d2;
                    pj(k)=pj(k)+td(xab1,yab1,xb,sig,k,n,tau)-td(xab2,yab2,xb,sig,k,n,tau);
                end
            end;
        else %% new solutions on PF when k is removed
            newpfpb=0;%%%%%becomes 1 when new pareto front point was in borderline area
            for i=1:l
                j=num(i);
                if (f1(j)==0) && (f10(j)==1)
                    newpfpb=1;
                    break
                end
            end
            if newpfpb ==0
                [xab(1),i]=max(A1(A1<xb(k,1))); xab(2)=min(A1(A1>xb(k,1)));
                if i>1
                    i=i-1;
                    j=num(i);
                    if f10(j) ==1
                        xab(1) = xab(1)-d1;
                    end
                end
                
                [yab(1),i]=max(B1(B1<xb(k,2)));yab(2)=min(B1(B1>xb(k,2)));
                if i>1
                    i = l+2-i;
                    j=num(i);
                    if f10(j) ==1
                        yab(1) = yab(1)-d2;
                    end
                end
                pj(k)=pj(k)+td(xab,yab,xb,sig,k,n,tau);
                %                 if td(xab,yab,xb,sig,k,n,tau)<0
                %                     error('negative value 7')
                %                 end
            else
                
                Ad1 = [-inf;sort(unique([A; A-d1; A+d1]));inf];
                Bd2 = [-inf;sort(unique([B; B-d2; B+d2]));inf];
                A12 = [-inf; sort(xb(f1&f2,1)); inf];
                B12 = [-inf; sort(xb(f1&f2,2)); inf];
                ll=sum(f1&f2);
                for i=1:ll+1
                    a1 = [A12(i)-d1,A12(i+1)+d1];
                    b1 = [B12(ll+2-i)+d2,B12(ll+3-i)-d2];
                    x1 = Ad1( (Ad1<=a1(2)) & (Ad1>=a1(1)) );
                    if b1(2)>b1(1)
                        y1 = Bd2( (Bd2<=b1(2)) & (Bd2>=b1(1)) );
                        pm=1;%%%%pm is short for plus or minus
                    else
                        y1 = Bd2( (Bd2<=b1(1)) & (Bd2>=b1(2)) );
                        pm=-1;%%%%pm is short for plus or minus
                    end
                    [m1, m2] =size(x1);
                    [n1, n2] =size(y1);
                    for p = 1:m1-1
                        for q = 1:n1-1
                            xab = [x1(p),x1(p+1)];
                            yab = [y1(q),y1(q+1)];
                            
                            xbn = xb;%%%%%%% xbn stores the coordinates of all points with kth alternative moving to a possible new location
                            xbn(k,1) = sum(xab)/2;
                            xbn(k,2) = sum(yab)/2;
                            f3=paretofront(r,xbn);
                            f30=borderline(r,f3,xbn,d1,d2);
                            if f3==f1
                                pj(k)=pj(k)+pm*td(xab,yab,xb,sig,k,n,tau);
                            elseif (f30(f3~=f1)&f10(f3~=f1))==1
                                pj(k)=pj(k)+pm*td(xab,yab,xb,sig,k,n,tau);
                            end
                        end
                    end
                end
                for i=1:ll
                    a1 = [A12(i)-d1,A12(i+2)+d1];
                    b1 = [B12(ll+2-i)-d2,B12(ll+2-i)+d2];
                    if a1(1)<a1(2)
                        x1 = Ad1( (Ad1<=a1(2)) & (Ad1>=a1(1)) );
                        pm=1;
                    else
                        x1 = Ad1( (Ad1<=a1(1)) & (Ad1>=a1(2)) );
                        pm=-1;
                    end
                    y1 = Bd2( (Bd2<=b1(2)) & (Bd2>=b1(1)) );
                    [m1, m2] =size(x1);
                    [n1, n2] =size(y1);
                    for p = 1:m1-1
                        for q = 1:n1-1
                            xab = [x1(p),x1(p+1)];
                            yab = [y1(q),y1(q+1)];
                            
                            xbn = xb;%%%%%%% xbn stores the coordinates of all points with kth alternative moving to a possible new location
                            xbn(k,1) = sum(xab)/2;
                            xbn(k,2) = sum(yab)/2;
                            f3=paretofront(r,xbn);
                            f30=borderline(r,f3,xbn,d1,d2);
                            if f3==f1
                                pj(k)=pj(k)+pm*td(xab,yab,xb,sig,k,n,tau);
                            elseif (f30(f3~=f1)&f10(f3~=f1))==1
                                pj(k)=pj(k)+pm*td(xab,yab,xb,sig,k,n,tau);
                            end
                        end
                    end
                end
            end
        end;
    else
        A=[A;inf];
        if f10(k)==0
            for i=1:l
                xab=[A(i),A(i+1)];yab=[B(l-i+1),inf]; %% if k is not borderlined, cant move out of the original area
                pj(k)=pj(k)+td(xab,yab,xb,sig,k,n,tau);
                %                 if td(xab,yab,xb,sig,k,n,tau)<0
                %                     error('negative value 10')
                %                 end
            end;
        else
            for i=1:l
                j=num(i);
                xab=[A(i),A(i+1)]-d1;yab=[B(l-i+1),inf]-d2;
                pj(k)=pj(k)+td(xab,yab,xb,sig,k,n,tau); %% stripe area doesn't need to minus small rectangles
                %                 if td(xab,yab,xb,sig,k,n,tau)<0
                %                     error('negative value 11')
                %                 end
                if f10(j)==0  %% if k is borderlined and jth alternative is not borderlined?i????????A??????j???????alternative?????
                    xab=[A(i)-d1,A(i)];yab=[B(l-i+1)-d2,B(l-i+1)];
                    pj(k)=pj(k)-td(xab,yab,xb,sig,k,n,tau); %% stripe area minus small rectangles
                    %                     if td(xab,yab,xb,sig,k,n,tau)<0
                    %                         error('negative value 12')
                    %                     end
                end;
            end;
        end;
    end;
%     if pj(k)<0
%         error('pj(%d) is %f, negative value 13, k is%f, f1 is %d, f10 is %d',k,pj(k),f1(k),f10(k))
%     elseif pj(k)>1
%         if l-sum(f1)+1>0
%             error('pj(%d) is %f, bigger than 1, f1 is %d, f10 is %d. l-sum(f1)+1 is %d, newpfpb is %d',k,pj(k),f1(k),f10(k),l-sum(f1)+1,newpfpb)
%         else
%             error('pj(%d) is %f, bigger than 1, f1 is %d, f10 is %d. l-sum(f1)+1 is %d',k,pj(k),f1(k),f10(k),l-sum(f1)+1)
%         end
%     end
end;
pj=1-pj;