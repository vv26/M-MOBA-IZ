function f00=borderline(r,f0,Mu,d1,d2)%% when i is in the borderline area, f00(i)=1
for i=1:r
    if f0(i)==0
        a=Mu(i,1)-d1< Mu([1:i-1,i+1:r],1);
        b=Mu(i,2)-d2< Mu([1:i-1,i+1:r],2);
        c=a|b;
        f00(i)=(sum(c)==r-1);
    else
        a=Mu(i,1)+d1< Mu([1:i-1,i+1:r],1);
        b=Mu(i,2)+d2< Mu([1:i-1,i+1:r],2);
        c=a|b;
        f00(i)=(sum(c)~=r-1);
    end
end