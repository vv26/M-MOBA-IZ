function f0=paretofront(r,Mu)
for k=1:r 
    a=Mu(k,1)<Mu([1:k-1,k+1:r],1);
    b=Mu(k,2)<Mu([1:k-1,k+1:r],2);
    c=a|b;
    cn=sum(c);
    f0(k)=(cn==(r-1));%% if point k is on the pareto front, f0(k)=1,0 otherwise
end;