function MI = MI(X1,Xtar,nbin)

assert(numel(X1) == numel(Xtar));
n=numel(X1);

X1=round((X1-min(X1))/(max(X1)-min(X1))*nbin+1);
Xtar=round((Xtar-min(Xtar))/(max(Xtar)-min(Xtar))*nbin+1);

Xv=[1:nbin+1];
Yv=[1:nbin+1];
values={Xv Yv};
count=hist3(cat(2,X1,Xtar),'ctrs',values); %first count the values
p=count/sum(sum(count)); 
[p1,n1]=hist(X1,Xv);
p1=p1/n;
[p2,n2]=hist(Xtar,Yv);
p2=p2/n;


for i=1:nbin+1
    for j=1:nbin+1
           % H(i,j,k)=log2(J3(i,j,k)/(p(i,j)*p3(k)));
PP(i,j)=p(i,j)*log2(p(i,j)/(p1(i)*p2(j)));
    end
end

MIV=nansum(PP(:));
MI=MIV;
