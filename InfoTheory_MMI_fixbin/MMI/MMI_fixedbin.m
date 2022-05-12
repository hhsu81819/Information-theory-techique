function MMI = MMIfix(X1,X2,Xtar,nbin)

%Writer: Hsin Hsu
assert(numel(X1) == numel(X2));
assert(numel(X1) == numel(Xtar));
n=numel(X1);

X1p=round((X1-min(X1))/(max(X1)-min(X1))*nbin+1);
X2p=round((X2-min(X2))/(max(X2)-min(X2))*nbin+1);
Xtarp=round((Xtar-min(Xtar))/(max(Xtar)-min(Xtar))*nbin+1);

J3=accumarray([X1p X2p Xtarp],1)/(n);


Xv=[1:nbin+1];
Yv=[1:nbin+1];
values={Xv Yv};
count=hist3(cat(2,X1p,X2p),'ctrs',values); %first count the values
pxy=count/sum(sum(count));

Xv=[1:nbin+1];
Yv=[1:nbin+1];
values={Xv Yv};
count=hist3(cat(2,X2p,Xtarp),'ctrs',values); %first count the values
pyz=count/sum(sum(count));

Xv=[1:nbin+1];
Yv=[1:nbin+1];
values={Xv Yv};
count=hist3(cat(2,X1p,Xtarp),'ctrs',values); %first count the values
pxz=count/sum(sum(count));

p3=squeeze(sum(sum(J3,1),2));
p2=squeeze(sum(sum(J3,1),3));
p1=squeeze(sum(sum(J3,2),3));

MMI=MultiMutualInformation(J3,pxy,p3);

