function [InfoTheo]=InfoTheoFun_mutual(X1,Xtar,nbin)
%=================================================================================
%Note that a trick of feeding two identical timeseries will give no-good answer
%This is due to the trncation error after fitting the linear model on the data
%The solution is to add a very-tiny nosie in the timeseries
%===================================================================================================

assert(numel(X1) == numel(Xtar));
n=numel(X1);

X1_=round((X1-min(X1))/(max(X1)-min(X1))*nbin+1);
Xtar_=round((Xtar-min(Xtar))/(max(Xtar)-min(Xtar))*nbin+1);

Xv=[1:nbin+1];
Yv=[1:nbin+1];
values={Xv Yv};
count=hist3(cat(2,X1_,Xtar_),'ctrs',values); %first count the values
p=count/sum(sum(count));
[p1,n1]=hist(X1_,Xv);
p1=p1/n;
[p2,n2]=hist(Xtar_,Yv);
p2=p2/n;
%Atom at 0 effect

for i=1:nbin+1
    for j=1:nbin+1
           % H(i,j,k)=log2(J3(i,j,k)/(p(i,j)*p3(k)));
PP(i,j)=p(i,j)*log2(p(i,j)/(p1(i)*p2(j)));
    end
end

MIV=nansum(PP(:));
MIxz=MIV;

H1_2D=-(p1.*log2(p1)); H1_2D=nansum(H1_2D);
H3_2D=-(p2.*log2(p2)); H3_2D=nansum(H3_2D);
Cri=Sigtest_mi(X1_,Xtar_,nbin,30);

%===================================================================================================
%                               2D part for I(X2,Y)
%                               2 Output: I(X2,Y), U2
%===================================================================================================
x=X1';
y=Xtar';
X = [ones(size(x')) x'];
b = regress(y',X); yy=b(2)*x+b(1);
z=y-yy; newz=cat(1,x,z);
z=sortrows(newz',2); J=sort(y);
z=z'; z(2,:)=J;
MIxz_non=MI(z(1,:)',z(2,:)',nbin);

%{
X = [ones(size(X1)) X1];
b = regress(Xtar,X);
yy=b(2)*X1+b(1);
z=Xtar-yy; newz=cat(2,X1,z);
z=sortrows(newz,2); J=sort(Xtar);
z(:,2)=J;
MIxz_non=MI(z(:,1),,nbin);
%}
MIxz_non=MIxz_non(1);
Cri_non=Sigtest_mi(z(:,1),z(:,2),nbin,30);

MIxz_lin=MIxz-MIxz_non;
InfoTheo=cat(1,MIxz,MIxz_non,MIxz_lin,H3_2D,H1_2D,Cri,Cri_non);





 


