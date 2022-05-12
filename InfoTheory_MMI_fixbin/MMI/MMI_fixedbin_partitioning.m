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

MMIxyz=MultiMutualInformation(J3,pxy,p3);
H1=nansum(-(p1.*log2(p1)));
H2=nansum(-(p2.*log2(p2)));
H3=nansum(-(p3.*log2(p3)));

MIxz=MutualInformation(pxz,p1,p3);             % I(Xs1;Xtar)
MIyz=MutualInformation(pyz,p2,p3);             % I(Xs2,Xtar)
MIxy=MutualInformation(pxy,p1,p2);          % I(Xs1,Xs2)
CMI_x_zy=MMIxyz-MIyz;          % I(Xs1;Xtar|Xs2)
CMI_y_zx=MMIxyz-MIxz;          % I(Xs2;Xtar|Xs1)
II=CMI_x_zy-MIxz;       % II=I(Xs2;Xs1;Xtar)
RMMI=min(MIxz,MIyz);              % R  : minimum shared information
                Rmin=max(0,-1*II);
                Ls=MIxy/(min(H1,H2));
RescaleR=(Rmin+Ls*(RMMI-Rmin));          % Output: Rs/H(Y)  : Rescaled Redundancy
S=II+RescaleR;                      % Output: S/H(Y)
Uy=MIyz-RescaleR;
Ux=MIxz-RescaleR;

Cri=MMI_Sigtest_Surrogate(X1,X2,Xtar,nbin,30);
MMI=cat(1,MMIxyz,H1,H2,H3,MIxz,MIyz,MIxy,Ux,Uy,RescaleR,S,Cri);


