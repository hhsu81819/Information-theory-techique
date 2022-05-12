function allMMI = MMIfix(X1,X2,Xtar,nbin)

%Writer: Hsin Hsu
%This is the main function for calculating MMI and its full decomposition
% Xi:timeseries of variable
% nbin: its cube cuts the 3d pdf into the corresponding number of pixels
%=== list of output
%   1   2    3    4    5    6    7    8  9  10 11 12
%1  MMI H(x) H(y) H(z) MIxz MIyz MIxy Ux Uy R  S  Critirion % For total_MMI
%2  MMI H(x) H(y) H(z) MIxz MIyz MIxy Ux Uy R  S  Critirion % For nonlinear_MMI 
%3  MMI H(x) H(y) H(z) MIxz MIyz MIxy Ux Uy R  S  Critirion % For linear_MMI
%===============================================================
%Note that a trick of feeding three identical timeseries will give no-good answer
%This is due to the trncation error after fitting the linear model on the data
%The solution is to add a very-tiny nosie in the timeseries

totMMI=MMI_fixedbin_partitioning(X1,X2,Xtar,nbin);

%===remove linear information====
X = [ones(size(X1)) X1 X2];
b = regress(Xtar,X);
yy=b(3)*X2+b(2)*X1+b(1);
z=Xtar-yy; newz=cat(2,X1,X2,z);
z=sortrows(newz,3); y=sort(Xtar);
z(:,3)=y;
%================================

nonMMI=MMI_fixedbin_partitioning(z(:,1),z(:,2),z(:,3),nbin);
linMMI=totMMI-nonMMI; linMMI(12)=nan;
allMMI=cat(1,totMMI',nonMMI',linMMI');


