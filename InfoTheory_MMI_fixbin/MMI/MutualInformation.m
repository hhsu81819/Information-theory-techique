function MI = MI(pxy,p1,p2)
bin=size(pxy,1);
p1=reshape(p1,[bin 1]);
p2=reshape(p2,[1 bin]);
LL=pxy./(p1*p2);

PP=pxy.*log2(LL);

MI=nansum(PP(:));
