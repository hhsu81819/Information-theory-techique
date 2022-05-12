function MMI = MMI(pxyz,pxy,p3)

%Writer: Hsin Hsu
bins=size(pxy,1);
ppe=repmat(p3,[1 bins bins]);
ppe=permute(ppe,[2 3 1]);
ppxy=repmat(pxy,[1 1 bins]);
LL=(pxyz./(ppxy.*ppe));
PP=pxyz.*log2(LL);
MMI=nansum(PP(:));

