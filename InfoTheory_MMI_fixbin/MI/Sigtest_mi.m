function Cri=Sigtest_mi(pdf_XS1,pdf_XTA,nbin,nTests)
%This fuction outputs all the parameters for information theory on three
%variables. The length of the data have to be th same. 
%Note that variable should be normalized before!

%X1=TS_pdf(:);
%X2=EC_pdf(:);
%Xtar=LH_pdf(:);
%{
%}
%X1=resample(X1,200,numel(X1(:)));
%X2=resample(X2,200,numel(X2(:)));
%Xtar=resample(Xtar,200,numel(Xtar(:)));
X1_shuff=reshape(pdf_XS1,[size(pdf_XS1,1)*size(pdf_XS1,2) 1]);
X3_shuff=reshape(pdf_XTA,[size(pdf_XS1,1)*size(pdf_XS1,2) 1]);

InfoT=MI(X1_shuff,X3_shuff,nbin);
I_real=InfoT;

for test=1:nTests
	load (['/homes/hhsu/02.InfoTheo/Function/InfoTheory_v2/tmp_sig/sX1_' num2str(test) '.mat'])
	rng(s);
        ran1=randperm(size(pdf_XS1,1));
	load (['/homes/hhsu/02.InfoTheo/Function/InfoTheory_v2/tmp_sig/sX3_' num2str(test) '.mat'])
        rng(s);
        ran3=randperm(size(pdf_XTA,1));

        X1_shuff = pdf_XS1(ran1,:);
        X3_shuff = pdf_XTA(ran3,:);

        X1_shuff=reshape(X1_shuff,[size(pdf_XS1,1)*size(pdf_XS1,2) 1]);
        X3_shuff=reshape(X3_shuff,[size(pdf_XS1,1)*size(pdf_XS1,2) 1]);

        InfoT=MI(X1_shuff,X3_shuff,nbin);
        MI_test(test)=InfoT;
end
MIv=cat(2,I_real,MI_test);
Cri=mean(MIv(:))+3*std(MIv(:));
