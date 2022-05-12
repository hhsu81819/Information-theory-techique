function Cri=SigtestInfo(pdf_XS1,pdf_XS2,pdf_XTA,nbin,nTests)
%This fuction outputs all the parameters for information theory on three
%variables. The length of the data have to be the same. 
%Note that variable should be normalized before!
%
%X1=TS_pdf(:);
%X2=EC_pdf(:);
%Xtar=LH_pdf(:);
%{
% Below is the process to save the random ways to permute timeseries
s = rng;
r = randperm(30)
for test=1:nTests
	s=rng;
	rng(s)
	r = randperm(30)
	save(['./tmp_sig/sX1_' num2str(test) '.mat'], 's')
	s=rng;
	rng(s)
	r = randperm(30)
        save(['./tmp_sig/sX2_' num2str(test) '.mat'], 's')
	s=rng;
	rng(s);
	r = randperm(30)
        save(['./tmp_sig/sX3_' num2str(test) '.mat'], 's')
end
%}

X1_shuff=reshape(pdf_XS1,[size(pdf_XS1,1)*size(pdf_XS1,2) 1]);
X2_shuff=reshape(pdf_XS2,[size(pdf_XS1,1)*size(pdf_XS1,2) 1]);
X3_shuff=reshape(pdf_XTA,[size(pdf_XS1,1)*size(pdf_XS1,2) 1]);

InfoT=MMI_fixedbin(X1_shuff,X2_shuff,X3_shuff,nbin);
I_real=InfoT;

for test=1:nTests
	load (['/homes/hhsu/02.InfoTheo/Function/InfoTheory_v2/tmp_sig/sX1_' num2str(test) '.mat'])
	rng(s);
        ran1=randperm(size(pdf_XS1,1));
	load (['/homes/hhsu/02.InfoTheo/Function/InfoTheory_v2/tmp_sig/sX2_' num2str(test) '.mat'])
        rng(s);
        ran2=randperm(size(pdf_XS2,1));
	load (['/homes/hhsu/02.InfoTheo/Function/InfoTheory_v2/tmp_sig/sX3_' num2str(test) '.mat'])
        rng(s);
        ran3=randperm(size(pdf_XTA,1));

        X1_shuff = pdf_XS1(ran1,:);
        X2_shuff = pdf_XS2(ran2,:);
        X3_shuff = pdf_XTA(ran3,:);

        X1_shuff=reshape(X1_shuff,[size(pdf_XS1,1)*size(pdf_XS1,2) 1]);
        X2_shuff=reshape(X2_shuff,[size(pdf_XS1,1)*size(pdf_XS1,2) 1]);
        X3_shuff=reshape(X3_shuff,[size(pdf_XS1,1)*size(pdf_XS1,2) 1]);

        InfoT=MMI_fixedbin(X1_shuff,X2_shuff,X3_shuff,nbin);
        MMI_test(test)=InfoT;
end
MMI=cat(2,I_real,MMI_test);
Cri=mean(MMI(:))+3*std(MMI(:));
