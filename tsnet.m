D=importdata('data.csv');
P=importdata('para.csv');
tsner=tsne(D.data,'Algorithm','exact','Perplexity',P.data(1),'NumDimensions',P.data(2),'NumPCAComponents',P.data(2));
csvwrite('tsnetoutput.csv',tsner);
clear all