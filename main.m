clear;
clc;
close all;
tic;

load 'Sandiego.mat';
mask=map;

hsize=5;         %filter window
sigma1=0.9;      %standard deviation
numsp=150;       % number of superpixels 
sigma2=0.9;      %standard deviation

[m,n,b]=size(data);
N=m*n;
Dat=normalize(data);
dat=reshape(Dat,m*n,b)';

% Taylor expansion for HSI
k=5;
q=zeros(b,k);
DC=zeros(b,N);
C=cell(k-1,1);
D=cell(k-1,1);
for i=1:N
    pt=dat(:,i);
    q=repmat(pt,1,k);
    DC(:,i)=sin(pt)+cos(pt);  %enhanced information
    for j=2:k
        vs=ones(b,1);
        for jj=1:j
            vs=vs.*q(:,jj);
        end
        C{j-1,1}(:,i)=vs;              %j-th energy information
    end
end
dc=reshape(DC',m,n,b);     %enhanced image
for j=1:k-1
    D{j,1}=reshape(C{j,1}',m,n,b);    %energy images
end


gausbg=FPN(dc,hsize,sigma1);    %3D hyperspectral feature pyramid
gfbg=wugz(gausbg,numsp,sigma2);  %Weighted irregular block filter

fw=zw(Dat);        %local difference weight
E=cell(k-1,1);
for i=1:k-1
    qz=quanz(D{i,1});      %global adaptive weight
    jq=fw.*qz.*D{i,1};
    E{i,1}=reshape(jq,m*n,b)';
end


%inverse Taylor expansion
ec=cell(k-1,1);
ht=zeros(m,n,b);
for j=2:k
    yj=sin(j*pi/2)+cos(j*pi/2);        %j-th derivative
    wj=-yj*(factorial(j))^(-1);        %coefficient of the jth energy
    ec{j-1,1}=wj*E{j-1,1};
    ht=ht+reshape(ec{j-1,1}',m,n,b);   %sum of energy images
end
X1=gfbg+ht;       %reconstructed HSI

%anomaly target extraction
fr=sum(abs(X1).^2,3);
res=(exp(-fr)+log(fr+1)).^2;
teibf=mat2gray(res);
figure;imshow(teibf);
toc;


%disp('Running ROC...');
det_map=reshape(teibf,N,1);
GT=reshape(mask,N,1);
mode_eq=1;
[AUC_D_F,AUC_D_tau,AUC_F_tau,AUC_TD,AUC_BS,AUC_SNPR,AUC_TDBS,AUC_ODP]=plot_3DROC(det_map,GT,mode_eq);