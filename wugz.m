function lvdata=wugz(data,numsp,sigma)
n_pc = 4;
[nl, ns, nb] = size(data); 

PC = fPCA_2D_SpecificV1(reshape(data, nl * ns, nb),1,1,0);
PC = PC(:,1:n_pc);
HIM = reshape(PC,nl,ns,n_pc);
x = PC';

% reshape the data into a vector
input_img = zeros(1, nl * ns * n_pc);

startpos = 1;
for i = 1 : nl
    for j = 1 : ns
        input_img(startpos : startpos + n_pc - 1) = HIM(i, j, :); % bands
        startpos = startpos + n_pc;
    end
end


numSuperpixels = numsp;  % the desired number of superpixels
compactness = 0.1; % compactness2 = 1-compactness, the clustering is according to: compactness*dxy+compactness2*dspectral
dist_type = 2; % distance type - 1:Euclidean£»2£ºSAD; 3:SID; 4:SAD-SID
seg_all = 1; % 1: all the pixles participate in the clustering£¬ 2: some pixels would not

[labels, numlabels, seedx, seedy] = hybridseg(input_img, nl, ns, n_pc, numSuperpixels, compactness, dist_type, seg_all);
Lab=reshape(labels,nl,ns);

% padding avoid edges 
Data=zeros(nl+2, ns+2, nb);   
Data(2:nl+1,2:ns+1,:)=data;
Data(1,2:ns+1,:)=data(1,:,:);
Data(nl+2,2:ns+1,:)=data(nl,:,:);
Data(:,1,:)=Data(:,2,:);
Data(:,ns+2,:)=Data(:,ns+1,:);

lvdata=zeros(nl, ns, nb);
for num=1:numlabels
    [idx,idy]=find(Lab==num);    %Find the address for each superpixel
    K=size(idx,1);
    ZK=K;
    Idx=idx;
    Idy=idy;
    idx=idx+1;
    idy=idy+1;
    for i=1:K
        gau2_xy=gaus_2D(i,ZK,idx,idy,sigma,Data);  %filter kernel
        tep=0;
        for kk=1:K
            tep=tep+Data(idx(kk),idy(kk),:).*gau2_xy(kk);
        end       
        lvdata(Idx(i),Idy(i),:)=tep;
    end
end
