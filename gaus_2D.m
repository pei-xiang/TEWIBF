function gau2_xy=gaus_2D(i,K,idx,idy,sigma,data)
num=i;
gauxy=zeros(K,1);
ys=zeros(K,1);
for j=1:K
    gauxy(j)=1/(2*pi*sigma^2).*exp(-((idx(j)-idx(num))^2+(idy(j)-idy(num))^2)/(2*sigma^2));  %spatial filter kernel
end
d1=data(idx(num),idy(num),:);
d1=d1(:);
for k=1:K
    d2=data(idx(k),idy(k),:);
    d2=d2(:);
    ys(k)=1./(1+exp(-acos(dot(d1,d2)/(norm(d1)*norm(d2)))));      %spectral weight
end
gau2_xy=gauxy.*ys;
end