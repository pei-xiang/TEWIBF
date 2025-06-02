function fw=zw(A)
[m,n,b]=size(A);
kdata=zeros(m+2,n+2,b);
kdata(2:m+1,2:n+1,:)=A;
kdata(1,2:n+1,:)=A(1,:,:);
kdata(m+2,2:n+1,:)=A(m,:,:);
kdata(:,1,:)=kdata(:,2,:);
kdata(:,n+2,:)=kdata(:,n-1,:);
Gx=zeros(m,n,b);
Gy=zeros(m,n,b);
Dx=zeros(m,n,b);
Dy=zeros(m,n,b);
for i=1:m
    for j=1:n
        Gx(i,j,:)=(kdata(i+1,j+2,:)-kdata(i+1,j,:))./2;
        Gy(i,j,:)=(kdata(i+2,j+1,:)-kdata(i,j+1,:))./2;
        Dx(i,j,:)=(kdata(i+2,j+2,:)-kdata(i,j,:))./2;
        Dy(i,j,:)=(kdata(i+2,j,:)-kdata(i,j+2,:))./2;
    end
end
td=Gx.^2+Gy.^2+Dx.^2+Dy.^2;
ft=sum(td.^2,3);
fw=1./(1+exp(-ft).^2);   % local difference weight
end