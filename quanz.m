function qz=quanz(A)
[lenth,with,bands]=size(A);
sp=reshape(A,lenth*with,bands);
sp1=sp';
mn1=mean(sp1,2);
qz=zeros(lenth,with);
for is=1:lenth
    for js=1:with
        xx=A(is,js,:);
        xx=xx(:);
        qz(is,js)=acos(dot(xx,mn1)/(norm(xx)*norm(mn1)));
    end
end
end