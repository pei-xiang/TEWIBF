function X=normalize(A)
[M,N,B]=size(A);
B1=reshape(A,M*N,B);
C=B1';
A1=min(C,[],2);   %�����ÿһ�е���Сֵ
A2=max(C,[],2);   %�����ÿһ�е����ֵ
X1=zeros(B,M*N);
for i=1:M*N
    X1(:,i)=(C(:,i)-A1)./(A2-A1);  %ÿһ�и�Ԫ�ؼ�ȥ��Сֵ���ٳ������ֵ����Сֵ֮��
end
X=reshape(X1',M,N,B);