
%%������ɫģ�� RollingGM(1,1)
%inputΪ5�Ա���,outputΪ1�����
function [PRE MRE lambda]=RollingGM(lambda,input,output)
syms a b;
c=[a b]';
A=input;
B=cumsum(A);%ԭʼ�����ۼ�
n=length(A);
for i=1:(n-1)
    C(i)=lambda*B(i)+(1-lambda)*B(i+1);%�����ۼӾ���
end
%�������������ֵ
D=A;D(1)=[];
D=D';
E=[-C;ones(1,n-1)];
c=inv(E*E')*E*D;
c=c';
a=c(1);b=c(2);
%���ݺ�������
F=[];F(1)=A(1);
for i=2:(n+2)
    F(i)=(A(1)-b/a)/exp(a*(i-1))+b/a;
end
G=[];G(1)=A(1);
for i=2:(n+2)
    G(i)=F(i)-F(i-1);%�õ�Ԥ�����������
end
PRE=G(1,1:n+1);
PRE_fit=G(1,1:n);
PRE_out=G(1,n+1);
MRE=mean(abs((PRE_fit'-input')./input'));
lambda;


