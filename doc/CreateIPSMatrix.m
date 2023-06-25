function [S,D]=CreateIPSMatrix(xyz,xyz_aero)
%function [S,D]=CreateIPSMatrix(xyz,xyz_aero)
%   ����IPS��������ṹ-��������Ĳ�ֵ����
%
%   �������:
%       xyz : �ṹ������������, ÿ��һ����, ���õ�ǰ��������
%       xyz_aero : ����������������, ÿ��һ����, ���õ�ǰ��������
%
%   �������:
%       S : ������ֵ����, ����=������������, ����=�ṹ��������
%       D : ������x�ĵ�����ֵ����, �ߴ���S��ͬ
%
%   �ο�����:
%       ������,����. �ߵ�����������ѧ�����. ��ѧ������. 2015:56-60
%

ns=size(xyz,     1);%�ṹ��������
na=size(xyz_aero,1);%������������

[xi,xj]=meshgrid(xyz(:,1),xyz(:,1));
[yi,yj]=meshgrid(xyz(:,2),xyz(:,2));

rij2=(xi-xj).^2+(yi-yj).^2;
K=rij2.*log(rij2+eps);
R=[ones(ns,1),xyz(:,1:2)];

%C=[zeros(3),R';R,K];
C=[K,R;R',zeros(3)];
Cinv=inv(C);

[xk,xj]=meshgrid(xyz_aero(:,1),xyz(:,1));
[yk,yj]=meshgrid(xyz_aero(:,2),xyz(:,2));
xk=xk';xj=xj';
yk=yk';yj=yj';

rkj2=(xk-xj).^2+(yk-yj).^2;
Kk=rkj2.*log(rkj2+eps);
%Kbar=[ones(na,1),xyz_aero(:,1:2),Kk];
Kbar=[Kk,ones(na,1),xyz_aero(:,1:2)];
Kbc=Kbar/C;
%S=KbC(:,4:end);
S=Kbc(:,1:ns);

Dk=2*(xk-xj).*(1+log(rkj2+eps));
%Khat=[zeros(na,1),ones(na,1),zeros(na,1),Dk];
Khat=[Dk,zeros(na,1),ones(na,1),zeros(na,1)];
KhC=Khat/C;

%D=KhC(:,4:end);
D=KhC(:,1:ns);

end