function [S,D]=CreateIPSMatrix(xyz,xyz_aero)
%function [S,D]=CreateIPSMatrix(xyz,xyz_aero)
%   采用IPS方法计算结构-气动网格的插值矩阵
%
%   输入参数:
%       xyz : 结构网格点坐标矩阵, 每行一个点, 仅用到前两列数据
%       xyz_aero : 气动网格点坐标矩阵, 每行一个点, 仅用到前两列数据
%
%   输出参数:
%       S : 变量插值矩阵, 行数=气动网格点个数, 列数=结构网格点个数
%       D : 变量对x的导数插值矩阵, 尺寸与S相同
%
%   参考文献:
%       赵永辉,黄锐. 高等气动弹性力学与控制. 科学出版社. 2015:56-60
%

ns=size(xyz,     1);%结构网格点个数
na=size(xyz_aero,1);%气动网格点个数

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