%% Paper5主梁位移计算代码（理论力学法）-喆泓大桥模型验证
% 吊杆悬链线
% 是以Py为基本未知量

warning off;       % 关闭报警信息
close all;         % 关闭开启的图窗
clear all;         % 清空变量
clc;               % 清空命令行

% load 'x_xi.txt';
% x=x_xi;
% load 'x_deltHanger.txt';
% x=x_deltHanger;
% load 'x_omega.txt';
% x=x_omega;
load 'x_Case4.txt';
x=x_Case4;
%------------------全桥参数（永宗大桥）------------------
% 主缆
Ec=2.05e8;      Ac=0.0715;     qc=6.2;
HA=29.0197;    HB=97.7050;       HC=99.9150;   HD=29.0197;     H_Mid=41.0320;
YA=22.4409;    YB=1;       YC=1;   YD=22.4409;   % 这一行2D时没有

% 吊杆
Eh=2.06e8;   Ah=0.00609;   qh=0.468981765;
nL_hanger=7;  nM_hanger=27;  nR_hanger=7;
gamma_c=0;  % 索夹自重

% 主梁
Eg=2.06e8;     Ig=4;    Ag=1.12;  Ag_JM=1.12;
Hg=22.4409*2;  n_zhizuo=4;

% 建模所用
Ig_xx=9.0987;  Ig_zz=4;  Ig_yy=207.04;  h_zz=36;  h_yy=3.5;  theta_x=0;
qg=200*ones(nL_hanger+nM_hanger+nR_hanger+n_zhizuo,1);

% 导入已知量
load 'xP.txt'; % xP是导入的吊点距主梁左端点的距离
load 'Coef_Eta4.txt'; % 主梁多项式系数
Coef_Eta=Coef_Eta4;

% 基本未知量
aM(1,:)=x(1); HM(1,:)=x(2); alpha_M(1,:)=x(3);
aL(1,:)=x(4); HL(1,:)=x(5); alpha_L(1,:)=x(6);
aR(1,:)=x(7); HR(1,:)=x(8); alpha_R(1,:)=x(9);
Hx=2*x(10);
for i=1:nL_hanger+nM_hanger+nR_hanger+n_zhizuo     % 主梁上的Py，其中有4个是支座反力，输入值为正。【x中11-53】
Py(i,:)=x(10+i); 
end
for i=1:nL_hanger 
ahL(i,:)=x(10+nL_hanger+nM_hanger+nR_hanger+n_zhizuo+i);    %【x中54-61】
end
for i=1:nM_hanger 
ahM(i,:)=x(10+2*nL_hanger+nM_hanger+nR_hanger+n_zhizuo+i);   %【x中62-84】
end
for i=1:nR_hanger 
ahR(i,:)=x(10+2*nL_hanger+2*nM_hanger+nR_hanger+n_zhizuo+i);    %【x中85-92】
end

% 主梁Z的计算
for i=1:nL_hanger+nM_hanger+nR_hanger+n_zhizuo
    ZG(i,:)=Coef_Eta(i,1)+Coef_Eta(i,2)*xP(i)+Coef_Eta(i,3)*xP(i)^2/2+Coef_Eta(i,4)*xP(i)^3/6+Coef_Eta(i,5)*xP(i)^4/24;
end
YGL=(YA-YB)*ones(nL_hanger,1); 
YGM=(YA-YB)*ones(nM_hanger,1);
YGR=(YA-YB)*ones(nR_hanger,1);

% 主缆悬链线段水平投影长度l的计算
for i=1:nL_hanger+nM_hanger+nR_hanger+n_zhizuo-1
    l(i,1)=xP(i+1)-xP(i);
end
lL=l(1:nL_hanger+1,1);
lM=l(nL_hanger+2:nL_hanger+nM_hanger+2,1);
lR=l(nL_hanger+nM_hanger+3:nL_hanger+nM_hanger+nR_hanger+3,1);

% 将吊杆力分开表示（用于主缆，论文中的P'）
for i=1:nL_hanger
    PLy(i,:)=Py(i+1,:);   % 用于主缆、应该>0
    ZGL(i,:)=ZG(i+1,:);   % 根据输入的主梁预拱度多项式计算得到
end
for i=1:nM_hanger
    PMy(i,:)=Py(i+nL_hanger+2,:);   % 用于主缆、应该>0
    ZGM(i,:)=ZG(i+nL_hanger+2,:);
end
for i=1:nR_hanger
    PRy(i,:)=Py(i+nL_hanger+nM_hanger+3,:);   % 用于主缆、应该>0
    ZGR(i,:)=ZG(i+nL_hanger+nM_hanger+3,:);
end

% 悬链线-左边跨
for i=1:nL_hanger+1
    % 主缆
    cL(i,:)=-HL(i,:)/qc; % 应该<0
    deltyL(i,:)=lL(i,:)*tan(alpha_L(i,:));   % 应该<0,数字先小后大
    deltzL(i,:)=cL(i,:)*(cosh(lL(i,:)/cL(i,:)/cos(alpha_L(i,:))+aL(i,:))-cosh(aL(i,:)));   % 应该<0,数字先小后大
    LL(i,:)=cL(i,:)*(sinh(lL(i,:)/cL(i,:)/cos(alpha_L(i,:))+aL(i,:))-sinh(aL(i,:)));
    SL(i,:)=cL(i,:)*(sinh(lL(i,:)/cL(i,:)/cos(alpha_L(i,:))+aL(i,:))-sinh(aL(i,:)))-HL(i,:)/2/Ec/Ac*(lL(i,:)/cos(alpha_L(i,:))+cL(i,:)/2*(sinh(2*(lL(i,:)/cL(i,:)/cos(alpha_L(i,:))+aL(i,:)))-sinh(2*aL(i,:)))); % 应该>0
    j=i;
    if j<=nL_hanger
        % 吊杆
        chL(j,:)=-PLy(j,:)/qh; % 应该<0
        sum_L(j,:)=sum(deltyL(1:j,:));  % 应该<0,左跨的主缆y向投影长度就是吊杆的长度
        lhL(j,:)=-sum_L(j,:);  % 应该>0
        deltzhL(j,:)=chL(j,:)*(cosh(lhL(j,:)/chL(j,:)+ahL(j,:))-cosh(ahL(j,:)));    % 应该>0,数字先小后大
        PLz(j,:)=PLy(j,:)*sinh(lhL(j,:)/chL(j,:)+ahL(j,:));   % 应该>0,利用吊杆悬链线下吊点夹角
        LhL(j,:)=chL(j,:)*(sinh(lhL(j,:)/chL(j,:)+ahL(j,:))-sinh(ahL(j,:)));
        ShL(j,:)=chL(j,:)*(sinh(lhL(j,:)/chL(j,:)+ahL(j,:))-sinh(ahL(j,:)))-PLy(j,:)/2/Eh/Ah*(lhL(j,:)+chL(j,:)/2*(sinh(2*(lhL(j,:)/chL(j,:)+ahL(j,:)))-sinh(2*ahL(j,:))));  % 应该>0
        gamma_hL(j,:)=qh*ShL(j,:);  % 应该>0
        PLz_dot2(j,:)=PLz(j,:)+gamma_hL(j,:)+gamma_c;  % 应该>0
        PLz2(j,:)=PLz_dot2(j,:)-gamma_hL(j,:)-gamma_c;
        deltzhL2(j,:)=-sum(deltzL(1:j,:))-ZGL(j,:);      % 应该>0,吊杆在z轴上的投影长度    deltzhL2(j,:)=-sum(deltzL(1:j,:))-ZGL(j,:); 
        PLz_dot(j,:)=PLy(j,:)*sinh(ahL(j,:));   %  应该>0,利用吊杆悬链线上吊点夹角
        % 主缆悬链线参数递推公式  
        alpha_L(j+1,:)=atan((HL(j,:)*sin(alpha_L(j,:))-PLy(j,:))/(HL(j,:)*cos(alpha_L(j,:))));  % 应该<0
        HL(j+1,:)=HL(j,:)*cos(alpha_L(j,:))*sec(alpha_L(j+1,:));  % 应该>0
        aL(j+1,:)=asinh((HL(j,:)*sinh(lL(j,:)/cL(j,:)/cos(alpha_L(j,:))+aL(j,:))-PLz_dot(j,:))/HL(j+1,:));  % 应该<0
    else
        break
    end
end
sum_deltzL=sum(deltzL);  % 应该<0
sum_deltyL=sum(deltyL);  % 应该<0

% 悬链线-中跨
for i=1:nM_hanger+1
    % 主缆
    cM(i,:)=-HM(i,:)/qc; % 应该<0
    deltyM(i,:)=lM(i,:)*tan(alpha_M(i,:));   % 先正后负（先大后小）且应该对称
    deltzM(i,:)=cM(i,:)*(cosh(lM(i,:)/cM(i,:)/cos(alpha_M(i,:))+aM(i,:))-cosh(aM(i,:)));   % 先正后负
    LM(i,:)=cM(i,:)*(sinh(lM(i,:)/cM(i,:)/cos(alpha_M(i,:))+aM(i,:))-sinh(aM(i,:)));
    SM(i,:)=cM(i,:)*(sinh(lM(i,:)/cM(i,:)/cos(alpha_M(i,:))+aM(i,:))-sinh(aM(i,:)))-HM(i,:)/2/Ec/Ac*(lM(i,:)/cos(alpha_M(i,:))+cM(i,:)/2*(sinh(2*(lM(i,:)/cM(i,:)/cos(alpha_M(i,:))+aM(i,:)))-sinh(2*aM(i,:))));  % 应该>0
    j=i;
    if j<=nM_hanger 
        % 吊杆
        chM(j,:)=-PMy(j,:)/qh; % 应该<0
        lhM(j,:)=YGM(j,:)-sum(deltyM(1:j,:));   % 应该>0
        deltzhM(j,:)=chM(j,:)*(cosh(lhM(j,:)/chM(j,:)+ahM(j,:))-cosh(ahM(j,:)));  % 应该>0
        PMz(j,:)=PMy(j,:)*sinh(lhM(j,:)/chM(j,:)+ahM(j,:));   % 应该>0，利用吊杆悬链线下吊点夹角
        LhM(j,:)=chM(j,:)*(sinh(lhM(j,:)/chM(j,:)+ahM(j,:))-sinh(ahM(j,:)));
        ShM(j,:)=chM(j,:)*(sinh(lhM(j,:)/chM(j,:)+ahM(j,:))-sinh(ahM(j,:)))-PMy(j,:)/2/Eh/Ah*(lhM(j,:)+chM(j,:)/2*(sinh(2*(lhM(j,:)/chM(j,:)+ahM(j,:)))-sinh(2*ahM(j,:))));  % 应该>0
        gamma_hM(j,:)=qh*ShM(j,:);  % 应该>0
        PMz_dot2(j,:)=PMz(j,:)+gamma_hM(j,:)+gamma_c;   % 应该>0
        PMz2(j,:)=PMz_dot2(j,:)-gamma_hM(j,:)-gamma_c;
        deltzhM2(j,:)=HB-sum(deltzM(1:j,:))-ZGM(j,:);   % 应该>0   deltzhM2(j,:)=HB-sum(deltzM(1:j,:))-ZGM(j,:);
        PMz_dot(j,:)=PMy(j,:)*sinh(ahM(j,:));    % 应该>0，利用吊杆悬链线上吊点夹角
        % 主缆悬链线参数递推公式
        alpha_M(j+1,:)=atan((HM(j,:)*sin(alpha_M(j,:))-PMy(j,:))/(HM(j,:)*cos(alpha_M(j,:))));   % 应该先正后负
        HM(j+1,:)=HM(j,:)*cos(alpha_M(j,:))*sec(alpha_M(j+1,:));   % 应该>0
        aM(j+1,:)=asinh((HM(j,:)*sinh(lM(j,:)/cM(j,:)/cos(alpha_M(j,:))+aM(j,:))-PMz_dot(j,:))/HM(j+1,:));  % 应该先正后负
    else
        break
    end
end
sum_deltzM=sum(deltzM);  % 应该=0，zB-zC
sum_deltzM_middle=sum(deltzM(1:ceil(nM_hanger/2),:));  % 应该>0， zB-zM
sum_deltyM=sum(deltyM);  % 应该=0，yB-yC

% 悬链线-右边跨
for i=1:nR_hanger+1
    % 主缆
    cR(i,:)=-HR(i,:)/qc; % 应该<0
    deltyR(i,:)=lR(i,:)*tan(alpha_R(i,:));   % 应该>0
    deltzR(i,:)=cR(i,:)*(cosh(lR(i,:)/cR(i,:)/cos(alpha_R(i,:))+aR(i,:))-cosh(aR(i,:)));  % 应该>0
    LR(i,:)=cR(i,:)*(sinh(lR(i,:)/cR(i,:)/cos(alpha_R(i,:))+aR(i,:))-sinh(aR(i,:)));
    SR(i,:)=cR(i,:)*(sinh(lR(i,:)/cR(i,:)/cos(alpha_R(i,:))+aR(i,:))-sinh(aR(i,:)))-HR(i,:)/2/Ec/Ac*(lR(i,:)/cos(alpha_R(i,:))+cR(i,:)/2*(sinh(2*(lR(i,:)/cR(i,:)/cos(alpha_R(i,:))+aR(i,:)))-sinh(2*aR(i,:))));   % 应该>0
    j=i;
    if j<=nR_hanger
        % 吊杆
        chR(j,:)=-PRy(j,:)/qh; % 应该<0
        lhR(j,:)=YGR(j,:)-sum(deltyR(1:j,:));   % 应该>0,数字先大后小
        deltzhR(j,:)=chR(j,:)*(cosh(lhR(j,:)/chR(j,:)+ahR(j,:))-cosh(ahR(j,:)));   % 应该>0,数字先大后小
        PRz(j,:)=PRy(j,:)*sinh(lhR(j,:)/chR(j,:)+ahR(j,:));   % 应该>0，利用吊杆悬链线下吊点夹角
        LhR(j,:)=chR(j,:)*(sinh(lhR(j,:)/chR(j,:)+ahR(j,:))-sinh(ahR(j,:)));
        ShR(j,:)=chR(j,:)*(sinh(lhR(j,:)/chR(j,:)+ahR(j,:))-sinh(ahR(j,:)))-PRy(j,:)/2/Eh/Ah*(lhR(j,:)+chR(j,:)/2*(sinh(2*(lhR(j,:)/chR(j,:)+ahR(j,:)))-sinh(2*ahR(j,:))));  % 应该>0
        gamma_hR(j,:)=qh*ShR(j,:);  % 应该>0，
        PRz_dot2(j,:)=PRz(j,:)+gamma_hR(j,:)+gamma_c;  % 应该>0
        PRz2(j,:)=PRz_dot2(j,:)-gamma_hR(j,:)-gamma_c;
        deltzhR2(j,:)=HC-sum(deltzR(1:j,:))-ZGR(j,:);  % 应该>0   deltzhR2(j,:)=HC-sum(deltzR(1:j,:))-ZGR(j,:); 
        PRz_dot(j,:)=PRy(j,:)*sinh(ahR(j,:));          % 应该>0
        % 主缆
        alpha_R(j+1,:)=atan((HR(j,:)*sin(alpha_R(j,:))-PRy(j,:))/(HR(j,:)*cos(alpha_R(j,:))));     % 应该>0
        HR(j+1,:)=HR(j,:)*cos(alpha_R(j,:))*sec(alpha_R(j+1,:));     % 应该>0
        aR(j+1,:)=asinh((HR(j,:)*sinh(lR(j,:)/cR(j,:)/cos(alpha_R(j,:))+aR(j,:))-PRz_dot(j,:))/HR(j+1,:));     % 应该>0   
    else
        break
    end
end
sum_deltzR=sum(deltzR);  % zC-zD
sum_deltyR=sum(deltyR);  % yC-yD

% 主缆坐标(以左锚点为全局坐标系原点)
XO=xP;
YO(1,:)=0+Hg/2;
for i=1:nL_hanger+1
    YO(i+1,:)=YO(1,:)+sum(deltyL(1:i,:));
end
for i=1:nM_hanger+1
    YO(i+nL_hanger+2,:)=YO(nL_hanger+2,:)+sum(deltyM(1:i,:));
end
for i=1:nR_hanger+1
    YO(i+nL_hanger+nM_hanger+3,:)=YO(nL_hanger+nM_hanger+3,:)+sum(deltyR(1:i,:));
end
ZO(1,:)=HA;
for i=1:nL_hanger+1
    ZO(i+1,:)=ZO(1,:)-sum(deltzL(1:i,:));
end
for i=1:nM_hanger+1
    ZO(i+nL_hanger+2,:)=ZO(nL_hanger+2,:)-sum(deltzM(1:i,:));
end
for i=1:nR_hanger+1
    ZO(i+nL_hanger+nM_hanger+3,:)=ZO(nL_hanger+nM_hanger+3,:)-sum(deltzR(1:i,:));
end

% 吊杆建模坐标---采用细分式建模
% 左跨
nz_hL=30;
for i=1:nL_hanger    
    for j=1:nz_hL
        lhL_1(i,j)=j*lhL(i,:)/nz_hL;
        zhL(i,j)=chL(i,:)*(cosh(lhL_1(i,j)/chL(i,:)+ahL(i,:))-cosh(ahL(i,:)));
        LhL_1(i,j)=chL(i,:)*(sinh(lhL_1(i,j)/chL(i,:)+ahL(i,:))-sinh(ahL(i,:)));
        ShL_1(i,j)=chL(i,:)*(sinh(lhL_1(i,j)/chL(i,:)+ahL(i,:))-sinh(ahL(i,:)))-PLy(i,:)/2/Eh/Ah*(lhL_1(i,j)+chL(i,:)/2*(sinh(2*(lhL_1(i,j)/chL(i,:)+ahL(i,:)))-sinh(2*ahL(i,:))));  % 应该>0
        LhL_2(i,1)=LhL_1(i,1);
        ShL_2(i,1)=ShL_1(i,1);
        if j>=2
            LhL_2(i,j)=LhL_1(i,j)-LhL_1(i,j-1);
            ShL_2(i,j)=ShL_1(i,j)-ShL_1(i,j-1);
        end
        epsilon_hanger_L(j,i)=(LhL_2(i,j)-ShL_2(i,j))/LhL_2(i,j);
        X_hL(j,i)=XO(i+1,:); % (这个是一列是一根吊杆上的代码)             % X_hL(i,j)=XO(i+1,:);  (这个是一行是一根吊杆上的代码)
        Z_hL(j,i)=ZO(i+1,:)-zhL(i,j);                                  % Z_hL(i,j)=ZO(i+1,:)-zhL(i,j);
        Y_hL(j,i)=YO(i+1,:)+lhL_1(i,j);                                % Y_hL(i,j)=YO(i+1,:)+lhL_1(i,j);
    end
end
% 中跨
nz_hM=30;
for i=1:nM_hanger    
    for j=1:nz_hM
        lhM_1(i,j)=j*lhM(i,:)/nz_hM;
        zhM(i,j)=chM(i,:)*(cosh(lhM_1(i,j)/chM(i,:)+ahM(i,:))-cosh(ahM(i,:)));
        LhM_1(i,j)=chM(i,:)*(sinh(lhM_1(i,j)/chM(i,:)+ahM(i,:))-sinh(ahM(i,:)));
        ShM_1(i,j)=chM(i,:)*(sinh(lhM_1(i,j)/chM(i,:)+ahM(i,:))-sinh(ahM(i,:)))-PMy(i,:)/2/Eh/Ah*(lhM_1(i,j)+chM(i,:)/2*(sinh(2*(lhM_1(i,j)/chM(i,:)+ahM(i,:)))-sinh(2*ahM(i,:))));  % 应该>0
        LhM_2(i,1)=LhM_1(i,1);
        ShM_2(i,1)=ShM_1(i,1);
        if j>=2
            LhM_2(i,j)=LhM_1(i,j)-LhM_1(i,j-1);
            ShM_2(i,j)=ShM_1(i,j)-ShM_1(i,j-1);
        end
        epsilon_hanger_M(j,i)=(LhM_2(i,j)-ShM_2(i,j))/LhM_2(i,j);  
        X_hM(j,i)=XO(i+nL_hanger+2,:);
        Z_hM(j,i)=ZO(i+nL_hanger+2,:)-zhM(i,j);
        Y_hM(j,i)=YO(i+nL_hanger+2,:)+lhM_1(i,j);
    end
end
% 右跨
nz_hR=30;
for i=1:nR_hanger    
    for j=1:nz_hR
        lhR_1(i,j)=j*lhR(i,:)/nz_hR;
        zhR(i,j)=chR(i,:)*(cosh(lhR_1(i,j)/chR(i,:)+ahR(i,:))-cosh(ahR(i,:)));
        LhR_1(i,j)=chR(i,:)*(sinh(lhR_1(i,j)/chR(i,:)+ahR(i,:))-sinh(ahR(i,:)));
        ShR_1(i,j)=chR(i,:)*(sinh(lhR_1(i,j)/chR(i,:)+ahR(i,:))-sinh(ahR(i,:)))-PRy(i,:)/2/Eh/Ah*(lhR_1(i,j)+chR(i,:)/2*(sinh(2*(lhR_1(i,j)/chR(i,:)+ahR(i,:)))-sinh(2*ahR(i,:))));  % 应该>0
        LhR_2(i,1)=LhR_1(i,1);
        ShR_2(i,1)=ShR_1(i,1);
        if j>=2
            LhR_2(i,j)=LhR_1(i,j)-LhR_1(i,j-1);
            ShR_2(i,j)=ShR_1(i,j)-ShR_1(i,j-1);
        end
        epsilon_hanger_R(j,i)=(LhR_2(i,j)-ShR_2(i,j))/LhR_2(i,j);
        X_hR(j,i)=XO(i+nL_hanger+nM_hanger+3,:);
        Z_hR(j,i)=ZO(i+nL_hanger+nM_hanger+3,:)-zhR(i,j);   % 按理说，第10行的Z应该=主梁对应坐标ZG
        Y_hR(j,i)=YO(i+nL_hanger+nM_hanger+3,:)+lhR_1(i,j);  % 按理说，第10行的Y应该都等于22.409
    end
end



%---------------------------------------------------------------------------------主梁(解析式不同）---------------------------------------------------------------------------------
% % 有曲线主梁段长度
% for i=1:nL_hanger+nM_hanger+nR_hanger+n_zhizuo-1
%     y_dot = @(x) sqrt(1+(Coef_Eta(i,2)+Coef_Eta(i,3)*x).^2);  
%     l_bend(i,:) = integral(y_dot, xP(i,:), xP(i+1,:));
%     l_bend_sum(i+1,:) = sum(l_bend(1:i,:));
%     Gg(i,:) = qg(i,:)*l_bend(i,:);
%     qg(i,:) = Gg(i,:)/(xP(i+1,:)-xP(i,:));
% end
% l_bend_sum(1,:)=0;
% Lg=sum(l_bend);
% Gg_all=sum(Gg);

% 梁段用线段（两点之间的距离公式）
for i=1:nL_hanger+nM_hanger+nR_hanger+n_zhizuo-1 
    l_bend(i,:) = sqrt((xP(i+1,:)-xP(i,:))^2+(ZG(i+1,:)-ZG(i,:))^2);
    l_bend_sum(i+1,:) = sum(l_bend(1:i,:));
    Gg(i,:) = qg(i,:)*l_bend(i,:);
    qg(i,:) = Gg(i,:)/(xP(i+1,:)-xP(i,:));
end
l_bend_sum(1,:)=0;
Lg=sum(l_bend);
Gg_all=sum(Gg);

% 积分常数关系
% 支反力（4个）
Pz(1,:)=2*Py(1,:);
Pz(nL_hanger+2,:)=2*Py(nL_hanger+2,:);
Pz(nL_hanger+nM_hanger+3,:)=2*Py(nL_hanger+nM_hanger+3,:);
Pz(nL_hanger+nM_hanger+nR_hanger+4,:)=2*Py(nL_hanger+nM_hanger+nR_hanger+4,:);
for i=1:nL_hanger
    Pz(i+1,:)=2*PLz(i,:);                        % 与主缆上的力方向相反
end
for i=1:nM_hanger
    Pz(i+nL_hanger+2,:)=2*PMz(i,:);              % 与主缆上的力方向相反
end
for i=1:nR_hanger
    Pz(i+nL_hanger+nM_hanger+3,:)=2*PRz(i,:);    % 与主缆上的力方向相反
end
% % ---------------------------------------------验证弯矩时所用---------------------------------------------
% load 'Pz_dan.txt';
% Pz=2*Pz_dan;
% % qg(:,:)=445.180;
%
% % ---------------------------------------------验证弯矩时所用---------------------------------------------

%% 用于检验弯矩、转角、挠度公式是否正确  【直梁段长度计算自重+全局坐标系】
sum_C1_1=0;
sum_C1_2=0;
sum_C1_3=0;
sum_C1_4=0;
sum_C1_5=0;
sum_C1_6=0;
for i=1:nL_hanger+nM_hanger+nR_hanger+n_zhizuo  %-1
    sum_C1_1=sum_C1_1+Pz(i,:)*(1/6*xP(end,:)^2-1/2*xP(end,:)*xP(i,:));
end
for i=1:nL_hanger+nM_hanger+nR_hanger+n_zhizuo-1  % -2
    sum_C1_2=sum_C1_2+1/2*Pz(i+1,:)*xP(i+1,:)^2;
    sum_C1_3=sum_C1_3+1/6/xP(end,:)*Pz(i+1,:)*xP(i+1,:)^3;
    sum_C1_4=sum_C1_4+Hx/xP(end,:)*((Coef_Eta(i+1,1)*xP(i+1,:)^2/2+Coef_Eta(i+1,2)*xP(i+1,:)^3/3+Coef_Eta(i+1,3)*xP(i+1,:)^4/8+Coef_Eta(i+1,4)*xP(i+1,:)^5/30+Coef_Eta(i+1,5)*xP(i+1,:)^6/144)-(Coef_Eta(i,1)*xP(i+1,:)^2/2+...
        Coef_Eta(i,2)*xP(i+1,:)^3/3+Coef_Eta(i,3)*xP(i+1,:)^4/8+Coef_Eta(i,4)*xP(i+1,:)^5/30+Coef_Eta(i,5)*xP(i+1,:)^6/144));
    sum_C1_5=sum_C1_5+Hx*((Coef_Eta(i+1,1)*xP(i+1,:)+Coef_Eta(i+1,2)*xP(i+1,:)^2/2+Coef_Eta(i+1,3)*xP(i+1,:)^3/6+Coef_Eta(i+1,4)*xP(i+1,:)^4/24+Coef_Eta(i+1,5)*xP(i+1,:)^5/120)-(Coef_Eta(i,1)*xP(i+1,:)+...
        Coef_Eta(i,2)*xP(i+1,:)^2/2+Coef_Eta(i,3)*xP(i+1,:)^3/6+Coef_Eta(i,4)*xP(i+1,:)^4/24+Coef_Eta(i,5)*xP(i+1,:)^5/120));
    sum_C1_6=sum_C1_6+qg(i,:)*(xP(i+1,:)-xP(i,:))*((xP(i+1,:)^3+xP(i+1,:)^2*xP(i,:)+xP(i,:)^2*xP(i+1,:)+xP(i,:)^3)/24/xP(end,:)-(1/6*xP(end,:)^2-xP(end,:)*(xP(i+1,:)+xP(i,:))/4)-(xP(i+1,:)^2+xP(i,:)*xP(i+1,:)+xP(i,:)^2)/6);
end
C(1,:)=qg(nL_hanger+nM_hanger+nR_hanger+n_zhizuo,:)*(xP(end,:)-xP(nL_hanger+nM_hanger+nR_hanger+n_zhizuo,:))^4/24/xP(end,:)+Hx*(Coef_Eta(nL_hanger+nM_hanger+nR_hanger+n_zhizuo,1)/2*xP(end,:) ...
    +Coef_Eta(nL_hanger+nM_hanger+nR_hanger+n_zhizuo,2)/6*xP(end,:)^2+Coef_Eta(nL_hanger+nM_hanger+nR_hanger+n_zhizuo,3)/24*xP(end,:)^3+Coef_Eta(nL_hanger+nM_hanger+nR_hanger+n_zhizuo,4)/120*xP(end,:)^4+Coef_Eta(nL_hanger+nM_hanger+nR_hanger+n_zhizuo,5)/720*xP(end,:)^5) ...
    -sum_C1_1-sum_C1_2+sum_C1_3+sum_C1_4-sum_C1_5-sum_C1_6;

D(1,:)=0;
for i=1:nL_hanger+nM_hanger+nR_hanger+n_zhizuo-1  % -2
    % C(i+1,:)=C(i,:)+1/2*Pz(i+1,:)*xP(i+1,:)^2+Hx*((Coef_Eta(i+1,1)*xP(i+1,:)+Coef_Eta(i+1,2)*xP(i+1,:)^2/2+Coef_Eta(i+1,3)*xP(i+1,:)^3/6+Coef_Eta(i+1,4)*xP(i+1,:)^4/24)...
    %     -(Coef_Eta(i,1)*xP(i+1,:)+Coef_Eta(i,2)*xP(i+1,:)^2/2+Coef_Eta(i,3)*xP(i+1,:)^3/6+Coef_Eta(i,4)*xP(i+1,:)^4/24))-qg(i,:)*(xP(i+1,:)-xP(i,:))*(xP(i+1,:)^2+xP(i+1,:)*xP(i,:)+xP(i,:)^2)/6;
    % D(i+1,:)=D(i,:)-1/6*Pz(i+1,:)*xP(i+1,:)^3-Hx*((Coef_Eta(i+1,1)*xP(i+1,:)^2/2+Coef_Eta(i+1,2)*xP(i+1,:)^3/3+Coef_Eta(i+1,3)*xP(i+1,:)^4/8+Coef_Eta(i+1,4)*xP(i+1,:)^5/30)...
    %     -(Coef_Eta(i,1)*xP(i+1,:)^2/2+Coef_Eta(i,2)*xP(i+1,:)^3/3+Coef_Eta(i,3)*xP(i+1,:)^4/8+Coef_Eta(i,4)*xP(i+1,:)^5/30))+qg(i,:)*(xP(i+1,:)-xP(i,:))*(xP(i+1,:)^3+xP(i+1,:)^2*xP(i,:)+xP(i+1,:)*xP(i,:)^2-xP(i,:)^3)/24;
    sum_C_1=0;
    sum_C_2=0;
    sum_C_3=0;
    sum_D_1=0;
    sum_D_2=0;
    sum_D_3=0;
    for j=1:i
        sum_C_1=sum_C_1+1/2*Pz(j+1,:)*xP(j+1,:)^2;
        sum_C_2=sum_C_2+Hx*((Coef_Eta(j+1,1)*xP(j+1,:)+Coef_Eta(j+1,2)*xP(j+1,:)^2/2+Coef_Eta(j+1,3)*xP(j+1,:)^3/6+Coef_Eta(j+1,4)*xP(j+1,:)^4/24+Coef_Eta(j+1,5)*xP(j+1,:)^5/120) ...
               -(Coef_Eta(j,1)*xP(j+1,:)+Coef_Eta(j,2)*xP(j+1,:)^2/2+Coef_Eta(j,3)*xP(j+1,:)^3/6+Coef_Eta(j,4)*xP(j+1,:)^4/24+Coef_Eta(j,5)*xP(j+1,:)^5/120));
        sum_C_3=sum_C_3+qg(j,:)*(xP(j+1,:)-xP(j,:))*(xP(j+1,:)^2+xP(j,:)*xP(j+1,:)+xP(j,:)^2)/6;

        sum_D_1=sum_D_1+1/6*Pz(j+1,:)*xP(j+1,:)^3;
        sum_D_2=sum_D_2+Hx*((Coef_Eta(j+1,1)*xP(j+1,:)^2/2+Coef_Eta(j+1,2)*xP(j+1,:)^3/3+Coef_Eta(j+1,3)*xP(j+1,:)^4/8+Coef_Eta(j+1,4)*xP(j+1,:)^5/30+Coef_Eta(j+1,5)*xP(j+1,:)^6/144) ...
            -(Coef_Eta(j,1)*xP(j+1,:)^2/2+Coef_Eta(j,2)*xP(j+1,:)^3/3+Coef_Eta(j,3)*xP(j+1,:)^4/8+Coef_Eta(j,4)*xP(j+1,:)^5/30+Coef_Eta(j,5)*xP(j+1,:)^6/144));
        sum_D_3=sum_D_3+qg(j,:)*(xP(j+1,:)-xP(j,:))*(xP(j+1,:)^3+xP(j+1,:)^2*xP(j,:)+xP(j+1,:)*xP(j,:)^2+xP(j,:)^3)/24;
    end
    C(i+1,:)=C(1,:)+sum_C_1+sum_C_2-sum_C_3;
    D(i+1,:)=D(1,:)-sum_D_1-sum_D_2+sum_D_3;
end
for i=1:nL_hanger+nM_hanger+nR_hanger+n_zhizuo %-1 
    sum_PM=0;
    sum_Ptheta=0;
    sum_Pomega=0;
    for j=1:i
        sum_PM=sum_PM+Pz(j,:)*(xP(i,:)-xP(j,:));
        sum_Ptheta=sum_Ptheta+Pz(j,:)*(1/2*xP(i,:)^2-xP(i,:)*xP(j,:));
        sum_Pomega=sum_Pomega+Pz(j,:)*(1/6*xP(i,:)^3-1/2*xP(i,:)^2*xP(j,:));
    end
    sum_g_M(1,:)=0;
    sum_g_theta(1,:)=0;
    sum_g_omega(1,:)=0;
    sum_g_M1=0;
    sum_g_theta1=0;
    sum_g_omega1=0;
    for j=2:i
        sum_g_M1=sum_g_M1+qg(j-1,:)*(xP(j,:)-xP(j-1,:))*(xP(i,:)-(xP(j,:)+xP(j-1,:))/2);
        sum_g_theta1=sum_g_theta1+qg(j-1,:)*(xP(j,:)-xP(j-1,:))*(1/2*xP(i,:)^2-xP(i,:)*(xP(j,:)+xP(j-1,:))/2);
        sum_g_omega1=sum_g_omega1+qg(j-1,:)*(xP(j,:)-xP(j-1,:))*(1/6*xP(i,:)^3-1/2*xP(i,:)^2*(xP(j,:)+xP(j-1,:))/2);
    end
    sum_g_M(i,:)=sum_g_M1;
    sum_g_theta(i,:)=sum_g_theta1;
    sum_g_omega(i,:)=sum_g_omega1;

    M(i,:)=-Hx*(Coef_Eta(i,1)+Coef_Eta(i,2)*xP(i)+Coef_Eta(i,3)*xP(i)^2/2+Coef_Eta(i,4)*xP(i)^3/6+Coef_Eta(i,5)*xP(i)^4/24)-sum_g_M(i,:)+sum_PM;    % M(i,:)=-Hx*ZG(i,:)-sum_gM(i,:)+sum_PM;
    theta(i,:)=-1/Eg/Ig*(-Hx*(Coef_Eta(i,1)*xP(i,:)+Coef_Eta(i,2)*xP(i,:)^2/2+Coef_Eta(i,3)*xP(i,:)^3/6+Coef_Eta(i,4)*xP(i,:)^4/24+Coef_Eta(i,5)*xP(i,:)^5/120)+sum_Ptheta-sum_g_theta(i,:)+C(i,:));
    omega(i,:)=-1/Eg/Ig*(-Hx*(Coef_Eta(i,1)*xP(i,:)^2/2+Coef_Eta(i,2)*xP(i,:)^3/6+Coef_Eta(i,3)*xP(i,:)^4/24+Coef_Eta(i,4)*xP(i,:)^5/120+Coef_Eta(i,5)*xP(i,:)^6/720)+sum_Pomega-sum_g_omega(i,:)+xP(i,:)*C(i,:)+D(i,:));
end

% 根据ZH求主梁0.1m的弯矩
m=2;
j=1;
for i=0:0.1:xP(end,:)
    j=j+1;
    if i<=xP(m)
        C_1m(j,:)=C(m);
        D_1m(j,:)=D(m);
    else
        m=m+1;
        C_1m(j,:)=C(m);
        D_1m(j,:)=D(m);
    end
end
% 根据ZH求主梁0.1m的弯矩
m=2;
j=1;
for i=0:0.1:xP(end,:) 
    sum_PM=0;
    sum_g_M1=0;
    sum_g_M2=0;
    if i<=xP(m)
        for k=1:m-1
            sum_PM=sum_PM+Pz(k,:)*(i-xP(k,:));
        end
        if m>=3
            for k=1:m-2
                sum_g_M1=sum_g_M1+qg(k,:)*(xP(k+1,:)-xP(k,:))*(i-(xP(k+1,:)+xP(k,:))/2);
            end
            sum_g_M2=qg(m-1,:)*(i-xP(m-1,:))*(i-(i+xP(m-1,:))/2);
        else 
            sum_g_M1=0;
            sum_g_M2=qg(m-1,:)*(i-xP(m-1,:))*(i-(i+xP(m-1,:))/2);
        end
    else 
        m=m+1;
        for k=1:m-1
            sum_PM=sum_PM+Pz(k,:)*(i-xP(k,:));
        end
        if m>=3
            for k=1:m-2
                sum_g_M1=sum_g_M1+qg(k,:)*(xP(k+1,:)-xP(k,:))*(i-(xP(k+1,:)+xP(k,:))/2);
            end
            sum_g_M2=qg(m-1,:)*(i-xP(m-1,:))*(i-(i+xP(m-1,:))/2);
        else 
            sum_g_M1=0;
            sum_g_M2=qg(m-1,:)*(i-xP(m-1,:))*(i-(i+xP(m-1,:))/2);
        end
    end
    M_01_Hx(j,:)=-Hx*(Coef_Eta(1,1)+Coef_Eta(1,2)*i+Coef_Eta(1,3)*i^2/2+Coef_Eta(1,4)*i^3/6+Coef_Eta(1,5)*i^4/24);
    M_01_g(j,:)=-sum_g_M1-sum_g_M2;
    M_01_P(j,:)=sum_PM;
    M_01(j,:)=M_01_Hx(j,:)+M_01_g(j,:)+M_01_P(j,:);
    
    j=j+1;
end
M_01=M_01/1000;

% 主梁坐标
XG=xP;
YG=Hg/2*ones(nL_hanger+nM_hanger+nR_hanger+n_zhizuo,1);
ZG=HA+ZG-omega;     % 梁的z与主缆的z相反，ZG=ZG+omega;


% 第二步 控制方程
% 高差 3个方程 【7个】
F(1,:)=(ZO(nL_hanger+2)-HB);          % 左跨高程闭合  sum_deltzL应该<0   F(1,:)=(sum_deltzL+abs((HA-HB))); 
F(2,:)=(YO(nL_hanger+2)-YB);          % 左跨横桥向闭合  sum_deltyL应该<0
F(3,:)=(ZO(nL_hanger+nM_hanger+3)-HC);                % 中跨高程闭合  sum_deltzM应该=0
F(4,:)=(ZO(nL_hanger+ceil(nM_hanger/2)+2)-H_Mid);        % 跨中中点高程闭合  sum_deltzM_middle应该>0
F(5,:)=(YO(nL_hanger+2)-YO(nL_hanger+nM_hanger+3));                       % 中跨横桥向闭合   sum_deltyM应该=0
F(6,:)=(ZO(nL_hanger+nM_hanger+nR_hanger+n_zhizuo)-HD);          % 右跨高程闭合  sum_deltzR应该>0
F(7,:)=(YO(nL_hanger+nM_hanger+nR_hanger+n_zhizuo)-YD);          % 左跨横桥向闭合 sum_deltyR应该>0

% 吊杆高程守恒 【nL_hanger+nM_hanger+nR_hanger个】
for i=1:nL_hanger
F(i+7,:)=(Y_hL(nz_hL,i)-YG(i+1,:));                % 【吊杆高程】 F中【8~14】
end
for i=1:nM_hanger
F(i+7+nL_hanger,:)=(Y_hM(nz_hM,i)-YG(i+nL_hanger+2,:));                % 【吊杆高程】 F中【15~41】
end
for i=1:nR_hanger
F(i+7+nL_hanger+nM_hanger,:)=(Y_hR(nz_hR,i)-YG(i+nL_hanger+nM_hanger+3,:));      % 【吊杆高程】 F中【42~48】
end

% 各吊点挠度=0  【nL_hanger+nM_hanger+nR_hanger+n_zhizuo-1个】 (吊杆数+支座数-2,左端点挠度为0已用于求D1,右端点用于求C1) 
for i=1:nL_hanger
F(i+7+nL_hanger+nM_hanger+nR_hanger,:)=(Z_hL(nz_hL,i)-ZG(i+1,:));                % 【吊杆高程】 F中【49~55】
end
for i=1:nM_hanger
F(i+7+2*nL_hanger+nM_hanger+nR_hanger,:)=(Z_hM(nz_hM,i)-ZG(i+nL_hanger+2,:));                % 【吊杆高程】 F中【56~82】
end
for i=1:nR_hanger
F(i+7+2*nL_hanger+2*nM_hanger+nR_hanger,:)=(Z_hR(nz_hR,i)-ZG(i+nL_hanger+nM_hanger+3,:));      % 【吊杆高程】 F中【83~89】
end

% 桥塔处主梁w=0
F(7+2*nL_hanger+2*nM_hanger+2*nR_hanger+1,:)=omega(nL_hanger+2,:);     % 【吊杆高程】 F中【90】
F(7+2*nL_hanger+2*nM_hanger+2*nR_hanger+2,:)=omega(nL_hanger+nM_hanger+3,:);     % 【吊杆高程】 F中【91】

% 主梁
F(6+2*nL_hanger+2*nM_hanger+2*nR_hanger+n_zhizuo,:)=(sum(Pz)-Gg_all);       % 主梁竖向力平衡 F中【92】
%F(7+2*nL_hanger+2*nM_hanger+2*nR_hanger+n_zhizuo,:)=(HL(1,:)*cos(alpha_L(1,:))-Hx/2);  % 左跨主缆水平分力与Hx相等  F中【93】
F(7+2*nL_hanger+2*nM_hanger+2*nR_hanger+n_zhizuo,:)=(HR(nR_hanger+1,:)*cos(alpha_R(nR_hanger+1,:))-Hx/2);  % 右跨主缆水平分力与Hx相等  F中【93】
F(8+2*nL_hanger+2*nM_hanger+2*nR_hanger+n_zhizuo,:)=M(end,:);    % F中【94】

% 桥塔
F(9+2*nL_hanger+2*nM_hanger+2*nR_hanger+n_zhizuo,:)=(HL(nL_hanger+1,:)*cos(alpha_L(nL_hanger+1,:))-HM(1,:)*cos(alpha_M(1,:))); % F中【95】
F(10+2*nL_hanger+2*nM_hanger+2*nR_hanger+n_zhizuo,:)=(HR(1,:)*cos(alpha_R(1,:))-HM(nM_hanger+1,:)*cos(alpha_M(nM_hanger+1,:))); % F中【96】
% 未知数个数： 10+2*nL_hanger+2*nM_hanger+2*nR_hanger+n_zhizuo
% 控制方程个数：10+2*nL_hanger+2*nM_hanger+2*nR_hanger+n_zhizuo

% 反向求左端点弯矩
sum_PM_offset=0;
sum_g_M_offset=0;
for j=1:nL_hanger+nM_hanger+nR_hanger+n_zhizuo
    sum_PM_offset=sum_PM_offset+Pz(j,:)*xP(j,:);
end
for j=2:nL_hanger+nM_hanger+nR_hanger+n_zhizuo
    sum_g_M_offset=sum_g_M_offset+qg(j-1,:)*(xP(j,:)-xP(j-1,:))*(xP(j,:)+xP(j-1,:))/2;
end
M0=-Hx*(ZG(end,:)-ZG(1,:))+sum_PM_offset-sum_g_M_offset;
F(11+2*nL_hanger+2*nM_hanger+2*nR_hanger+n_zhizuo,:)=M0;


% 主缆初应变
for i=1:nL_hanger+1
    epsilon_cable(i,:)=(LL(i,:)-SL(i,:))/LL(i,:);
    Phi_cL(i,:)=atan((ZO(i+1)-ZO(i))/(XO(i+1)-XO(i)));
end
for i=1:nM_hanger+1
    epsilon_cable(i+nL_hanger+1,:)=(LM(i,:)-SM(i,:))/LM(i,:);
    Phi_cM(i,:)=atan((ZO(i+nL_hanger+2)-ZO(i+nL_hanger+1))/(XO(i+nL_hanger+2)-XO(i+nL_hanger+1)));
end
for i=1:nR_hanger+1
    epsilon_cable(i+nL_hanger+nM_hanger+2,:)=(LR(i,:)-SR(i,:))/LR(i,:);
    Phi_cR(i,:)=atan((ZO(i+nL_hanger+nM_hanger+3)-ZO(i+nL_hanger+nM_hanger+2))/(XO(i+nL_hanger+nM_hanger+3)-XO(i+nL_hanger+nM_hanger+2)));
end
% 支反力
R_support(1,:)=Pz(1,:)-Hx*tan(Phi_cL(1,:));
R_support(2,:)=Pz(nL_hanger+2,:);
R_support(3,:)=Pz(nL_hanger+nM_hanger+3,:);
R_support(4,:)=Pz(end,:)-Hx*tan(-Phi_cR(end,:));

% 主梁的初应变
% 主梁Z的计算
size_xP_1m=xP(end,:);
xP_1m(1,:)=0;
for i=1:size_xP_1m
    xP_1m(i+1,:)=xP_1m(i,:)+1;
end
XG_1m=xP_1m;
YG_1m=zeros(size_xP_1m+1,1);
for i=1:size_xP_1m+1
    ZG_1m(i,:)=Coef_Eta(1,1)+Coef_Eta(1,2)*xP_1m(i)+Coef_Eta(1,3)*xP_1m(i)^2/2+Coef_Eta(1,4)*xP_1m(i)^3/6+Coef_Eta(1,5)*xP_1m(i)^4/24;
end
for i=1:size_xP_1m
    k_1m(i,:)=(ZG_1m(i+1,:)-ZG_1m(i,:))/(xP_1m(i+1,:)-xP_1m(i,:));
    epsilon_beam(i,:)=Hx*sqrt(1+k_1m(i,:)^2)/Eg/Ag;
end
ZG_1m=HA+ZG_1m;

% 吊杆的初应变+建模坐标
mergedData_XYZ_h = [];
mergedData_epsilon_h = [];
% 左跨
mergedData_XYZ_hL = [];
mergedData_epsilon_hL = [];
for i=1:nL_hanger
    for j=1:nz_hL
        mergedData_XYZ_hL = [mergedData_XYZ_hL; X_hL(j,i), Y_hL(j,i), Z_hL(j,i)];
        mergedData_epsilon_hL = [mergedData_epsilon_hL; epsilon_hanger_L(j,i)];
    end
    PL_cable(i,:)=sqrt(PLy(i,:)^2+PLz_dot(i,:)^2);    
    PL_beam(i,:)=sqrt(PLy(i,:)^2+PLz(i,:)^2);
end
% 中跨
mergedData_XYZ_hM = [];
mergedData_epsilon_hM = [];
for i=1:nM_hanger
    for j=1:nz_hM
        mergedData_XYZ_hM = [mergedData_XYZ_hM; X_hM(j,i), Y_hM(j,i), Z_hM(j,i)];
        mergedData_epsilon_hM = [mergedData_epsilon_hM; epsilon_hanger_M(j,i)];
    end
    PM_cable(i,:)=sqrt(PMy(i,:)^2+PMz_dot(i,:)^2);    
    PM_beam(i,:)=sqrt(PMy(i,:)^2+PMz(i,:)^2);
end
% 右跨
mergedData_XYZ_hR = [];
mergedData_epsilon_hR = [];
for i=1:nR_hanger
    for j=1:nz_hR
        mergedData_XYZ_hR = [mergedData_XYZ_hR; X_hR(j,i), Y_hR(j,i), Z_hR(j,i)];
        mergedData_epsilon_hR = [mergedData_epsilon_hR; epsilon_hanger_R(j,i)];
    end
    PR_cable(i,:)=sqrt(PRy(i,:)^2+PRz_dot(i,:)^2);    
    PR_beam(i,:)=sqrt(PRy(i,:)^2+PRz(i,:)^2);
end
% 合并
mergedData_XYZ_h = [mergedData_XYZ_hL;mergedData_XYZ_hM;mergedData_XYZ_hR];
mergedData_epsilon_h = [mergedData_epsilon_hL;mergedData_epsilon_hM;mergedData_epsilon_hR];

% 其它数据
P_cable=[PL_cable;PM_cable;PR_cable];
P_beam=[PL_beam;PM_beam;PR_beam];

P_cable_z=[PLz_dot;PMz_dot;PRz_dot];
P_beam_z=[PLz;PMz;PRz];
P_delt_z=P_cable_z-P_beam_z;

G_Hanger=[gamma_hL;gamma_hM;gamma_hR];
M(1,:)=M0;
%% 导入到ansys建模文件夹中

% 导出吊杆建模坐标-上游
subfolderPath = fullfile(pwd, 'Ansys_Modeling_coordinates_01m'); 
fullFileName_F = fullfile(subfolderPath, 'Cor_hanger_F.txt'); 
fileID_F = fopen(fullFileName_F, 'w');
% 写入列标题  
fprintf(fileID_F, '%s\t%s\t%s\t%s\n', '0', '1', '2', '3');  
% 遍历向量并写入数据  
for k=1:length(mergedData_XYZ_h)
    fprintf(fileID_F, '%d\t%.12f\t%.12f\t%.12f\n', k, mergedData_XYZ_h(k,1), mergedData_XYZ_h(k,3), mergedData_XYZ_h(k,2));  
end
% 关闭文件  
fclose(fileID_F);

% 导出吊杆建模坐标-下游
fullFileName_B = fullfile(subfolderPath, 'Cor_hanger_B.txt'); 
fileID_B = fopen(fullFileName_B, 'w');
% 写入列标题  
fprintf(fileID_B, '%s\t%s\t%s\t%s\n', '0', '1', '2', '3');  
% 遍历向量并写入数据  
for k=1:length(mergedData_XYZ_h)
    fprintf(fileID_B, '%d\t%.12f\t%.12f\t%.12f\n', k, mergedData_XYZ_h(k,1), mergedData_XYZ_h(k,3), -mergedData_XYZ_h(k,2)); 
end
% 关闭文件  
fclose(fileID_B);

% 导出吊杆应变
fullFileName = fullfile(subfolderPath, 'Real_Steel.txt'); 
fileID = fopen(fullFileName, 'w');
% 写入列标题  
fprintf(fileID, '%s\t%s\t%s\n', '0', '1', '2'); 
% 遍历向量并写入数据   
for i = 1:length(mergedData_epsilon_h)  
    fprintf(fileID, '%d\t%.12f\t%.12f\n', i, Ah, mergedData_epsilon_h(i));  
end  
% 关闭文件  
fclose(fileID);


% 主缆建模坐标
fullFileName= fullfile(subfolderPath, 'Cor_ZhuLan_B.txt');
fileID = fopen(fullFileName, 'w'); 
% 写入列标题  
fprintf(fileID, '%s\t%s\t%s\t%s\n', '0', '1', '2', '3');  
% 遍历向量并写入数据  
for i = 2:length(XO)-1  
    fprintf(fileID, '%d\t%.12f\t%.12f\t%.12f\n', i-1, XO(i), ZO(i), -YO(i));  
end  
% 关闭文件  
fclose(fileID);

fullFileName= fullfile(subfolderPath, 'Cor_ZhuLan_F.txt');
fileID = fopen(fullFileName, 'w'); 
% 写入列标题  
fprintf(fileID, '%s\t%s\t%s\t%s\n', '0', '1', '2', '3');  
% 遍历向量并写入数据  
for i = 2:length(XO)-1  
    fprintf(fileID, '%d\t%.12f\t%.12f\t%.12f\n', i-1, XO(i), ZO(i), YO(i));  
end  
% 关闭文件  
fclose(fileID);


% 主缆应变
fullFileName= fullfile(subfolderPath, 'Real_ZhuLan.txt');
fileID = fopen(fullFileName, 'w'); 
% 写入列标题  
fprintf(fileID, '%s\t%s\t%s\n', '0', '1', '2');  
% 遍历向量并写入数据  
for i = 1:length(epsilon_cable)  
    fprintf(fileID, '%d\t%.12f\t%.12f\n', i, Ac, epsilon_cable(i));  
end  
% 关闭文件  
fclose(fileID);

% 主梁建模坐标
fullFileName= fullfile(subfolderPath, 'Cor_Grid_beam_01m.txt');
fileID = fopen(fullFileName, 'w'); 
% 写入列标题  
fprintf(fileID, '%s\t%s\t%s\t%s\n', '0', '1', '2', '3');  
% 遍历向量并写入数据  
for i = 1:length(ZG_1m)  
    fprintf(fileID, '%d\t%.12f\t%.12f\t%.12f\n', i, XG_1m(i), ZG_1m(i), 0);  
end  
% 关闭文件  
fclose(fileID);

% 鱼骨梁建模坐标 
fullFileName= fullfile(subfolderPath, 'yugu_shang.txt');
fileID = fopen(fullFileName, 'w'); 
% 写入列标题  
fprintf(fileID, '%s\t%s\t%s\t%s\n', '0', '1', '2', '3');  
% 遍历向量并写入数据  
for i = 1:length(ZG)  
    fprintf(fileID, '%d\t%.12f\t%.12f\t%.12f\n', i, XG(i), ZG(i), YG(i));  
end  
% 关闭文件  
fclose(fileID);
 
fullFileName= fullfile(subfolderPath, 'yugu_xia.txt');
fileID = fopen(fullFileName, 'w'); 
% 写入列标题  
fprintf(fileID, '%s\t%s\t%s\t%s\n', '0', '1', '2', '3');  
% 遍历向量并写入数据  
for i = 1:length(ZG)  
    fprintf(fileID, '%d\t%.12f\t%.12f\t%.12f\n', i, XG(i), ZG(i), -YG(i));  
end  
% 关闭文件  
fclose(fileID);

% 主梁应变
fullFileName= fullfile(subfolderPath, 'Real_Beam_01m.txt');
fileID = fopen(fullFileName, 'w'); 
% 写入列标题  
fprintf(fileID, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', '0', '1', '2', '3', '4', '5', '6', '7', '8');  
% 遍历向量并写入数据  
for i = 1:length(epsilon_beam)  
    fprintf(fileID, '%d\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\n', i, Ag_JM, Ig_zz, Ig_yy, h_zz, h_yy, theta_x, -epsilon_beam(i), Ig_xx);     % epsilon_beam(i)
end  
% 关闭文件  
fclose(fileID);

