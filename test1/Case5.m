warning off;
close all;
clear all;
clc;

% false  true
options = optimoptions('fsolve','InitDamping',1e-100,...
                       'ScaleProblem','jacobian',...
                       'Diagnostics','on',...`
                       'UseParallel',false,...
                       'FiniteDifferenceStepSize', eps^(1/3),...
                       'FiniteDifferenceType','central',...
                       'TolFun',1.0e-10,...
                       'TolX',1.0e-10,...
                       'MaxFunEvals',8e30,... 
                       'MaxIter',200000,...
                       'Display','iter',...
                       'Algorithm','Levenberg-Marquardt');

% Import initial value
load 'x_Case4.txt';
x=x_Case4;
chuzhi_3D_hangerXLX=x;


% Call the main function
[x,fval,exitflag,output]=fsolve(@(YY)CQ_self_anchored(YY),chuzhi_3D_hangerXLX,options);
    
if exitflag > 0    
    message_KL = sprintf('Solved, solver stopped');  
else  
    message_KL = sprintf('No solution found');  
end




%% Main Function
function F = CQ_self_anchored(x)

%------------------全桥参数（永宗大桥）------------------
% Cable
Ec=2.05e8;      Ac=0.0715;     qc=6.2;
HA=29.0197;    HB=97.7050;       HC=99.9150;   HD=29.0197;     H_Mid=41.0320;
YA=22.4409;    YB=1;       YC=1;   YD=22.4409;

% Hanger
Eh=2.06e8;   Ah=0.00609;   qh=0.468981765;
nL_hanger=7;  nM_hanger=27;  nR_hanger=7;

% Stiffening girder
Eg=2.06e8;     Ig=4;    Ag=1.12;
Hg=22.4409*2;  n_zhizuo=4;
qg=200*ones(nL_hanger+nM_hanger+nR_hanger+n_zhizuo,1);

% Import known quantities
load 'xP.txt';            % The distance from the hanger point to the left end of the stiffening girder
load 'Coef_Eta4.txt';     % polynomial coefficient of the stiffening girder
Coef_Eta=Coef_Eta4;

% Basic unknown parameters
aM(1,:)=x(1); HM(1,:)=x(2); alpha_M(1,:)=x(3);
aL(1,:)=x(4); HL(1,:)=x(5); alpha_L(1,:)=x(6);
aR(1,:)=x(7); HR(1,:)=x(8); alpha_R(1,:)=x(9);
Hx=2*x(10);
for i=1:nL_hanger+nM_hanger+nR_hanger+n_zhizuo
Py(i,:)=x(10+i); 
end
for i=1:nL_hanger 
ahL(i,:)=x(10+nL_hanger+nM_hanger+nR_hanger+n_zhizuo+i);
end
for i=1:nM_hanger 
ahM(i,:)=x(10+2*nL_hanger+nM_hanger+nR_hanger+n_zhizuo+i);
end
for i=1:nR_hanger 
ahR(i,:)=x(10+2*nL_hanger+2*nM_hanger+nR_hanger+n_zhizuo+i);
end

% Z-coordinate calculation of the stiffening girder
for i=1:nL_hanger+nM_hanger+nR_hanger+n_zhizuo
    ZG(i,:)=Coef_Eta(i,1)+Coef_Eta(i,2)*xP(i)+Coef_Eta(i,3)*xP(i)^2/2+Coef_Eta(i,4)*xP(i)^3/6+Coef_Eta(i,5)*xP(i)^4/24;
end
YGL=(YA-YB)*ones(nL_hanger,1); 
YGM=(YA-YB)*ones(nM_hanger,1);
YGR=(YA-YB)*ones(nR_hanger,1);

% Horizontal projection length of main cable suspension line segment
for i=1:nL_hanger+nM_hanger+nR_hanger+n_zhizuo-1
    l(i,1)=xP(i+1)-xP(i);
end
lL=l(1:nL_hanger+1,1);
lM=l(nL_hanger+2:nL_hanger+nM_hanger+2,1);
lR=l(nL_hanger+nM_hanger+3:nL_hanger+nM_hanger+nR_hanger+3,1);

for i=1:nL_hanger
    PLy(i,:)=Py(i+1,:);
    ZGL(i,:)=ZG(i+1,:);
end
for i=1:nM_hanger
    PMy(i,:)=Py(i+nL_hanger+2,:);
    ZGM(i,:)=ZG(i+nL_hanger+2,:);
end
for i=1:nR_hanger
    PRy(i,:)=Py(i+nL_hanger+nM_hanger+3,:);
    ZGR(i,:)=ZG(i+nL_hanger+nM_hanger+3,:);
end

%% Cable
% Left side span
for i=1:nL_hanger+1
    cL(i,:)=-HL(i,:)/qc; 
    deltyL(i,:)=lL(i,:)*tan(alpha_L(i,:)); 
    deltzL(i,:)=cL(i,:)*(cosh(lL(i,:)/cL(i,:)/cos(alpha_L(i,:))+aL(i,:))-cosh(aL(i,:)));
    LL(i,:)=cL(i,:)*(sinh(lL(i,:)/cL(i,:)/cos(alpha_L(i,:))+aL(i,:))-sinh(aL(i,:)));
    SL(i,:)=cL(i,:)*(sinh(lL(i,:)/cL(i,:)/cos(alpha_L(i,:))+aL(i,:))-sinh(aL(i,:)))-HL(i,:)/2/Ec/Ac*(lL(i,:)/cos(alpha_L(i,:))+cL(i,:)/2*(sinh(2*(lL(i,:)/cL(i,:)/cos(alpha_L(i,:))+aL(i,:)))-sinh(2*aL(i,:))));
    j=i;
    if j<=nL_hanger
        chL(j,:)=-PLy(j,:)/qh;
        lhL(j,:)=-sum(deltyL(1:j,:));  
        deltzhL(j,:)=chL(j,:)*(cosh(lhL(j,:)/chL(j,:)+ahL(j,:))-cosh(ahL(j,:))); 
        PLz(j,:)=PLy(j,:)*sinh(lhL(j,:)/chL(j,:)+ahL(j,:)); 
        LhL(j,:)=chL(j,:)*(sinh(lhL(j,:)/chL(j,:)+ahL(j,:))-sinh(ahL(j,:)));
        ShL(j,:)=chL(j,:)*(sinh(lhL(j,:)/chL(j,:)+ahL(j,:))-sinh(ahL(j,:)))-PLy(j,:)/2/Eh/Ah*(lhL(j,:)+chL(j,:)/2*(sinh(2*(lhL(j,:)/chL(j,:)+ahL(j,:)))-sinh(2*ahL(j,:))));
        gamma_hL(j,:)=qh*ShL(j,:); 
        PLz_dot2(j,:)=PLz(j,:)+gamma_hL(j,:)+gamma_c;
        PLz2(j,:)=PLz_dot2(j,:)-gamma_hL(j,:)-gamma_c;
        PLz_dot(j,:)=PLy(j,:)*sinh(ahL(j,:)); 
        % Parametric recursive formula 
        alpha_L(j+1,:)=atan((HL(j,:)*sin(alpha_L(j,:))-PLy(j,:))/(HL(j,:)*cos(alpha_L(j,:)))); 
        HL(j+1,:)=HL(j,:)*cos(alpha_L(j,:))*sec(alpha_L(j+1,:)); 
        aL(j+1,:)=asinh((HL(j,:)*sinh(lL(j,:)/cL(j,:)/cos(alpha_L(j,:))+aL(j,:))-PLz_dot(j,:))/HL(j+1,:)); 
    else
        break
    end
end

% Main span
for i=1:nM_hanger+1
    cM(i,:)=-HM(i,:)/qc;
    deltyM(i,:)=lM(i,:)*tan(alpha_M(i,:));
    deltzM(i,:)=cM(i,:)*(cosh(lM(i,:)/cM(i,:)/cos(alpha_M(i,:))+aM(i,:))-cosh(aM(i,:)));
    LM(i,:)=cM(i,:)*(sinh(lM(i,:)/cM(i,:)/cos(alpha_M(i,:))+aM(i,:))-sinh(aM(i,:)));
    SM(i,:)=cM(i,:)*(sinh(lM(i,:)/cM(i,:)/cos(alpha_M(i,:))+aM(i,:))-sinh(aM(i,:)))-HM(i,:)/2/Ec/Ac*(lM(i,:)/cos(alpha_M(i,:))+cM(i,:)/2*(sinh(2*(lM(i,:)/cM(i,:)/cos(alpha_M(i,:))+aM(i,:)))-sinh(2*aM(i,:))));
    j=i;
    if j<=nM_hanger 
        chM(j,:)=-PMy(j,:)/qh; 
        lhM(j,:)=YGM(j,:)-sum(deltyM(1:j,:));
        deltzhM(j,:)=chM(j,:)*(cosh(lhM(j,:)/chM(j,:)+ahM(j,:))-cosh(ahM(j,:)));
        PMz(j,:)=PMy(j,:)*sinh(lhM(j,:)/chM(j,:)+ahM(j,:));
        LhM(j,:)=chM(j,:)*(sinh(lhM(j,:)/chM(j,:)+ahM(j,:))-sinh(ahM(j,:)));
        ShM(j,:)=chM(j,:)*(sinh(lhM(j,:)/chM(j,:)+ahM(j,:))-sinh(ahM(j,:)))-PMy(j,:)/2/Eh/Ah*(lhM(j,:)+chM(j,:)/2*(sinh(2*(lhM(j,:)/chM(j,:)+ahM(j,:)))-sinh(2*ahM(j,:))));
        gamma_hM(j,:)=qh*ShM(j,:); 
        PMz_dot2(j,:)=PMz(j,:)+gamma_hM(j,:)+gamma_c;
        PMz2(j,:)=PMz_dot2(j,:)-gamma_hM(j,:)-gamma_c;
        PMz_dot(j,:)=PMy(j,:)*sinh(ahM(j,:));
        % Parametric recursive formula
        alpha_M(j+1,:)=atan((HM(j,:)*sin(alpha_M(j,:))-PMy(j,:))/(HM(j,:)*cos(alpha_M(j,:))));
        HM(j+1,:)=HM(j,:)*cos(alpha_M(j,:))*sec(alpha_M(j+1,:));
        aM(j+1,:)=asinh((HM(j,:)*sinh(lM(j,:)/cM(j,:)/cos(alpha_M(j,:))+aM(j,:))-PMz_dot(j,:))/HM(j+1,:));
    else
        break
    end
end

% Right side span
for i=1:nR_hanger+1
    cR(i,:)=-HR(i,:)/qc;
    deltyR(i,:)=lR(i,:)*tan(alpha_R(i,:));
    deltzR(i,:)=cR(i,:)*(cosh(lR(i,:)/cR(i,:)/cos(alpha_R(i,:))+aR(i,:))-cosh(aR(i,:)));
    LR(i,:)=cR(i,:)*(sinh(lR(i,:)/cR(i,:)/cos(alpha_R(i,:))+aR(i,:))-sinh(aR(i,:)));
    SR(i,:)=cR(i,:)*(sinh(lR(i,:)/cR(i,:)/cos(alpha_R(i,:))+aR(i,:))-sinh(aR(i,:)))-HR(i,:)/2/Ec/Ac*(lR(i,:)/cos(alpha_R(i,:))+cR(i,:)/2*(sinh(2*(lR(i,:)/cR(i,:)/cos(alpha_R(i,:))+aR(i,:)))-sinh(2*aR(i,:))));
    j=i;
    if j<=nR_hanger
        chR(j,:)=-PRy(j,:)/qh;
        lhR(j,:)=YGR(j,:)-sum(deltyR(1:j,:));
        deltzhR(j,:)=chR(j,:)*(cosh(lhR(j,:)/chR(j,:)+ahR(j,:))-cosh(ahR(j,:)));
        PRz(j,:)=PRy(j,:)*sinh(lhR(j,:)/chR(j,:)+ahR(j,:));
        LhR(j,:)=chR(j,:)*(sinh(lhR(j,:)/chR(j,:)+ahR(j,:))-sinh(ahR(j,:)));
        ShR(j,:)=chR(j,:)*(sinh(lhR(j,:)/chR(j,:)+ahR(j,:))-sinh(ahR(j,:)))-PRy(j,:)/2/Eh/Ah*(lhR(j,:)+chR(j,:)/2*(sinh(2*(lhR(j,:)/chR(j,:)+ahR(j,:)))-sinh(2*ahR(j,:))));
        gamma_hR(j,:)=qh*ShR(j,:);
        PRz_dot2(j,:)=PRz(j,:)+gamma_hR(j,:)+gamma_c;
        PRz2(j,:)=PRz_dot2(j,:)-gamma_hR(j,:)-gamma_c;
        PRz_dot(j,:)=PRy(j,:)*sinh(ahR(j,:));
        % Parametric recursive formula
        alpha_R(j+1,:)=atan((HR(j,:)*sin(alpha_R(j,:))-PRy(j,:))/(HR(j,:)*cos(alpha_R(j,:))));
        HR(j+1,:)=HR(j,:)*cos(alpha_R(j,:))*sec(alpha_R(j+1,:)); 
        aR(j+1,:)=asinh((HR(j,:)*sinh(lR(j,:)/cR(j,:)/cos(alpha_R(j,:))+aR(j,:))-PRz_dot(j,:))/HR(j+1,:)); 
    else
        break
    end
end
sum_deltzR=sum(deltzR);  % zC-zD
sum_deltyR=sum(deltyR);  % yC-yD

% Coordinates of the cable 
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


%% Hanger - Catenary model
% Left side span
nz_hL=10;
for i=1:nL_hanger    
    for j=1:nz_hL
        lhL_1(i,j)=j*lhL(i,:)/nz_hL;
        zhL(i,j)=chL(i,:)*(cosh(lhL_1(i,j)/chL(i,:)+ahL(i,:))-cosh(ahL(i,:)));
        LhL_1(i,j)=chL(i,:)*(sinh(lhL_1(i,j)/chL(i,:)+ahL(i,:))-sinh(ahL(i,:)));
        ShL_1(i,j)=chL(i,:)*(sinh(lhL_1(i,j)/chL(i,:)+ahL(i,:))-sinh(ahL(i,:)))-PLy(i,:)/2/Eh/Ah*(lhL_1(i,j)+chL(i,:)/2*(sinh(2*(lhL_1(i,j)/chL(i,:)+ahL(i,:)))-sinh(2*ahL(i,:))));
        LhL_2(i,1)=LhL_1(i,1);
        ShL_2(i,1)=ShL_1(i,1);
        if j>=2
            LhL_2(i,j)=LhL_1(i,j)-LhL_1(i,j-1);
            ShL_2(i,j)=ShL_1(i,j)-ShL_1(i,j-1);
        end
        epsilon_hanger_L(j,i)=(LhL_2(i,j)-ShL_2(i,j))/LhL_2(i,j);
        X_hL(j,i)=XO(i+1,:);
        Z_hL(j,i)=ZO(i+1,:)-zhL(i,j);
        Y_hL(j,i)=YO(i+1,:)+lhL_1(i,j);
    end
end
% Main span
nz_hM=10;
for i=1:nM_hanger    
    for j=1:nz_hM
        lhM_1(i,j)=j*lhM(i,:)/nz_hM;
        zhM(i,j)=chM(i,:)*(cosh(lhM_1(i,j)/chM(i,:)+ahM(i,:))-cosh(ahM(i,:)));
        LhM_1(i,j)=chM(i,:)*(sinh(lhM_1(i,j)/chM(i,:)+ahM(i,:))-sinh(ahM(i,:)));
        ShM_1(i,j)=chM(i,:)*(sinh(lhM_1(i,j)/chM(i,:)+ahM(i,:))-sinh(ahM(i,:)))-PMy(i,:)/2/Eh/Ah*(lhM_1(i,j)+chM(i,:)/2*(sinh(2*(lhM_1(i,j)/chM(i,:)+ahM(i,:)))-sinh(2*ahM(i,:))));
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
% Right side span
nz_hR=10;
for i=1:nR_hanger    
    for j=1:nz_hR
        lhR_1(i,j)=j*lhR(i,:)/nz_hR;
        zhR(i,j)=chR(i,:)*(cosh(lhR_1(i,j)/chR(i,:)+ahR(i,:))-cosh(ahR(i,:)));
        LhR_1(i,j)=chR(i,:)*(sinh(lhR_1(i,j)/chR(i,:)+ahR(i,:))-sinh(ahR(i,:)));
        ShR_1(i,j)=chR(i,:)*(sinh(lhR_1(i,j)/chR(i,:)+ahR(i,:))-sinh(ahR(i,:)))-PRy(i,:)/2/Eh/Ah*(lhR_1(i,j)+chR(i,:)/2*(sinh(2*(lhR_1(i,j)/chR(i,:)+ahR(i,:)))-sinh(2*ahR(i,:)))); 
        LhR_2(i,1)=LhR_1(i,1);
        ShR_2(i,1)=ShR_1(i,1);
        if j>=2
            LhR_2(i,j)=LhR_1(i,j)-LhR_1(i,j-1);
            ShR_2(i,j)=ShR_1(i,j)-ShR_1(i,j-1);
        end
        epsilon_hanger_R(j,i)=(LhR_2(i,j)-ShR_2(i,j))/LhR_2(i,j);
        X_hR(j,i)=XO(i+nL_hanger+nM_hanger+3,:);
        Z_hR(j,i)=ZO(i+nL_hanger+nM_hanger+3,:)-zhR(i,j);
        Y_hR(j,i)=YO(i+nL_hanger+nM_hanger+3,:)+lhR_1(i,j);
    end
end


%% Stiffening girder
for i=1:nL_hanger+nM_hanger+nR_hanger+n_zhizuo-1
    y_dot = @(x) sqrt(1+(Coef_Eta(i,2)+Coef_Eta(i,3)*x).^2);  
    l_bend(i,:) = integral(y_dot, xP(i,:), xP(i+1,:));
    l_bend_sum(i+1,:) = sum(l_bend(1:i,:));
    Gg(i,:) = qg(i,:)*l_bend(i,:);
    qg(i,:) = Gg(i,:)/(xP(i+1,:)-xP(i,:));
end
l_bend_sum(1,:)=0;
Lg=sum(l_bend);
Gg_all=sum(Gg);

Pz(1,:)=2*Py(1,:);
Pz(nL_hanger+2,:)=2*Py(nL_hanger+2,:);
Pz(nL_hanger+nM_hanger+3,:)=2*Py(nL_hanger+nM_hanger+3,:);
Pz(nL_hanger+nM_hanger+nR_hanger+4,:)=2*Py(nL_hanger+nM_hanger+nR_hanger+4,:);
for i=1:nL_hanger
    Pz(i+1,:)=2*PLz(i,:); 
end
for i=1:nM_hanger
    Pz(i+nL_hanger+2,:)=2*PMz(i,:); 
end
for i=1:nR_hanger
    Pz(i+nL_hanger+nM_hanger+3,:)=2*PRz(i,:);
end

% Internal force of the stiffening girder
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
for i=1:nL_hanger+nM_hanger+nR_hanger+n_zhizuo-1
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

    M(i,:)=-Hx*(Coef_Eta(i,1)+Coef_Eta(i,2)*xP(i)+Coef_Eta(i,3)*xP(i)^2/2+Coef_Eta(i,4)*xP(i)^3/6+Coef_Eta(i,5)*xP(i)^4/24)-sum_g_M(i,:)+sum_PM;
    theta(i,:)=-1/Eg/Ig*(-Hx*(Coef_Eta(i,1)*xP(i,:)+Coef_Eta(i,2)*xP(i,:)^2/2+Coef_Eta(i,3)*xP(i,:)^3/6+Coef_Eta(i,4)*xP(i,:)^4/24+Coef_Eta(i,5)*xP(i,:)^5/120)+sum_Ptheta-sum_g_theta(i,:)+C(i,:));
    omega(i,:)=-1/Eg/Ig*(-Hx*(Coef_Eta(i,1)*xP(i,:)^2/2+Coef_Eta(i,2)*xP(i,:)^3/6+Coef_Eta(i,3)*xP(i,:)^4/24+Coef_Eta(i,4)*xP(i,:)^5/120+Coef_Eta(i,5)*xP(i,:)^6/720)+sum_Pomega-sum_g_omega(i,:)+xP(i,:)*C(i,:)+D(i,:));
end

XG=xP;
YG=Hg/2*ones(nL_hanger+nM_hanger+nR_hanger+n_zhizuo,1);
ZG=HA+ZG-omega;

for i=1:nL_hanger
    deltzhL2(i,:)=-sum(deltzL(1:i,:))-ZG(i+1,:)+HA;
end
for i=1:nM_hanger
    deltzhM2(i,:)=HB-sum(deltzM(1:i,:))-ZG(i+nL_hanger+2,:);
end
for i=1:nR_hanger
    deltzhR2(i,:)=HC-sum(deltzR(1:i,:))-ZG(i+nL_hanger+nM_hanger+3,:);
end

% Constraint equations
F(1,:)=(ZO(nL_hanger+2)-HB);
F(2,:)=(YO(nL_hanger+2)-YB);
F(3,:)=(ZO(nL_hanger+nM_hanger+3)-HC);
F(4,:)=(ZO(nL_hanger+ceil(nM_hanger/2)+2)-H_Mid);
F(5,:)=(YO(nL_hanger+2)-YO(nL_hanger+nM_hanger+3));
F(6,:)=(ZO(nL_hanger+nM_hanger+nR_hanger+n_zhizuo)-HD);
F(7,:)=(YO(nL_hanger+nM_hanger+nR_hanger+n_zhizuo)-YD);
for i=1:nL_hanger
F(i+7,:)=(deltzhL2(i,:)-deltzhL(i,:));
end
for i=1:nM_hanger
F(i+7+nL_hanger,:)=(deltzhM2(i,:)-deltzhM(i,:));
end
for i=1:nR_hanger
F(i+7+nL_hanger+nM_hanger,:)=(deltzhR2(i,:)-deltzhR(i,:));
end
for i=1:nL_hanger+nM_hanger+nR_hanger+n_zhizuo-2
    F(i+7+nL_hanger+nM_hanger+nR_hanger,:)=omega(i+1,:);
end
F(6+2*nL_hanger+2*nM_hanger+2*nR_hanger+n_zhizuo,:)=(sum(Pz)-Gg_all);
F(7+2*nL_hanger+2*nM_hanger+2*nR_hanger+n_zhizuo,:)=(HR(nR_hanger+1,:)*cos(alpha_R(nR_hanger+1,:))-Hx/2);
F(8+2*nL_hanger+2*nM_hanger+2*nR_hanger+n_zhizuo,:)=M(end,:);
F(9+2*nL_hanger+2*nM_hanger+2*nR_hanger+n_zhizuo,:)=(HL(nL_hanger+1,:)*cos(alpha_L(nL_hanger+1,:))-HM(1,:)*cos(alpha_M(1,:))); 
F(10+2*nL_hanger+2*nM_hanger+2*nR_hanger+n_zhizuo,:)=(HR(1,:)*cos(alpha_R(1,:))-HM(nM_hanger+1,:)*cos(alpha_M(nM_hanger+1,:))); 
end 
