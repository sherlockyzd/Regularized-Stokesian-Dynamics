clear;clc;close all;
load FarfieldScalar.mat
MX_xishu=load('../near-mid/MX.dat');
l=1;
lambda=MX_xishu(l,1);


% cita=Lb0_list(1:32);
% cita2=Lb0_list(17:27);
% cita3=Lb0_list(20:45);
cita=Lb0_list;
cita2=Lb0_list;
cita3=Lb0_list;
% X_M11ff=X11M(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% X_M12ff=X12M(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% X_M22ff=X11M(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% X_M21ff=X12M(Ls*(l-1)+20:Ls*(l-1)+45)./6;
X_M11ff=X11M(Ls*(l-1)+1:Ls*l)./6;
X_M12ff=X12M(Ls*(l-1)+1:Ls*l)./6;
X_M22ff=X11M(Ls*(l-1)+1:Ls*l)./6;
X_M21ff=X12M(Ls*(l-1)+1:Ls*l)./6;
x_mMX=5E1;
s=cita2+2;
sff=cita3+2;
%%
g1 = 6/5*lambda^2/(1+lambda)^3;
g2 = 3/25*lambda*(1 + 17*lambda - 9*lambda^2)/(1+lambda)^3;
g3 = 1/350*(5+272*lambda-831*lambda^2+1322*lambda^3-415*lambda^4)/(1+lambda)^3;
g4 = 6/5*lambda^3/(1+lambda)^3;
g5 = 3/25*lambda^2*( -4+ 17*lambda - 4*lambda^2)/(1+lambda)^3;
g6 = 1/350*lambda*(-65+832*lambda-1041*lambda^2+832*lambda^3-65*lambda^4)/(1+lambda)^3;
X_M11=g1./cita+g2*log(1./cita)+MX_xishu(l,2)+g3*cita.*log(1./cita);
X_M12=(g4./cita+g5*log(1./cita)+0.125*(1+lambda)^3*MX_xishu(l,3)+g6*cita.*log(1./cita))*8/(1+lambda)^3;
%%
load ../near-mid/MX.mat P
ln=180;

f=zeros(ln,1);
for k=1:ln
    kk=k-1;
   for q=0:kk
    f(k)=f(k)+P(200+2,200+kk-q,200+q)*lambda^(q+mod(kk,2));
   end
   f(k)=f(k)*2^(kk);
end
X_M11W=zeros(size(cita2));
X_M12W=X_M11W;
for kk=0:(ln/2-1)
%     kk=nn-1;
    X_M11W=X_M11W+f(2*kk+1)*(1+lambda)^(-2*kk)*s.^(-2*kk);
    X_M12W=X_M12W+8/(1+lambda)^3*f(2*kk+1+1)*(1+lambda)^(-2*kk-1)*s.^(-2*kk-1);
end
%%
lambda=1/lambda;
g1 = 6/5*lambda^2/(1+lambda)^3;
g2 = 3/25*lambda*(1 + 17*lambda - 9*lambda^2)/(1+lambda)^3;
g3 = 1/350*(5+272*lambda-831*lambda^2+1322*lambda^3-415*lambda^4)/(1+lambda)^3;
g4 = 6/5*lambda^3/(1+lambda)^3;
g5 = 3/25*lambda^2*( -4+ 17*lambda - 4*lambda^2)/(1+lambda)^3;
g6 = 1/350*lambda*(-65+832*lambda-1041*lambda^2+832*lambda^3-65*lambda^4)/(1+lambda)^3;
X_M22=g1./cita+g2*log(1./cita)+MX_xishu(l,5)+g3*cita.*log(1./cita);
X_M21=(g4./cita+g5*log(1./cita)+0.125*(1+lambda)^3*MX_xishu(l,4)+g6*cita.*log(1./cita))*8/(1+lambda)^3;

f=zeros(ln,1);
for k=1:ln
    kk=k-1;
   for q=0:kk
    f(k)=f(k)+P(200+2,200+kk-q,200+q)*lambda^(q+mod(kk,2));
   end
   f(k)=f(k)*2^(kk);
end
X_M22W=zeros(size(cita2));
X_M21W=X_M22W;
for kk=0:(ln/2-1)
%     kk=nn-1;
    X_M22W=X_M22W+f(2*kk+1)*(1+lambda)^(-2*kk)*s.^(-2*kk);
    X_M21W=X_M21W+8/(1+lambda)^3*f(2*kk+1+1)*(1+lambda)^(-2*kk-1)*s.^(-2*kk-1);
end
%%
X_M11=X_M11*(1+lambda)^3*5/36.d0;
X_M12=X_M12*(1+lambda)^3*5/36.d0;
X_M11W=X_M11W*(1+lambda)^3*5/36.d0;
X_M12W=X_M12W*(1+lambda)^3*5/36.d0;
X_M22=X_M22*(1+lambda)^3*5/36.d0;
X_M21=X_M21*(1+lambda)^3*5/36.d0;
X_M22W=X_M22W*(1+lambda)^3*5/36.d0;
X_M21W=X_M21W*(1+lambda)^3*5/36.d0;
%%
xi=Lb0_list;
invxi=1.0./xi;
loginvxi=log(invxi);
xiloginvxi=xi.*loginvxi;
x11a=-1.23041+0.25*invxi+1.8918*xi+9.0*loginvxi/40.0+3.0*xiloginvxi/112.0;
x12a=-x11a+0.00312-0.0011*xi;
xm=invxi/6.0+3.0*loginvxi/20.0+47.0*xiloginvxi/280.0 ...
         -0.740815+0.706802*xi;


%%
figure(1)
strlam=num2str(MX_xishu(l,1),'%5.2f');
sgtitle(['\lambda=',strlam],'FontName','Times New Roman','fontsize',15);
subplot(2,2,1)
semilogx(cita,X_M11,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,X_M11W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,X_M11ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
% xlabel('x~{\chi}_{\alpha}^{2}(3)')
ylabel('X_{11}^{M}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_mMX]);
legend({'Nearly touching','Widely separated','Far field'},'fontsize',8, ...
    'Location','best','FontName','Times New Roman')
subplot(2,2,2)
semilogx(cita,X_M12,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,X_M12W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,X_M12ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('X_{12}^{M}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_mMX]);
subplot(2,2,3)
semilogx(cita,X_M21,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,X_M21W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,X_M21ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('X_{21}^{M}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_mMX]);
subplot(2,2,4)
semilogx(cita,X_M22,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,X_M22W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,X_M22ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('X_{22}^{M}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_mMX]);

%%
figure(3)
% subplot(2,2,1)
semilogx(cita,xm,'mo','markersize',3,'MarkerFaceColor','m','linewidth',1.0),hold on;
semilogx(cita,0.5*(X_M12+X_M11)-0.5*(X_M12ff+X_M11ff),'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,0.5*(X_M12W+X_M11W)-0.5*(X_M12ff+X_M11ff),'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,0.5*(X_M12ff+X_M11ff),'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
% xlabel('x~{\chi}_{\alpha}^{2}(3)')
ylabel('X_{11}^{M}','FontName','Times New Roman','fontsize',10)
xlim([1E-6,1e1]);
ylim([0,1e1])
legend({'Nearly touching','Widely separated','Far field'},'fontsize',12, ...
    'Location','best','FontName','Times New Roman')
%%
X_M11_test1=X_M11-X_M11ff;
X_M11_test2=X_M11W-X_M11ff;
X_M11_test3=[X_M11_test1(1:24);X_M11_test2(25:50)];
X_M12_test1=X_M12-X_M12ff;
X_M12_test2=X_M12W-X_M12ff;
X_M12_test3=[X_M12_test1(1:24);X_M12_test2(25:50)];
save MX_tb X_M11_test3 X_M12_test3
X_M11_test3mf=[X_M11(1:24);X_M11W(25:50)];
X_M12_test3mf=[X_M12(1:24);X_M12W(25:50)];
save MX_tbmf X_M11_test3mf X_M12_test3mf