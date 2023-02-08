clear;clc;close all;
load FarfieldScalar.mat
HY_xishu=load('../near-mid/HY.dat');
l=1;
lambda=HY_xishu(l,1);

% cita=Lb0_list(1:32);
% cita2=Lb0_list(17:27);
% cita3=Lb0_list(20:45);
% Y_H11ff=Y11H(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% Y_H12ff=Y12H(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% Y_H22ff=Y11H(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% Y_H21ff=Y12H(Ls*(l-1)+20:Ls*(l-1)+45)./6;
cita=Lb0_list;
cita2=Lb0_list;
cita3=Lb0_list;
Y_H11ff=Y11H(Ls*(l-1)+1:Ls*l)./6;
Y_H12ff=Y12H(Ls*(l-1)+1:Ls*l)./6;
Y_H22ff=Y11H(Ls*(l-1)+1:Ls*l)./6;
Y_H21ff=Y12H(Ls*(l-1)+1:Ls*l)./6;

x_mHY=5E1;
s=cita2+2;
sff=cita3+2;
%%
g2 = 1/10*lambda*(2-lambda)/(1+lambda)^2;
g3 = 1/500*(16-61*lambda + 180*lambda^2 + 2*lambda^3)/(1+lambda)^2;
g5 = 1/20*(lambda^2 + 7*lambda^3) / (1+lambda)^2;
g6 = 1/1000*(43*lambda + 147*lambda^2 - 185*lambda^3 +221*lambda^4)/(1+lambda)^2;
Y_H11=g2*log(1./cita)+HY_xishu(l,2)+g3*cita.*log(1./cita);
Y_H12=(g5*log(1./cita)+0.125*(1+lambda)^3*HY_xishu(l,3)+g6*cita.*log(1./cita))*8/(1+lambda)^3;
%%
load ../near-mid/HY.mat P
ln=180;

f=zeros(ln,1);
for k=1:ln
    kk=k-1;
   for q=0:kk
    f(k)=f(k)+P(200+2,200+kk-q,200+q)*lambda^(q+mod(kk,2));
   end
   f(k)=-3/8*f(k)*2^(kk);
end
Y_H11W=zeros(size(cita2));
Y_H12W=Y_H11W;
for kk=0:(ln/2-1)
%     kk=nn-1;
    Y_H11W=Y_H11W+f(2*kk+1)*(1+lambda)^(-2*kk)*s.^(-2*kk);
    Y_H12W=Y_H12W+8/(1+lambda)^3*f(2*kk+1+1)*(1+lambda)^(-2*kk-1)*s.^(-2*kk-1);
end
%%
lambda=1/lambda;
g2 = 1/10*lambda*(2-lambda)/(1+lambda)^2;
g3 = 1/500*(16-61*lambda + 180*lambda^2 + 2*lambda^3)/(1+lambda)^2;
g5 = 1/20*(lambda^2 + 7*lambda^3) / (1+lambda)^2;
g6 = 1/1000*(43*lambda + 147*lambda^2 - 185*lambda^3 +221*lambda^4)/(1+lambda)^2;
Y_H22=g2*log(1./cita)+HY_xishu(l,5)+g3*cita.*log(1./cita);
Y_H21=(g5*log(1./cita)+0.125*(1+lambda)^3*HY_xishu(l,4)+g6*cita.*log(1./cita))*8/(1+lambda)^3;

f=zeros(ln,1);
for k=1:ln
    kk=k-1;
   for q=0:kk
    f(k)=f(k)+P(200+2,200+kk-q,200+q)*lambda^(q+mod(kk,2));
   end
   f(k)=-3/8*f(k)*2^(kk);
end
Y_H22W=zeros(size(cita2));
Y_H21W=Y_H22W;
for kk=0:(ln/2-1)
%     kk=nn-1;
    Y_H22W=Y_H22W+f(2*kk+1)*(1+lambda)^(-2*kk)*s.^(-2*kk);
    Y_H21W=Y_H21W+8/(1+lambda)^3*f(2*kk+1+1)*(1+lambda)^(-2*kk-1)*s.^(-2*kk-1);
end
%%
Y_H11=Y_H11*(1+lambda)^3/6.d0;
Y_H12=Y_H12*(1+lambda)^3/6.d0;
Y_H11W=Y_H11W*(1+lambda)^3/6.d0;
Y_H12W=Y_H12W*(1+lambda)^3/6.d0;
Y_H22=Y_H22*(1+lambda)^3/6.d0;
Y_H21=Y_H21*(1+lambda)^3/6.d0;
Y_H22W=Y_H22W*(1+lambda)^3/6.d0;
Y_H21W=Y_H21W*(1+lambda)^3/6.d0;
%%
figure(1)
strlam=num2str(HY_xishu(l,1),'%5.2f');
hh=sgtitle(['\lambda=',strlam]);
set(hh,'FontName','Times New Roman','fontsize',15)

subplot(2,2,1)
semilogx(cita,Y_H11,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_H11W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_H11ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
% xlabel('x~{\chi}_{\alpha}^{2}(3)')
ylabel('Y_{11}^{H}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_mHY]);
legend({'Nearly touching','Widely separated','Far field'},'fontsize',8, ...
    'Location','best','FontName','Times New Roman')
subplot(2,2,2)
semilogx(cita,Y_H12,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_H12W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_H12ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Y_{12}^{H}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_mHY]);
subplot(2,2,3)
semilogx(cita,Y_H21,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_H21W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_H21ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Y_{21}^{H}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_mHY]);
subplot(2,2,4)
semilogx(cita,Y_H22,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_H22W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_H22ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Y_{22}^{H}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_mHY]);

%%
figure(3)
% subplot(2,2,1)
semilogx(cita,Y_H11,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_H11W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_H11ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
% xlabel('x~{\chi}_{\alpha}^{2}(3)')
ylabel('Y_{11}^{H}','FontName','Times New Roman','fontsize',10)
xlim([1E-6,1e1]);
% ylim([1.3,1.5])
legend({'Nearly touching','Widely separated','Far field'},'fontsize',12, ...
    'Location','best','FontName','Times New Roman')
%%
Y_H11_test1=Y_H11-Y_H11ff;
Y_H11_test2=Y_H11W-Y_H11ff;
Y_H11_test3=[Y_H11_test1(1:24);Y_H11_test2(25:50)];
Y_H12_test1=Y_H12-Y_H12ff;
Y_H12_test2=Y_H12W-Y_H12ff;
Y_H12_test3=[Y_H12_test1(1:24);Y_H12_test2(25:50)];
save HY_tb Y_H11_test3 Y_H12_test3
Y_H11_test3mf=[Y_H11(1:24);Y_H11W(25:50)];
Y_H12_test3mf=[Y_H12(1:24);Y_H12W(25:50)];
save HY_tbmf Y_H11_test3mf Y_H12_test3mf