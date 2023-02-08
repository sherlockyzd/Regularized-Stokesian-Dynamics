clear;clc;close all;
load FarfieldScalar.mat
GY_xishu=load('../near-mid/GY.dat');
l=1;
lambda=GY_xishu(l,1);


% cita=Lb0_list(1:32);
% cita2=Lb0_list(17:27);
% cita3=Lb0_list(20:45);
% Y_G11ff=Y11G(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% Y_G12ff=Y12G(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% Y_G22ff=Y11G(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% Y_G21ff=Y12G(Ls*(l-1)+20:Ls*(l-1)+45)./6;
cita=Lb0_list;
cita2=Lb0_list;
cita3=Lb0_list;
Y_G11ff=Y11G(Ls*(l-1)+1:Ls*l)./6;
Y_G12ff=Y12G(Ls*(l-1)+1:Ls*l)./6;
Y_G22ff=Y11G(Ls*(l-1)+1:Ls*l)./6;
Y_G21ff=Y12G(Ls*(l-1)+1:Ls*l)./6;

x_mGY=5E1;
s=cita2+2;
sff=cita3+2;
%%
g2 = 1/10*lambda*(4-lambda+7*lambda^2)/(1+lambda)^3;
g3 = 1/500*(32 - 179*lambda+532*lambda^2-356*lambda^3+221*lambda^4)/(1+lambda)^3;
Y_G11=g2*log(1./cita)+GY_xishu(l,2)+g3*cita.*log(1./cita);
Y_G12=(-g2*log(1./cita)+0.25*(1+lambda)^2*GY_xishu(l,3)-g3*cita.*log(1./cita))*4/(1+lambda)^2;
%%
load ../near-mid/GY.mat P
ln=180;

f=zeros(ln,1);
for k=1:ln
    kk=k-1;
   for q=0:kk
    f(k)=f(k)+P(200+2,200+kk-q,200+q)*lambda^(q);
   end
   f(k)=0.75*f(k)*2^(kk);
end
Y_G11W=zeros(size(cita2));
Y_G12W=Y_G11W;
for kk=0:(ln/2-1)
%     kk=nn-1;
    Y_G11W=Y_G11W+f(2*kk+1+1)*(1+lambda)^(-2*kk-1)*s.^(-2*kk-1);
    Y_G12W=Y_G12W-4/(1+lambda)^2*f(2*kk+1)*(1+lambda)^(-2*kk)*s.^(-2*kk);
end
%%
lambda=1/lambda;
g2 = 1/10*lambda*(4-lambda+7*lambda^2)/(1+lambda)^3;
g3 = 1/500*(32 - 179*lambda+532*lambda^2-356*lambda^3+221*lambda^4)/(1+lambda)^3;
Y_G22=g2*log(1./cita)+GY_xishu(l,5)+g3*cita.*log(1./cita);
Y_G21=(-g2*log(1./cita)+0.25*(1+lambda)^2*GY_xishu(l,4)-g3*cita.*log(1./cita))*4/(1+lambda)^2;

f=zeros(ln,1);
for k=1:ln
    kk=k-1;
   for q=0:kk
    f(k)=f(k)+P(200+2,200+kk-q,200+q)*lambda^(q);
   end
   f(k)=0.75*f(k)*2^(kk);
end
Y_G22W=zeros(size(cita2));
Y_G21W=Y_G22W;
for kk=0:(ln/2-1)
%     kk=nn-1;
    Y_G22W=Y_G22W+f(2*kk+1+1)*(1+lambda)^(-2*kk-1)*s.^(-2*kk-1);
    Y_G21W=Y_G21W-4/(1+lambda)^2*f(2*kk+1)*(1+lambda)^(-2*kk)*s.^(-2*kk);
end
%%
Y_G11=Y_G11*(1+lambda)^2/6.d0;
Y_G12=Y_G12*(1+lambda)^2/6.d0;
Y_G11W=Y_G11W*(1+lambda)^2/6.d0;
Y_G12W=Y_G12W*(1+lambda)^2/6.d0;
Y_G22=Y_G22*(1+lambda)^2/6.d0;
Y_G21=Y_G21*(1+lambda)^2/6.d0;
Y_G22W=Y_G22W*(1+lambda)^2/6.d0;
Y_G21W=Y_G21W*(1+lambda)^2/6.d0;
%%
figure(1)
strlam=num2str(GY_xishu(l,1),'%5.2f');
hh=sgtitle(['\lambda=',strlam]);
set(hh,'FontName','Times New Roman','fontsize',15)
subplot(2,2,1)
semilogx(cita,Y_G11,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_G11W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_G11ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
% xlabel('x~{\chi}_{\alpha}^{2}(3)')
ylabel('Y_{11}^{G}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_mGY]);
legend({'Nearly touching','Widely separated','Far field'},'fontsize',8, ...
    'Location','best','FontName','Times New Roman')
subplot(2,2,2)
semilogx(cita,Y_G12,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_G12W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_G12ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Y_{12}^{G}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_mGY]);
subplot(2,2,3)
semilogx(cita,Y_G21,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_G21W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_G21ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Y_{21}^{G}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_mGY]);
subplot(2,2,4)
semilogx(cita,Y_G22,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_G22W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_G22ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Y_{22}^{G}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_mGY]);

%%
figure(3)

% subplot(2,2,1)
semilogx(cita,Y_G11,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_G11W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_G11ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
% xlabel('x~{\chi}_{\alpha}^{2}(3)')
ylabel('Y_{11}^{G}','FontName','Times New Roman','fontsize',10)
xlim([1E-6,1e1]);
% ylim([1.3,1.5])
legend({'Nearly touching','Widely separated','Far field'},'fontsize',12, ...
    'Location','best','FontName','Times New Roman')
%%
Y_G11_test1=Y_G11-Y_G11ff;
Y_G11_test2=Y_G11W-Y_G11ff;
Y_G11_test3=[Y_G11_test1(1:24);Y_G11_test2(25:50)];
Y_G12_test1=Y_G12-Y_G12ff;
Y_G12_test2=Y_G12W-Y_G12ff;
Y_G12_test3=[Y_G12_test1(1:24);Y_G12_test2(25:50)];
save GY_tb Y_G11_test3 Y_G12_test3
Y_G11_test3mf=[Y_G11(1:24);Y_G11W(25:50)];
Y_G12_test3mf=[Y_G12(1:24);Y_G12W(25:50)];
save GY_tbmf Y_G11_test3mf Y_G12_test3mf