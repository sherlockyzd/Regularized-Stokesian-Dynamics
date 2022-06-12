clear;clc;close all;
load FarfieldScalar.mat
GX_xishu=load('../near-mid/GX.dat');
l=1;
lambda=GX_xishu(l,1);

% cita=Lb0_list(1:32);
% cita2=Lb0_list(17:27);
% cita3=Lb0_list(20:45);
% X_G11ff=X11G(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% X_G12ff=X12G(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% X_G22ff=X11G(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% X_G21ff=X12G(Ls*(l-1)+20:Ls*(l-1)+45)./6;
cita=Lb0_list;
cita2=Lb0_list;
cita3=Lb0_list;
X_G11ff=X11G(Ls*(l-1)+1:Ls*l)./6;
X_G12ff=X12G(Ls*(l-1)+1:Ls*l)./6;
X_G22ff=X11G(Ls*(l-1)+1:Ls*l)./6;
X_G21ff=X12G(Ls*(l-1)+1:Ls*l)./6;
x_mGX=5E1;
s=cita2+2;
sff=cita3+2;
%%
g1 = 3*lambda^2/(1+lambda)^3;
g2 = 3/10*lambda*(1+12*lambda-4*lambda^2)/(1+lambda)^3;
g3 = 1/140*(5+181*lambda-453*lambda^2+566*lambda^3-65*lambda^4)/(1+lambda)^3;
X_G11=g1./cita+g2*log(1./cita)+GX_xishu(l,2)+g3*cita.*log(1./cita);
X_G12=(-g1./cita-g2*log(1./cita)+0.25*(1+lambda)^2*GX_xishu(l,3)-g3*cita.*log(1./cita))*4/(1+lambda)^2;
%%
load ../near-mid/GX.mat P
ln=180;

f=zeros(ln,1);
for k=1:ln
    kk=k-1;
   for q=0:kk
    f(k)=f(k)+P(200+2,200+kk-q,200+q)*lambda^(q);
   end
   f(k)=0.75*f(k)*2^(kk);
end
X_G11W=zeros(size(cita2));
X_G12W=X_G11W;
for kk=0:(ln/2-1)
%     kk=nn-1;
    X_G11W=X_G11W+f(2*kk+1+1)*(1+lambda)^(-2*kk-1)*s.^(-2*kk-1);
    X_G12W=X_G12W-4/(1+lambda)^2*f(2*kk+1)*(1+lambda)^(-2*kk)*s.^(-2*kk);
end
%%
lambda=1/lambda;
g1 = 3*lambda^2/(1+lambda)^3;
g2 = 3/10*lambda*(1+12*lambda-4*lambda^2)/(1+lambda)^3;
g3 = 1/140*(5+181*lambda-453*lambda^2+566*lambda^3-65*lambda^4)/(1+lambda)^3;
X_G22=g1./cita+g2*log(1./cita)+GX_xishu(l,5)+g3*cita.*log(1./cita);
X_G21=(-g1./cita-g2*log(1./cita)+0.25*(1+lambda)^2*GX_xishu(l,4)-g3*cita.*log(1./cita))*4/(1+lambda)^2;

f=zeros(ln,1);
for k=1:ln
    kk=k-1;
   for q=0:kk
    f(k)=f(k)+P(200+2,200+kk-q,200+q)*lambda^(q);
   end
   f(k)=0.75*f(k)*2^(kk);
end
X_G22W=zeros(size(cita2));
X_G21W=X_G22W;
for kk=0:(ln/2-1)
%     kk=nn-1;
    X_G22W=X_G22W+f(2*kk+1+1)*(1+lambda)^(-2*kk-1)*s.^(-2*kk-1);
    X_G21W=X_G21W-4/(1+lambda)^2*f(2*kk+1)*(1+lambda)^(-2*kk)*s.^(-2*kk);
end
%%
X_G11=X_G11*(1+lambda)^2/6.d0;
X_G12=X_G12*(1+lambda)^2/6.d0;
X_G11W=X_G11W*(1+lambda)^2/6.d0;
X_G12W=X_G12W*(1+lambda)^2/6.d0;
X_G22=X_G22*(1+lambda)^2/6.d0;
X_G21=X_G21*(1+lambda)^2/6.d0;
X_G22W=X_G22W*(1+lambda)^2/6.d0;
X_G21W=X_G21W*(1+lambda)^2/6.d0;
%%
figure(1)
strlam=num2str(GX_xishu(l,1),'%5.2f');
hh=sgtitle(['\lambda=',strlam]);
set(hh,'FontName','Times New Roman','fontsize',15)
subplot(2,2,1)
semilogx(cita,X_G11,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,X_G11W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,X_G11ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
% xlabel('x~{\chi}_{\alpha}^{2}(3)')
ylabel('X_{11}^{G}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_mGX]);
legend({'Nearly touching','Widely separated','Far field'},'fontsize',8, ...
    'Location','best','FontName','Times New Roman')
subplot(2,2,2)
semilogx(cita,X_G12,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,X_G12W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,X_G12ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('X_{12}^{G}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_mGX]);
subplot(2,2,3)
semilogx(cita,X_G21,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,X_G21W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,X_G21ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('X_{21}^{G}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_mGX]);
subplot(2,2,4)
semilogx(cita,X_G22,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,X_G22W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,X_G22ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('X_{22}^{G}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_mGX]);

%%
figure(3)
% subplot(2,2,1)
semilogx(cita,X_G11,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,X_G11W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,X_G11ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
% xlabel('x~{\chi}_{\alpha}^{2}(3)')
ylabel('X_{11}^{G}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,1e1]);
ylim([0,1e1])
legend({'Nearly touching','Widely separated','Far field'},'fontsize',12, ...
    'Location','best','FontName','Times New Roman')
%%
X_G11_test1=X_G11-X_G11ff;
X_G11_test2=X_G11W-X_G11ff;
X_G11_test3=[X_G11_test1(1:24);X_G11_test2(25:50)];
X_G12_test1=X_G12-X_G12ff;
X_G12_test2=X_G12W-X_G12ff;
X_G12_test3=[X_G12_test1(1:24);X_G12_test2(25:50)];
save GX_tb X_G11_test3 X_G12_test3
X_G11_test3mf=[X_G11(1:24);X_G11W(25:50)];
X_G12_test3mf=[X_G12(1:24);X_G12W(25:50)];
save GX_tbmf X_G11_test3mf X_G12_test3mf