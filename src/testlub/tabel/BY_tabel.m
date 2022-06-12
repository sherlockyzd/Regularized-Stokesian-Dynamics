clear;clc;close all;

load FarfieldScalar.mat
BY_xishu=load('../near-mid/BY.dat');
l=1;
lambda=BY_xishu(l,1);
% [lammbdal,ss,X11A,X12A,Y11A,Y12A,Y11B,Y12B,X11C,X12C,Y11C,Y12C,X11G,X12G,...
%     Y11G,Y12G,Y11H,Y12H,X11M,X12M,Y11M,Y12M,Z11M,Z12M]=textread('tablefar.txt',...
%     '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',1);
% Ls=45;
% Lb0_list=zeros(45,1);
% Lb0_list(1:5)=[1:5]*2.0d0*1e-5;
% Lb0_list(6:10)=[1:5]*2.0d0*1e-4;
% Lb0_list(11:15)=[1:5]*2.0d0*1e-3;
% Lb0_list(16:20)=[1:5]*2.0d0*1e-2;
% Lb0_list(21:25)=[1:5]*2.0d0*1e-1;
% Lb0_list(26:45)=[1:20]*0.5d0+1.d0;
% cita=Lb0_list(1:32);
% cita2=Lb0_list(17:27);
% cita3=Lb0_list(20:45);
% Y_B11ff=Y11B(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% Y_B12ff=Y12B(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% Y_B22ff=Y11B(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% Y_B21ff=Y12B(Ls*(l-1)+20:Ls*(l-1)+45)./6;
cita=Lb0_list;
cita2=Lb0_list;
cita3=Lb0_list;
Y_B11ff=Y11B(Ls*(l-1)+1:Ls*l)./6;
Y_B12ff=Y12B(Ls*(l-1)+1:Ls*l)./6;
Y_B22ff=Y11B(Ls*(l-1)+1:Ls*l)./6;
Y_B21ff=Y12B(Ls*(l-1)+1:Ls*l)./6;
x_max=5E1;
s=cita2+2;
sff=cita3+2;
%%
g2 = -1/5*lambda*(4+lambda)/(1+lambda)^2;
g3 = -1/250*(32-33*lambda+83*lambda^2+43*lambda^3)/(1+lambda)^2;
Y_B11=g2*log(1./cita)+BY_xishu(l,2)+g3*cita.*log(1./cita);
Y_B12=-(g2*log(1./cita)-0.25*(1+lambda)^2*BY_xishu(l,3)+g3*cita.*log(1./cita))*4/(1+lambda)^2;



%%
load ../near-mid/BY.mat Q
ln=150;

f=zeros(ln,1);
for k=1:ln
    kk=k-1;
   for q=0:kk
    f(k)=f(k)+Q(200+1,200+kk-q,200+q)*lambda^q;
   end
   f(k)=f(k)*2^(kk+1);
end
Y_B11W=zeros(size(cita2));
Y_B12W=Y_B11W;
for kk=0:(ln/2-1)
%     kk=nn-1;
    Y_B11W=Y_B11W+f(2*kk+1+1)*(1+lambda)^(-2*kk-1)*s.^(-2*kk-1);
    Y_B12W=Y_B12W-4/(1+lambda)^2*f(2*kk+1)*(1+lambda)^(-2*kk)*s.^(-2*kk);
end

% betaf=1.5/(1+lambda)./sff+(1+lambda^2)./(2*(1+lambda)^3)./sff.^3;
% Y_B11ff=1./(1-betaf.^2);
% Y_B12ff=betaf./(betaf.^2-1);
%%
lambda=1/lambda;
g2 = -1/5*lambda*(4+lambda)/(1+lambda)^2;
g3 = -1/250*(32-33*lambda+83*lambda^2+43*lambda^3)/(1+lambda)^2;
Y_B22=g2*log(1./cita)+BY_xishu(l,5)+g3*cita.*log(1./cita);
Y_B21=-(g2*log(1./cita)-0.25*(1+lambda)^2*BY_xishu(l,4)+g3*cita.*log(1./cita))*4/(1+lambda)^2;


f=zeros(ln,1);
for k=1:ln
    kk=k-1;
   for q=0:kk
    f(k)=f(k)+Q(200+1,200+kk-q,200+q)*lambda^q;
   end
   f(k)=f(k)*2^(kk+1);
end
Y_B22W=zeros(size(cita2));
Y_B21W=Y_B22W;
for kk=0:(ln/2-1)
%     kk=nn-1;
    Y_B22W=Y_B22W+f(2*kk+1+1)*(1+lambda)^(-2*kk-1)*s.^(-2*kk-1);
    Y_B21W=Y_B21W-4/(1+lambda)^2*f(2*kk+1)*(1+lambda)^(-2*kk)*s.^(-2*kk);
end

% alphaf=3/(1+lambda)./sff-(1+lambda^2)./(4*(1+lambda)^3)./sff.^3;
% betaf=1.5/(1+lambda)./sff+(1+lambda^2)./(2*(1+lambda)^3)./sff.^3;
% Y_B22ff=1./(1-betaf.^2);
% Y_B21ff=betaf./(betaf.^2-1);
%%
Y_B11=Y_B11*(1+lambda)^2/6.d0;
Y_B12=Y_B12*(1+lambda)^2/6.d0;
Y_B11W=Y_B11W*(1+lambda)^2/6.d0;
Y_B12W=Y_B12W*(1+lambda)^2/6.d0;
Y_B22=Y_B22*(1+lambda)^2/6.d0;
Y_B21=Y_B21*(1+lambda)^2/6.d0;
Y_B22W=Y_B22W*(1+lambda)^2/6.d0;
Y_B21W=Y_B21W*(1+lambda)^2/6.d0;
%%
figure(2)
strlam=num2str(BY_xishu(l,1),'%5.2f');
sgtitle(['\lambda=',strlam],'FontName','Times New Roman','fontsize',15);
subplot(2,2,1)
semilogx(cita,Y_B11,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_B11W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_B11ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Y_{11}^{B}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_max]);
legend({'Nearly touching','Widely separated','Far field'},'fontsize',8, ...
    'Location','best','FontName','Times New Roman')
subplot(2,2,2)
semilogx(cita,Y_B12,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_B12W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_B12ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Y_{12}^{B}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_max]);
subplot(2,2,3)
semilogx(cita,Y_B21,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_B21W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_B21ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Y_{21}^{B}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_max]);
subplot(2,2,4)
semilogx(cita,Y_B22,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_B22W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_B22ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Y_{22}^{B}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_max]);
%%
figure(3)
% subplot(2,2,1)
semilogx(cita,Y_B11,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_B11W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_B11ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
% xlabel('x~{\chi}_{\alpha}^{2}(3)')
ylabel('Y_{11}^{B}','FontName','Times New Roman','fontsize',10)
xlim([1E-6,5e1]);
legend({'Nearly touching','Widely separated','Far field'},'fontsize',12, ...
    'Location','best','FontName','Times New Roman')
%%
Y_B11_test1=Y_B11-Y_B11ff;
Y_B11_test2=Y_B11W-Y_B11ff;
Y_B11_test3=[Y_B11_test1(1:24);Y_B11_test2(25:50)];
Y_B12_test1=Y_B12-Y_B12ff;
Y_B12_test2=Y_B12W-Y_B12ff;
Y_B12_test3=[Y_B12_test1(1:24);Y_B12_test2(25:50)];
save BY_tb Y_B11_test3 Y_B12_test3
Y_B11_test3mf=[Y_B11(1:24);Y_B11W(25:50)];
Y_B12_test3mf=[Y_B12(1:24);Y_B12W(25:50)];
save BY_tbmf Y_B11_test3mf Y_B12_test3mf