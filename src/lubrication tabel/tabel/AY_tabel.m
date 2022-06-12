clear;clc;close all;

load FarfieldScalar.mat
AY_xishu=load('../near-mid/AY.dat');
l=1;
lambda=AY_xishu(l,1);

% Ls=50;
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
cita=Lb0_list;
cita2=Lb0_list;
cita3=Lb0_list;
% Y_A11ff=Y11A(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% Y_A12ff=Y12A(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% Y_A22ff=Y11A(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% Y_A21ff=Y12A(Ls*(l-1)+20:Ls*(l-1)+45)./6;
Y_A11ff=Y11A(Ls*(l-1)+1:Ls*l)./6;
Y_A12ff=Y12A(Ls*(l-1)+1:Ls*l)./6;
Y_A22ff=Y11A(Ls*(l-1)+1:Ls*l)./6;
Y_A21ff=Y12A(Ls*(l-1)+1:Ls*l)./6;
x_max=5E1;
s=cita2+2;
sff=cita3+2;
%%
g2y = 4/15*lambda*(2+lambda+2*lambda^2)/(1+lambda)^3;
g3y = 2/375*(16-45*lambda+58*lambda^2-45*lambda^3+16*lambda^4)/(1+lambda)^3;
Y_A11=g2y*log(1./cita)+AY_xishu(l,2)+g3y*cita.*log(1./cita);
Y_A12=-(g2y*log(1./cita)-0.5*(1+lambda)*AY_xishu(l,3)+g3y*cita.*log(1./cita))*2/(1+lambda);



%%
load  ../near-mid/AY.mat P
ln=150;

f=zeros(ln,1);
for k=1:ln
    kk=k-1;
   for q=0:kk
    f(k)=f(k)+P(200+1,200+kk-q,200+q)*lambda^q;
   end
   f(k)=f(k)*2^kk;
end
Y_A11W=zeros(size(cita2));
Y_A12W=Y_A11W;
for nn=1:ln/2
    kk=nn-1;
    Y_A11W=Y_A11W+f(2*kk+1)*(1+lambda)^(-2*kk)*s.^(-2*kk);
    Y_A12W=Y_A12W-2/(1+lambda)*f(2*nn)*(1+lambda)^(-2*kk-1)*s.^(-2*kk-1);
end
% 
% betaf=1.5/(1+lambda)./sff+(1+lambda^2)./(2*(1+lambda)^3)./sff.^3;
% Y_A11ff=1./(1-betaf.^2);
% Y_A12ff=betaf./(betaf.^2-1);
%%
lambda=1/lambda;
g2y = 4/15*lambda*(2+lambda+2*lambda^2)/(1+lambda)^3;
g3y = 2/375*(16-45*lambda+58*lambda^2-45*lambda^3+16*lambda^4)/(1+lambda)^3;
Y_A22=g2y*log(1./cita)+AY_xishu(l,5)+g3y*cita.*log(1./cita);
Y_A21=-(g2y*log(1./cita)-0.5*(1+lambda)*AY_xishu(l,4)+g3y*cita.*log(1./cita))*2/(1+lambda);
f=zeros(ln,1);
for k=1:ln
    kk=k-1;
   for q=0:kk
    f(k)=f(k)+P(200+1,200+kk-q,200+q)*lambda^q;
   end
   f(k)=f(k)*2^kk;
end
Y_A22W=zeros(size(cita2));
Y_A21W=Y_A22W;
for nn=1:ln/2
    kk=nn-1;
    Y_A22W=Y_A22W+f(2*kk+1)*(1+lambda)^(-2*kk)*s.^(-2*kk);
    Y_A21W=Y_A21W-2/(1+lambda)*f(2*nn)*(1+lambda)^(-2*kk-1)*s.^(-2*kk-1);
end
% 
% alphaf=3/(1+lambda)./sff-(1+lambda^2)./(4*(1+lambda)^3)./sff.^3;
% betaf=1.5/(1+lambda)./sff+(1+lambda^2)./(2*(1+lambda)^3)./sff.^3;
% Y_A22ff=1./(1-betaf.^2);
% Y_A21ff=betaf./(betaf.^2-1);
%%
Y_A11=Y_A11*(1+lambda)/2.d0;
Y_A12=Y_A12*(1+lambda)/2.d0;
Y_A11W=Y_A11W*(1+lambda)/2.d0;
Y_A12W=Y_A12W*(1+lambda)/2.d0;
Y_A22=Y_A22*(1+lambda)/2.d0;
Y_A21=Y_A21*(1+lambda)/2.d0;
Y_A22W=Y_A22W*(1+lambda)/2.d0;
Y_A21W=Y_A21W*(1+lambda)/2.d0;
%%
figure(2)
strlam=num2str(AY_xishu(l,1),'%5.2f');
sgtitle(['\lambda=',strlam],'FontName','Times New Roman','fontsize',15);
subplot(2,2,1)
semilogx(cita,Y_A11,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_A11W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_A11ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Y_{11}^{A}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_max]);
legend({'Nearly touching','Widely separated','Far field'},'fontsize',8, ...
    'Location','best','FontName','Times New Roman')
subplot(2,2,2)
semilogx(cita,Y_A12,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_A12W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_A12ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Y_{12}^{A}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_max]);
subplot(2,2,3)
semilogx(cita,Y_A21,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_A21W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_A21ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Y_{21}^{A}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_max]);
subplot(2,2,4)
semilogx(cita,Y_A22,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_A22W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_A22ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Y_{22}^{A}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_max]);
%%
figure(3)
% subplot(2,2,1)
semilogx(cita,Y_A11,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_A11W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_A11ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
% xlabel('x~{\chi}_{\alpha}^{2}(3)')
ylabel('Y_{11}^{A}','FontName','Times New Roman','fontsize',10)
xlim([5E-7,5e1]);
legend({'Nearly touching','Widely separated','Far field'},'fontsize',12, ...
    'Location','best','FontName','Times New Roman')
%%
Y_A11_test1=Y_A11-Y_A11ff;
Y_A11_test2=Y_A11W-Y_A11ff;
Y_A11_test3=[Y_A11_test1(1:24);Y_A11_test2(25:50)];
Y_A12_test1=Y_A12-Y_A12ff;
Y_A12_test2=Y_A12W-Y_A12ff;
Y_A12_test3=[Y_A12_test1(1:24);Y_A12_test2(25:50)];
save AY_tb Y_A11_test3 Y_A12_test3
Y_A11_test3mf=[Y_A11(1:24);Y_A11W(25:50)];
Y_A12_test3mf=[Y_A12(1:24);Y_A12W(25:50)];
save AY_tbmf Y_A11_test3mf Y_A12_test3mf