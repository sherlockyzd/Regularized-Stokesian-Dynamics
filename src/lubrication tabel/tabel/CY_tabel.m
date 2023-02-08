clear;clc;close all;
load FarfieldScalar.mat
CY_xishu=load('../near-mid/CY.dat');
l=1;
lambda=CY_xishu(l,1);
% cita=Lb0_list(1:32);
% cita2=Lb0_list(17:27);
% cita3=Lb0_list(20:45);
% Y_C11ff=Y11C(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% Y_C12ff=Y12C(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% Y_C22ff=Y11C(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% Y_C21ff=Y12C(Ls*(l-1)+20:Ls*(l-1)+45)./6;
cita=Lb0_list;
cita2=Lb0_list;
cita3=Lb0_list;
Y_C11ff=Y11C(Ls*(l-1)+1:Ls*l)./6;
Y_C12ff=Y12C(Ls*(l-1)+1:Ls*l)./6;
Y_C22ff=Y11C(Ls*(l-1)+1:Ls*l)./6;
Y_C21ff=Y12C(Ls*(l-1)+1:Ls*l)./6;
x_max=5E1;
s=cita2+2;
sff=cita3+2;
%%
g2 = 2/5*lambda/(1+lambda);
g3 = 1/125*(8+6*lambda+33*lambda^2)/(1+lambda);
g4 = 4/5*lambda^2/(1+lambda)^4;
g5 = 2/125*lambda*(43-24*lambda+43*lambda^2)/(1+lambda)^4;
Y_C11=g2*log(1./cita)+CY_xishu(l,2)+g3*cita.*log(1./cita);
Y_C12=(g4*log(1./cita)+CY_xishu(l,3)+g5*cita.*log(1./cita));
% Y_C11=Y_C11*2/3;
% Y_C12=Y_C12*2/3;


%%
load ../near-mid/CY.mat Q
ln=150;

f=zeros(ln,1);
for k=1:ln
    kk=k-1;
   for q=0:kk
    f(k)=f(k)+Q(200+1,200+kk-q,200+q)*lambda^(q+mod(kk,2));
   end
   f(k)=f(k)*2^(kk);
end
Y_C11W=zeros(size(cita2));
Y_C12W=Y_C11W;
for kk=0:(ln/2-1)
%     kk=nn-1;
    Y_C11W=Y_C11W+f(2*kk+1)*(1+lambda)^(-2*kk)*s.^(-2*kk);
    Y_C12W=Y_C12W+8/(1+lambda)^3*f(2*kk+2)*(1+lambda)^(-2*kk-1)*s.^(-2*kk-1);
end
%%
lambda=1/lambda;
g2 = 2/5*lambda/(1+lambda);
g3 = 1/125*(8+6*lambda+33*lambda^2)/(1+lambda);
g4 = 4/5*lambda^2/(1+lambda)^4;
g5 = 2/125*lambda*(43-24*lambda+43*lambda^2)/(1+lambda)^4;
Y_C22=g2*log(1./cita)+CY_xishu(l,5)+g3*cita.*log(1./cita);
Y_C21=(g4*log(1./cita)+CY_xishu(l,4)+g5*cita.*log(1./cita));

f=zeros(ln,1);
for k=1:ln
    kk=k-1;
   for q=0:kk
    f(k)=f(k)+Q(200+1,200+kk-q,200+q)*lambda^(q+mod(kk,2));
   end
   f(k)=f(k)*2^(kk);
end
Y_C22W=zeros(size(cita2));
Y_C21W=Y_C22W;
for kk=0:(ln/2-1)
%     kk=nn-1;
    Y_C22W=Y_C22W+f(2*kk+1)*(1+lambda)^(-2*kk)*s.^(-2*kk);
    Y_C21W=Y_C21W+8/(1+lambda)^3*f(2*kk+2)*(1+lambda)^(-2*kk-1)*s.^(-2*kk-1);
end
% alphaf=3/(1+lambda)./sff-(1+lambda^2)./(4*(1+lambda)^3)./sff.^3;
% betaf=1.5/(1+lambda)./sff+(1+lambda^2)./(2*(1+lambda)^3)./sff.^3;
% Y_C22ff=1./(1-betaf.^2);
% Y_C21ff=betaf./(betaf.^2-1);
%%
Y_C11=Y_C11*(1+lambda)^3/6.d0;
Y_C12=Y_C12*(1+lambda)^3/6.d0;
Y_C11W=Y_C11W*(1+lambda)^3/6.d0;
Y_C12W=Y_C12W*(1+lambda)^3/6.d0;
Y_C22=Y_C22*(1+lambda)^3/6.d0;
Y_C21=Y_C21*(1+lambda)^3/6.d0;
Y_C22W=Y_C22W*(1+lambda)^3/6.d0;
Y_C21W=Y_C21W*(1+lambda)^3/6.d0;
%%
figure(2)
strlam=num2str(CY_xishu(l,1),'%5.2f');
hh=sgtitle(['\lambda=',strlam]);
set(hh,'FontName','Times New Roman','fontsize',15)
subplot(2,2,1)
semilogx(cita,Y_C11,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_C11W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_C11ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Y_{11}^{C}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_max]);
legend({'Nearly touching','Widely separated','Far field'},'fontsize',8, ...
    'Location','best','FontName','Times New Roman')
subplot(2,2,2)
semilogx(cita,Y_C12,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_C12W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_C12ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Y_{12}^{C}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_max]);
subplot(2,2,3)
semilogx(cita,Y_C21,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_C21W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_C21ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Y_{21}^{C}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_max]);
subplot(2,2,4)
semilogx(cita,Y_C22,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_C22W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_C22ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Y_{22}^{C}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_max]);
%%
figure(3)
% subplot(2,2,1)
semilogx(cita,Y_C11,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_C11W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_C11ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
% xlabel('x~{\chi}_{\alpha}^{2}(3)')
ylabel('Y_{11}^{C}','FontName','Times New Roman','fontsize',10)
xlim([1E-6,5e1]);
legend({'Nearly touching','Widely separated','Far field'},'fontsize',12, ...
    'Location','best','FontName','Times New Roman')
%%
Y_C11_test1=Y_C11-Y_C11ff;
Y_C11_test2=Y_C11W-Y_C11ff;
Y_C11_test3=[Y_C11_test1(1:24);Y_C11_test2(25:50)];
Y_C12_test1=Y_C12-Y_C12ff;
Y_C12_test2=Y_C12W-Y_C12ff;
Y_C12_test3=[Y_C12_test1(1:24);Y_C12_test2(25:50)];
save CY_tb Y_C11_test3 Y_C12_test3
Y_C11_test3mf=[Y_C11(1:24);Y_C11W(25:50)];
Y_C12_test3mf=[Y_C12(1:24);Y_C12W(25:50)];
save CY_tbmf Y_C11_test3mf Y_C12_test3mf