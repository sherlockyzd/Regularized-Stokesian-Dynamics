clear;clc;close all;
load FarfieldScalar.mat
CX_xishu=load('../near-mid/CX.dat');
l=1;
lambda=CX_xishu(l,1);

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
% X_C11ff=X11C(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% X_C12ff=X12C(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% X_C22ff=X11C(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% X_C21ff=X12C(Ls*(l-1)+20:Ls*(l-1)+45)./6;
cita=Lb0_list;
cita2=Lb0_list;
cita3=Lb0_list;
X_C11ff=X11C(Ls*(l-1)+1:Ls*l)./6;
X_C12ff=X12C(Ls*(l-1)+1:Ls*l)./6;
X_C22ff=X11C(Ls*(l-1)+1:Ls*l)./6;
X_C21ff=X12C(Ls*(l-1)+1:Ls*l)./6;

x_mCX=5E1;
s=cita2+2;
sff=cita3+2;
%%
X_C11=CX_xishu(l,2)-lambda^2/(4*(1+lambda))*cita.*log(1./cita);
X_C12=(CX_xishu(l,3)+lambda^2/(4*(1+lambda))*cita.*log(1./cita))*8/(1+lambda)^3;
%%
load ../near-mid/CX.mat Q
ln=150;

f=zeros(ln,1);
for k=1:ln
    kk=k-1;
   for q=0:kk
    f(k)=f(k)+Q(200+1,200+kk-q,200+q)*lambda^(q+mod(kk,2));
   end
   f(k)=f(k)*2^(kk);
end
X_C11W=zeros(size(cita2));
X_C12W=X_C11W;
for kk=0:(ln/2-1)
%     kk=nn-1;
    X_C11W=X_C11W+f(2*kk+1)*(1+lambda)^(-2*kk)*s.^(-2*kk);
    X_C12W=X_C12W-8/(1+lambda)^3*f(2*kk+2)*(1+lambda)^(-2*kk-1)*s.^(-2*kk-1);
end
%%
lambda=1/lambda;
X_C22=CX_xishu(l,5)-lambda^2/(4*(1+lambda))*cita.*log(1./cita);
X_C21=(CX_xishu(l,4)+lambda^2/(4*(1+lambda))*cita.*log(1./cita))*8/(1+lambda)^3;

f=zeros(ln,1);
for k=1:ln
    kk=k-1;
   for q=0:kk
    f(k)=f(k)+Q(200+1,200+kk-q,200+q)*lambda^(q+mod(kk,2));
   end
   f(k)=f(k)*2^(kk);
end
X_C22W=zeros(size(cita2));
X_C21W=X_C22W;
for kk=0:(ln/2-1)
%     kk=nn-1;
    X_C22W=X_C22W+f(2*kk+1)*(1+lambda)^(-2*kk)*s.^(-2*kk);
    X_C21W=X_C21W-8/(1+lambda)^3*f(2*kk+2)*(1+lambda)^(-2*kk-1)*s.^(-2*kk-1);
end
%%
X_C11=X_C11*(1+lambda)^3/6.d0;
X_C12=X_C12*(1+lambda)^3/6.d0;
X_C11W=X_C11W*(1+lambda)^3/6.d0;
X_C12W=X_C12W*(1+lambda)^3/6.d0;
X_C22=X_C22*(1+lambda)^3/6.d0;
X_C21=X_C21*(1+lambda)^3/6.d0;
X_C22W=X_C22W*(1+lambda)^3/6.d0;
X_C21W=X_C21W*(1+lambda)^3/6.d0;
%%
figure(1)
strlam=num2str(CX_xishu(l,1),'%5.2f');
sgtitle(['\lambda=',strlam],'FontName','Times New Roman','fontsize',15);
subplot(2,2,1)
semilogx(cita,X_C11,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,X_C11W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,X_C11ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
% xlabel('x~{\chi}_{\alpha}^{2}(3)')
ylabel('X_{11}^{C}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_mCX]);
legend({'Nearly touching','Widely separated','Far field'},'fontsize',8, ...
    'Location','best','FontName','Times New Roman')
subplot(2,2,2)
semilogx(cita,X_C12,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,X_C12W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,X_C12ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('X_{12}^{C}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_mCX]);
subplot(2,2,3)
semilogx(cita,X_C21,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,X_C21W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,X_C21ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('X_{21}^{C}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_mCX]);
subplot(2,2,4)
semilogx(cita,X_C22,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,X_C22W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,X_C22ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('X_{22}^{C}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,x_mCX]);

%%
figure(3)
% subplot(2,2,1)
semilogx(cita,X_C11,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,X_C11W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,X_C11ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
% xlabel('x~{\chi}_{\alpha}^{2}(3)')
ylabel('X_{11}^{C}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,1e1]);
ylim([1.3,1.5])
legend({'Nearly touching','Widely separated','Far field'},'fontsize',12, ...
    'Location','best','FontName','Times New Roman')
%%
X_C11_test1=X_C11-X_C11ff;
X_C11_test2=X_C11W-X_C11ff;
X_C11_test3=[X_C11_test1(1:24);X_C11_test2(25:50)];
X_C12_test1=X_C12-X_C12ff;
X_C12_test2=X_C12W-X_C12ff;
X_C12_test3=[X_C12_test1(1:24);X_C12_test2(25:50)];
save CX_tb X_C11_test3 X_C12_test3
X_C11_test3mf=[X_C11(1:24);X_C11W(25:50)];
X_C12_test3mf=[X_C12(1:24);X_C12W(25:50)];
save CX_tbmf X_C11_test3mf X_C12_test3mf