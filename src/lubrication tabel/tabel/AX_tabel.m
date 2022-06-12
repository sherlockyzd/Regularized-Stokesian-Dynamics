clear;clc;close all;
format long
fid=fopen('tablefar.txt');
[lammbdal,ss,X11A,X12A,Y11A,Y12A,Y11B,Y12B,X11C,X12C,Y11C,Y12C,X11G,X12G,...
    Y11G,Y12G,Y11H,Y12H,X11M,X12M,Y11M,Y12M,Z11M,Z12M]=textread('../FarfieldForce/tablefar.txt',...
    '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',1);
Ls=60;
Lb0_list=zeros(Ls,1);
Lb0_list(1:10)=[1:10]*1.0d0*1e-5;
Lb0_list(11:15)=[1:5]*2.0d0*1e-4;
Lb0_list(16:20)=[1:5]*2.0d0*1e-3;
Lb0_list(21:25)=[1:5]*2.0d0*1e-2;
Lb0_list(26:50)=[1:25]*0.1d0+0.1d0;
Lb0_list(51:60)=[1:10]*0.5d0+2.6d0;

save FarfieldScalar.mat
load FarfieldScalar.mat

AX_xishu=load('../near-mid/AX.dat');

l=1;
lambda=AX_xishu(l,1);

% cita=Lb0_list(1:32);
% cita2=Lb0_list(17:27);
% cita3=Lb0_list(20:45);
cita=Lb0_list;
cita2=Lb0_list;
cita3=Lb0_list;
% X_A11ff=X11A(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% X_A12ff=X12A(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% X_A22ff=X11A(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% X_A21ff=X12A(Ls*(l-1)+20:Ls*(l-1)+45)./6;
X_A11ff=X11A(Ls*(l-1)+1:Ls*(l-1)+Ls)./6;
X_A12ff=X12A(Ls*(l-1)+1:Ls*(l-1)+Ls)./6;
X_A22ff=X11A(Ls*(l-1)+1:Ls*(l-1)+Ls)./6;
X_A21ff=X12A(Ls*(l-1)+1:Ls*(l-1)+Ls)./6;


X_Aax=5E1;
s=cita2+2;
sff=cita3+2;
%%
g1 = 2*lambda^2/(1+lambda)^3;
g2 = 1/5*lambda*(1+7*lambda+lambda^2)/(1+lambda)^3;
g3 = 1/42*(1+18*lambda-29*lambda^2+18*lambda^3+lambda^4)/(1+lambda)^3;
X_A11=g1./cita+g2*log(1./cita)+AX_xishu(l,2)+g3*cita.*log(1./cita);
X_A12=-2/(1+lambda)*(g1./cita+g2*log(1./cita)-0.5*(1+lambda)*AX_xishu(l,3)+g3*cita.*log(1./cita));
%%
load ../near-mid/AX.mat P
ln=150;
f=zeros(ln,1);
for k=1:ln
   for q=1:k
    f(k)=f(k)+P(200+1,200+k-q,200+q-1)*lambda^(q-1);
   end
   f(k)=f(k)*2^(k-1);
end

X_A11W=zeros(size(cita2));
X_A12W=X_A11W;
for nn=1:ln/2
    kk=nn-1;
    X_A11W=X_A11W+f(2*kk+1)*(1+lambda)^(-2*kk)*s.^(-2*kk);
    X_A12W=X_A12W-2/(1+lambda)*f(2*kk+2)*(1+lambda)^(-2*kk-1)*s.^(-2*kk-1);
end
alphaf=3/(1+lambda)./sff-(1+lambda^2)./(4*(1+lambda)^3)./sff.^3;
X_A11ff0=1./(1-alphaf.^2);
X_A12ff0=alphaf./(alphaf.^2-1);

%%
lambda=1/lambda;
g1 = 2*lambda^2/(1+lambda)^3;
g2 = 1/5*lambda*(1+7*lambda+lambda^2)/(1+lambda)^3;
g3 = 1/42*(1+18*lambda-29*lambda^2+18*lambda^3+lambda^4)/(1+lambda)^3;
X_A22=g1./cita+g2*log(1./cita)+AX_xishu(l,5)+g3*cita.*log(1./cita);
X_A21=-2/(1+lambda)*(g1./cita+g2*log(1./cita)-0.5*(1+lambda)*AX_xishu(l,4)+g3*cita.*log(1./cita));

f=zeros(ln,1);
for k=1:ln
   for q=1:k
    f(k)=f(k)+P(200+1,200+k-q,200+q-1)*lambda^(q-1);
   end
   f(k)=f(k)*2^(k-1);
end

X_A22W=zeros(size(cita2));
X_A21W=X_A22W;
for nn=1:ln/2
    kk=nn-1;
    X_A22W=X_A22W+f(2*kk+1)*(1+lambda)^(-2*kk)*s.^(-2*kk);
    X_A21W=X_A21W-2/(1+lambda)*f(2*kk+2)*(1+lambda)^(-2*kk-1)*s.^(-2*kk-1);
end


%%
X_A11=X_A11*(1+lambda)/2.d0;
X_A12=X_A12*(1+lambda)/2.d0;
X_A11W=X_A11W*(1+lambda)/2.d0;
X_A12W=X_A12W*(1+lambda)/2.d0;
X_A22=X_A22*(1+lambda)/2.d0;
X_A21=X_A21*(1+lambda)/2.d0;
X_A22W=X_A22W*(1+lambda)/2.d0;
X_A21W=X_A21W*(1+lambda)/2.d0;

%%
figure(1)
strlam=num2str(AX_xishu(l,1),'%5.2f');
% hh=suptitle(['\lambda=',strlam]);
% set(hh,'FontName','Times New Roman','fontsize',15)
sgtitle(['\lambda=',strlam],'FontName','Times New Roman','fontsize',15);
% set(hh,)
subplot(2,2,1)
semilogx(cita,X_A11,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,X_A11W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,X_A11ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
% xlabel('x~{\chi}_{\alpha}^{2}(3)')
ylabel('X_{11}^{A}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,X_Aax]);
legend({'Nearly touching','Widely separated','Far field'},'fontsize',8, ...
    'Location','best','FontName','Times New Roman')
subplot(2,2,2)
semilogx(cita,X_A12,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,X_A12W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,X_A12ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('X_{12}^{A}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,X_Aax]);
subplot(2,2,3)
semilogx(cita,X_A21,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,X_A21W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,X_A21ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('X_{21}^{A}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,X_Aax]);
subplot(2,2,4)
semilogx(cita,X_A22,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,X_A22W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,X_A22ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('X_{22}^{A}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,X_Aax]);

%%
figure(3)
% subplot(2,2,1)
semilogx(cita,X_A11,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,X_A11W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,X_A11ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
% xlabel('x~{\chi}_{\alpha}^{2}(3)')
ylabel('X_{11}^{A}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,1e1]);
ylim([0,1e2]);
legend({'Nearly touching','Widely separated','Far field'},'fontsize',12, ...
    'Location','best','FontName','Times New Roman')
%%
X_A11_test1=X_A11-X_A11ff;
X_A11_test2=X_A11W-X_A11ff;
X_A11_test3=[X_A11_test1(1:24);X_A11_test2(25:50)];
X_A12_test1=X_A12-X_A12ff;
X_A12_test2=X_A12W-X_A12ff;
X_A12_test3=[X_A12_test1(1:24);X_A12_test2(25:50)];
save AX_tb X_A11_test3 X_A12_test3

X_A11_test3mf=[X_A11(1:24);X_A11W(25:50)];
X_A12_test3mf=[X_A12(1:24);X_A12W(25:50)];
save AX_tbmf X_A11_test3mf X_A12_test3mf

