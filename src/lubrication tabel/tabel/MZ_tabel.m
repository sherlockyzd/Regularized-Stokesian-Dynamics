clear;clc;close all;
load FarfieldScalar.mat
MZ_xishu=load('../near-mid/MZ.dat');
l=1;
lambda=MZ_xishu(l,1);

% cita=Lb0_list(1:32);
% cita2=Lb0_list(17:27);
% cita3=Lb0_list(20:45);
% Z_M11ff=Z11M(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% Z_M12ff=Z12M(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% Z_M22ff=Z11M(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% Z_M21ff=Z12M(Ls*(l-1)+20:Ls*(l-1)+45)./6;
cita=Lb0_list;
cita2=Lb0_list;
cita3=Lb0_list;
Z_M11ff=Z11M(Ls*(l-1)+1:Ls*l)./6;
Z_M12ff=Z12M(Ls*(l-1)+1:Ls*l)./6;
Z_M22ff=Z11M(Ls*(l-1)+1:Ls*l)./6;
Z_M21ff=Z12M(Ls*(l-1)+1:Ls*l)./6;

Z_MMZ=5E1;
s=cita2+2;
sff=cita3+2;
%%
g3 =-3/10*(lambda^2+lambda^4)/(1+lambda)^3;
Z_M11=MZ_xishu(l,2)+g3*cita.*log(1./cita);
Z_M12=(0.125*(1+lambda)^3*MZ_xishu(l,3)-g3*cita.*log(1./cita))*8/(1+lambda)^3;
%%
load ../near-mid/MZ.mat P
ln=180;

f=zeros(ln,1);
for k=1:ln
    kk=k-1;
   for q=0:kk
    f(k)=f(k)+P(200+2,200+kk-q,200+q)*lambda^(q+mod(kk,2));
   end
   f(k)=f(k)*2^(kk);
end
Z_M11W=zeros(size(cita2));
Z_M12W=Z_M11W;
for kk=0:(ln/2-1)
%     kk=nn-1;
    Z_M11W=Z_M11W+f(2*kk+1)*(1+lambda)^(-2*kk)*s.^(-2*kk);
    Z_M12W=Z_M12W-8/(1+lambda)^3*f(2*kk+1+1)*(1+lambda)^(-2*kk-1)*s.^(-2*kk-1);
end
%%
lambda=1/lambda;
g3 =-3/10*(lambda^2+lambda^4)/(1+lambda)^3;
Z_M22=MZ_xishu(l,5)+g3*cita.*log(1./cita);
Z_M21=(0.125*(1+lambda)^3*MZ_xishu(l,4)-g3*cita.*log(1./cita))*8/(1+lambda)^3;

f=zeros(ln,1);
for k=1:ln
    kk=k-1;
   for q=0:kk
    f(k)=f(k)+P(200+2,200+kk-q,200+q)*lambda^(q+mod(kk,2));
   end
   f(k)=f(k)*2^(kk);
end
Z_M22W=zeros(size(cita2));
Z_M21W=Z_M22W;
for kk=0:(ln/2-1)
%     kk=nn-1;
    Z_M22W=Z_M22W+f(2*kk+1)*(1+lambda)^(-2*kk)*s.^(-2*kk);
    Z_M21W=Z_M21W-8/(1+lambda)^3*f(2*kk+1+1)*(1+lambda)^(-2*kk-1)*s.^(-2*kk-1);
end
%%
Z_M11=Z_M11*(1+lambda)^3*5/36.d0;
Z_M12=Z_M12*(1+lambda)^3*5/36.d0;
Z_M11W=Z_M11W*(1+lambda)^3*5/36.d0;
Z_M12W=Z_M12W*(1+lambda)^3*5/36.d0;
Z_M22=Z_M22*(1+lambda)^3*5/36.d0;
Z_M21=Z_M21*(1+lambda)^3*5/36.d0;
Z_M22W=Z_M22W*(1+lambda)^3*5/36.d0;
Z_M21W=Z_M21W*(1+lambda)^3*5/36.d0;
%%
figure(1)
strlam=num2str(MZ_xishu(l,1),'%5.2f');
hh=sgtitle(['\lambda=',strlam]);
set(hh,'FontName','Times New Roman','fontsize',15)
subplot(2,2,1)
semilogx(cita,Z_M11,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Z_M11W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Z_M11ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
% xlabel('x~{\chi}_{\alpha}^{2}(3)')
ylabel('Z_{11}^{M}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,Z_MMZ]);
legend({'Nearly touching','Widely separated','Far field'},'fontsize',8, ...
    'Location','best','FontName','Times New Roman')
subplot(2,2,2)
semilogx(cita,Z_M12,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Z_M12W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Z_M12ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Z_{12}^{M}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,Z_MMZ]);
subplot(2,2,3)
semilogx(cita,Z_M21,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Z_M21W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Z_M21ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Z_{21}^{M}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,Z_MMZ]);
subplot(2,2,4)
semilogx(cita,Z_M22,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Z_M22W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Z_M22ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Z_{22}^{M}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,Z_MMZ]);

%%
figure(3)
strlam=num2str(MZ_xishu(l,1),'%5.2f');
hh=sgtitle(['\lambda=',strlam]);
set(hh,'FontName','Times New Roman','fontsize',15)
% subplot(2,2,1)
semilogx(cita,Z_M11,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Z_M11W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Z_M11ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
% xlabel('x~{\chi}_{\alpha}^{2}(3)')
ylabel('Z_{11}^{M}','FontName','Times New Roman','fontsize',10)
xlim([1E-6,1e1]);
ylim([1,1.5])
legend({'Nearly touching','Widely separated','Far field'},'fontsize',12, ...
    'Location','best','FontName','Times New Roman')
%%
Z_M11_test1=Z_M11-Z_M11ff;
Z_M11_test2=Z_M11W-Z_M11ff;
Z_M11_test3=[Z_M11_test1(1:24);Z_M11_test2(25:50)];
Z_M12_test1=Z_M12-Z_M12ff;
Z_M12_test2=Z_M12W-Z_M12ff;
Z_M12_test3=[Z_M12_test1(1:24);Z_M12_test2(25:50)];
save MZ_tb Z_M11_test3 Z_M12_test3
Z_M11_test3mf=[Z_M11(1:24);Z_M11W(25:50)];
Z_M12_test3mf=[Z_M12(1:24);Z_M12W(25:50)];
save MZ_tbmf Z_M11_test3mf Z_M12_test3mf