clear;clc;close all;
load FarfieldScalar.mat
MY_xishu=load('../near-mid/MY.dat');
l=1;
lambda=MY_xishu(l,1);

% cita=Lb0_list(1:32);
% cita2=Lb0_list(17:27);
% cita3=Lb0_list(20:45);
% Y_M11ff=Y11M(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% Y_M12ff=Y12M(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% Y_M22ff=Y11M(Ls*(l-1)+20:Ls*(l-1)+45)./6;
% Y_M21ff=Y12M(Ls*(l-1)+20:Ls*(l-1)+45)./6;
cita=Lb0_list;
cita2=Lb0_list;
cita3=Lb0_list;
Y_M11ff=Y11M(Ls*(l-1)+1:Ls*l)./6;
Y_M12ff=Y12M(Ls*(l-1)+1:Ls*l)./6;
Y_M22ff=Y11M(Ls*(l-1)+1:Ls*l)./6;
Y_M21ff=Y12M(Ls*(l-1)+1:Ls*l)./6;

Y_MMY=5E1;
s=cita2+2;
sff=cita3+2;
%%
g2 = 6/25*lambda*(1-lambda+4*lambda^2)/(1+lambda)^3;
g3 = 1/625*(24 - 201*lambda+882*lambda^2-1182*lambda^3+591*lambda^4)/(1+lambda)^3;
g5 = 3/50*lambda^2*(7 -10*lambda +7*lambda^2)/(1+lambda)^3;
g6 = 3/2500*lambda*(221-728*lambda+1902*lambda^2-728*lambda^3+221*lambda^4)/(1+lambda)^3;
Y_M11=g2*log(1./cita)+MY_xishu(l,2)+g3*cita.*log(1./cita);
Y_M12=(g5*log(1./cita)+0.125*(1+lambda)^3*MY_xishu(l,3)+g6*cita.*log(1./cita))*8/(1+lambda)^3;
%%
load ../near-mid/MY.mat P
ln=180;

f=zeros(ln,1);
for k=1:ln
    kk=k-1;
   for q=0:kk
    f(k)=f(k)+P(200+2,200+kk-q,200+q)*lambda^(q+mod(kk,2));
   end
   f(k)=f(k)*2^(kk);
end
Y_M11W=zeros(size(cita2));
Y_M12W=Y_M11W;
for kk=0:(ln/2-1)
%     kk=nn-1;
    Y_M11W=Y_M11W+f(2*kk+1)*(1+lambda)^(-2*kk)*s.^(-2*kk);
    Y_M12W=Y_M12W+8/(1+lambda)^3*f(2*kk+1+1)*(1+lambda)^(-2*kk-1)*s.^(-2*kk-1);
end
%%
lambda=1/lambda;
g2 = 6/25*lambda*(1-lambda+4*lambda^2)/(1+lambda)^3;
g3 = 1/625*(24 - 201*lambda+882*lambda^2-1182*lambda^3+591*lambda^4)/(1+lambda)^3;
g5 = 3/50*lambda^2*(7 -10*lambda +7*lambda^2)/(1+lambda)^3;
g6 = 3/2500*lambda*(221-728*lambda+1902*lambda^2-728*lambda^3+221*lambda^4)/(1+lambda)^3;
Y_M22=g2*log(1./cita)+MY_xishu(l,5)+g3*cita.*log(1./cita);
Y_M21=(g5*log(1./cita)+0.125*(1+lambda)^3*MY_xishu(l,4)+g6*cita.*log(1./cita))*8/(1+lambda)^3;

f=zeros(ln,1);
for k=1:ln
    kk=k-1;
   for q=0:kk
    f(k)=f(k)+P(200+2,200+kk-q,200+q)*lambda^(q+mod(kk,2));
   end
   f(k)=f(k)*2^(kk);
end
Y_M22W=zeros(size(cita2));
Y_M21W=Y_M22W;
for kk=0:(ln/2-1)
%     kk=nn-1;
    Y_M22W=Y_M22W+f(2*kk+1)*(1+lambda)^(-2*kk)*s.^(-2*kk);
    Y_M21W=Y_M21W+8/(1+lambda)^3*f(2*kk+1+1)*(1+lambda)^(-2*kk-1)*s.^(-2*kk-1);
end
%%
Y_M11=Y_M11*(1+lambda)^3*5/36.d0;
Y_M12=Y_M12*(1+lambda)^3*5/36.d0;
Y_M11W=Y_M11W*(1+lambda)^3*5/36.d0;
Y_M12W=Y_M12W*(1+lambda)^3*5/36.d0;
Y_M22=Y_M22*(1+lambda)^3*5/36.d0;
Y_M21=Y_M21*(1+lambda)^3*5/36.d0;
Y_M22W=Y_M22W*(1+lambda)^3*5/36.d0;
Y_M21W=Y_M21W*(1+lambda)^3*5/36.d0;
%%
figure(1)
strlam=num2str(MY_xishu(l,1),'%5.2f');
hh=sgtitle(['\lambda=',strlam]);
set(hh,'FontName','Times New Roman','fontsize',15)
subplot(2,2,1)
semilogx(cita,Y_M11,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_M11W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_M11ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
% xlabel('x~{\chi}_{\alpha}^{2}(3)')
ylabel('Y_{11}^{M}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,Y_MMY]);
legend({'Nearly touching','Widely separated','Far field'},'fontsize',8, ...
    'Location','best','FontName','Times New Roman')
subplot(2,2,2)
semilogx(cita,Y_M12,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_M12W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_M12ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Y_{12}^{M}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,Y_MMY]);
subplot(2,2,3)
semilogx(cita,Y_M21,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_M21W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_M21ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Y_{21}^{M}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,Y_MMY]);
subplot(2,2,4)
semilogx(cita,Y_M22,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_M22W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_M22ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
ylabel('Y_{22}^{M}','FontName','Times New Roman','fontsize',10)
xlim([1E-5,Y_MMY]);

%%
figure(3)
% subplot(2,2,1)
semilogx(cita,Y_M11,'bo','markersize',3,'linewidth',1.0),hold on;
semilogx(cita2,Y_M11W,'ro','markersize',3,'MarkerFaceColor','r','linewidth',1.3),hold on;
semilogx(cita3,Y_M11ff,'ko','markersize',2,'MarkerFaceColor','k','linewidth',1)
% xlabel('x~{\chi}_{\alpha}^{2}(3)')
ylabel('Y_{11}^{M}','FontName','Times New Roman','fontsize',10)
xlim([1E-6,1e1]);
% ylim([1.3,1.5])
legend({'Nearly touching','Widely separated','Far field'},'fontsize',12, ...
    'Location','best','FontName','Times New Roman')
%%
Y_M11_test1=Y_M11-Y_M11ff;
Y_M11_test2=Y_M11W-Y_M11ff;
Y_M11_test3=[Y_M11_test1(1:24);Y_M11_test2(25:50)];
Y_M12_test1=Y_M12-Y_M12ff;
Y_M12_test2=Y_M12W-Y_M12ff;
Y_M12_test3=[Y_M12_test1(1:24);Y_M12_test2(25:50)];
save MY_tb Y_M11_test3 Y_M12_test3
Y_M11_test3mf=[Y_M11(1:24);Y_M11W(25:50)];
Y_M12_test3mf=[Y_M12(1:24);Y_M12W(25:50)];
save MY_tbmf Y_M11_test3mf Y_M12_test3mf