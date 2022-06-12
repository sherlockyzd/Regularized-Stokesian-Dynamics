clear;clc;
format long
load AX_tbmf.mat
load AY_tbmf.mat
load BY_tbmf.mat
load CX_tbmf.mat
load CY_tbmf.mat
load GX_tbmf.mat
load GY_tbmf.mat
load HY_tbmf.mat
load MX_tbmf.mat
load MY_tbmf.mat
load MZ_tbmf.mat
load FarfieldScalar.mat
Lb0_list=Lb0_list+2.0;
Lb0_list(51:end)=[];
table=[Lb0_list,X_A11_test3mf,X_A12_test3mf,Y_A11_test3mf,Y_A12_test3mf,Y_B11_test3mf,Y_B12_test3mf,X_C11_test3mf,X_C12_test3mf,...
    Y_C11_test3mf,Y_C12_test3mf,X_G11_test3mf,X_G12_test3mf,Y_G11_test3mf,Y_G12_test3mf,Y_H11_test3mf,Y_H12_test3mf,...
    X_M11_test3mf,X_M12_test3mf,Y_M11_test3mf,Y_M12_test3mf,Z_M11_test3mf,Z_M12_test3mf];

fid=fopen('lubscalarmf.txt','w');
str=['                  '];
fprintf(fid,[ '%s\n'],['   s^','            ','  X11A',str,'  X12A',str,'  Y11A',str,'  Y12A',str,...
    '  Y11B',str,'  Y12B',str,'  X11C',str,'  X12C',str,'  Y11C',str,'  Y12C',str,...
    '  X11G',str,'  X12G',str,'  Y11G',str,'  Y12G',str,'  Y11H',str,'  Y12H',str,'  X11M',str,...
    '  X12M',str, '  Y11M',str,'  Y12M',str,'  Z11M',str,'  Z12M']);
for ii=1:45
fprintf(fid,[ '%15.7E%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E', ...
             '%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E\r\n'],table(ii,:));
end

% fprintf(fid,'%3d %11.4f\r\n',results);
fclose(fid);
% dlmwrite('lubscalar.txt',table,'delimiter',' ','precision',5)    %precision保留三位有效数字