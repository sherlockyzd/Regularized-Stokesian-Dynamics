clear;clc;
format long
load AX_tb.mat
load AY_tb.mat
load BY_tb.mat
load CX_tb.mat
load CY_tb.mat
load GX_tb.mat
load GY_tb.mat
load HY_tb.mat
load MX_tb.mat
load MY_tb.mat
load MZ_tb.mat
load FarfieldScalar.mat
Lb0_list=Lb0_list+2.0;
Lb0_list(51:end)=[];
table=[Lb0_list,X_A11_test3,X_A12_test3,Y_A11_test3,Y_A12_test3,Y_B11_test3,Y_B12_test3,X_C11_test3,X_C12_test3,...
    Y_C11_test3,Y_C12_test3,X_G11_test3,X_G12_test3,Y_G11_test3,Y_G12_test3,Y_H11_test3,Y_H12_test3,...
    X_M11_test3,X_M12_test3,Y_M11_test3,Y_M12_test3,Z_M11_test3,Z_M12_test3];

fid=fopen('lubscalar.txt','w');
str=['                  '];
fprintf(fid,[ '%s\n'],['   s^','            ','  X11A',str,'  X12A',str,'  Y11A',str,'  Y12A',str,...
    '  Y11B',str,'  Y12B',str,'  X11C',str,'  X12C',str,'  Y11C',str,'  Y12C',str,...
    '  X11G',str,'  X12G',str,'  Y11G',str,'  Y12G',str,'  Y11H',str,'  Y12H',str,'  X11M',str,...
    '  X12M',str, '  Y11M',str,'  Y12M',str,'  Z11M',str,'  Z12M']);
for ii=1:50
fprintf(fid,[ '%15.7E%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E', ...
             '%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E%24.15E\r\n'],table(ii,:));
end

% fprintf(fid,'%3d %11.4f\r\n',results);
fclose(fid);
% dlmwrite('lubscalar.txt',table,'delimiter',' ','precision',5)    %precision保留三位有效数字