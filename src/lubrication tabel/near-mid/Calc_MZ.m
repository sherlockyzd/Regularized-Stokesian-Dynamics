%Zhengdong Yu
%Institute of Hydraulics, Department of Hydraulic Engineering, Tsinghua University
%April 2022
%yuzd17@mails.tsinghua.edu.cn

%This function calculates the order one corrections to the near-field form
%of the lubrication resistance functions given in:
%D. J. Jeffrey. The calculation of the low Reynolds number resistance
%functions for two unequal spheres. Phys. Fluids A Fluid Dyn., 4(1):16?29,
%Jan 1992. doi: 10.1063/1.858494.

%The user must imput the particle size ratio lambda

%The output is of the form ZM_11, ZM_12, ZM_22, ZM_21

%This code relies on a recursive formula for finding the coefficients. The
%output from a perious run can be saved in a matrix MZ.dat to be reused on
%future runs, greatly cutting the computational time.
%%
clc
clear all;
fclose all;
% pack
% format long

% function [] =     Calc_MZ()
warning off;
l = [1.0,0.5,0.25,0.2,0.125,0.1,0.01,2,3,4,5,10,20,100]; %Lambda
m = 180; %Number of terms to include in the recursive sum
Nl = length(l);
M = zeros(Nl, 4);

i = 0; %Set i to 0 if you wish to load in a matrix MZ.dat from a previous run. 
if i == 1
    n = 500;
    Pcheck = zeros(n, n, n);
    P = zeros(n, n, n);
    
    Qcheck = zeros(n, n, n);
    Q = zeros(n, n, n);
    
    
    Vcheck = zeros(n, n, n);
    V = zeros(n, n, n);
else
    load('MZ')
end



for k = 1:Nl
    lambda = l(k);
    g3 =-3/10*(lambda^2+lambda^4)/(1+lambda)^3;
    invlambda = 1/lambda;
    invg3 =-3/10*(invlambda^2+invlambda^4)/(1+invlambda)^3;
    
    M_11Z =  - g3 + 1;
    M_22Z = - invg3 + 1;
    for n = 2:2:m
        n
        fm  = 0;
        fm2 = 0;
        for q = 0:n
            [temp, Vcheck, V, Pcheck, P, Qcheck, Q] = P_MZ(2, n-q, q, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            tempf=(1+lambda)^(q-n);
            invtempf=(1+invlambda)^(q-n);
            fm = fm + tempf*temp*(lambda/(1+lambda))^q;
            fm2 = fm2 + invtempf*temp*(invlambda/(1+invlambda))^q;
        end
        temp = fm + 4*g3/(n*(n+2));
        M_11Z = M_11Z + temp;
        temp = fm2 + 4*invg3/(n*(n+2));
        M_22Z = M_22Z + temp;
    end
    
    M_12Z =  2*g3;
    M_21Z =  2*invg3;
    
    for n = 1:2:m
        n
        fm = 0;
        fm2 = 0;
        for q = 0:n
            [temp, Vcheck, V, Pcheck, P, Qcheck, Q] = P_MZ(2, n-q, q, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            tempf=(1+lambda)^(q+1-n);
            invtempf=(1+invlambda)^(q+1-n);
            fm = fm + tempf*temp*(lambda/(1+lambda))^(q+1);
            fm2 = fm2 + invtempf*temp*(invlambda/(1+invlambda))^(q+1);
        end
        temp = fm+ 4*g3/(n*(n+2));
        M_12Z = M_12Z - temp;
        temp = fm2+ 4*invg3/(n*(n+2));
        M_21Z = M_21Z - temp;
        
    end
    M_12Z = M_12Z*8/(1+lambda)^3;
    M_21Z = M_21Z*8/(1+invlambda)^3;
    M(k, :) = [M_11Z, M_12Z, M_21Z, M_22Z];
end

save('MZ','P','Pcheck','V','Vcheck','Q','Qcheck')

csvwrite('MZ.dat',[l',M]);

% end

function [val, Vcheck, V, Pcheck, P, Qcheck, Q] = P_MZ(n, p, q, val, Vcheck, V, Pcheck, P, Qcheck, Q)
if Pcheck(n+200, p+200, q+200) == 1
    val = P(n+200, p+200, q+200);
else
    
    if p == 0 && q == 0
        if n == 2
            val = val + 1;
        end
    else
        for s = 2:q
            [temp1, Vcheck, V, Pcheck, P, Qcheck, Q] = P_MZ(s, q-s, p-n+1, 0,  Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp2, Vcheck, V, Pcheck, P, Qcheck, Q] = P_MZ(s, q-s, p-n-1, 0,  Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp3, Vcheck, V, Pcheck, P, Qcheck, Q] = V_MZ(s, q-s-2, p-n+1, 0,  Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp4, Vcheck, V, Pcheck, P, Qcheck, Q] = Q_MZ(s, q-s-1, p-n+1, 0,  Vcheck, V, Pcheck, P, Qcheck, Q);
            val = val+ 1/(n+1)*nchoosek(n+s, n+2)*(...
                (n+0.5)*((n*s+16)*(n+s)-2*(n*s+4)*(n*s+1))/(s*(2*s-1)*(n+s))*temp1...
                +n*(n-0.5)*temp2...
                +n*(4*n^2-1)/(2*(2*s+1))*temp3...
                -2*(4*n*n-1)*temp4);
        end
    end
    Pcheck(n+200, p+200, q+200) = 1;
    P(n+200, p+200, q+200) = val;
end
end


function [val, Vcheck, V, Pcheck, P, Qcheck, Q] = V_MZ(n, p, q, val,  Vcheck, V, Pcheck, P, Qcheck, Q)
if Vcheck(n+200, p+200, q+200) == 1
    val = V(n+200, p+200, q+200);
else
    
    if p == 0 && q == 0
        if n == 2
            val = val + 1;
        end
    else
        for s = 2:q
            [temp, Vcheck, V, Pcheck, P, Qcheck, Q] = P_MZ(s, q-s, p-n-1, 0,  Vcheck, V, Pcheck, P, Qcheck, Q);
            val = val + 1/(n+1)*nchoosek(n+s, n+2)*temp;
        end
        val = val*2*n/(2*n+3);
        [temp, Vcheck, V, Pcheck, P, Qcheck, Q] = P_MZ(n, p, q, 0,  Vcheck, V, Pcheck, P, Qcheck, Q);
        val = val + temp;
    end
    Vcheck(n+200, p+200, q+200) = 1;
    V(n+200, p+200, q+200) = val;
end
end

function [val, Vcheck, V, Pcheck, P, Qcheck, Q] = Q_MZ(n, p, q, val,  Vcheck, V, Pcheck, P, Qcheck, Q)
if Qcheck(n+200, p+200, q+200) == 1
    val = Q(n+200, p+200, q+200);
else
    if p == 0 && q == 0
        val = val;
    else
        for s = 2:q
            [temp, Vcheck, V, Pcheck, P, Qcheck, Q] = Q_MZ(s, q-s-1, p-n, 0,  Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp1, Vcheck, V, Pcheck, P, Qcheck, Q] = P_MZ(s, q-s, p-n, 0,  Vcheck, V, Pcheck, P, Qcheck, Q);
            val = val + 1/(n+1)*nchoosek(n+s, n+2)*( s*temp...
                - 2/(n*s)*temp1);
        end
    end
    Qcheck(n+200, p+200, q+200) = 1;
    Q(n+200, p+200, q+200) = val;
end
end


