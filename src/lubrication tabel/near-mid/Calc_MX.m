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

%The output is of the form XM_11, XM_12, XM_22, XM_21

%This code relies on a recursive formula for finding the coefficients. The
%output from a perious run can be saved in a matrix MX.dat to be reused on
%future runs, greatly cutting the computational time.
%%
clear;clc;
format long


% function [] =   Calc_MX()
warning off;
l = [1.0,0.5,0.25,0.2,0.125,0.1,0.01,2,3,4,5,10,20,100]; %Lambda
m = 180; %Number of terms to include in the recursive sum
Nl = length(l);
A = zeros(Nl, 4);

i = 1; %Set i to 0 if you wish to load in a matrix MX.dat from a previous run. 
if i == 0
    n = 500;
    Pcheck = zeros(n, n, n);
    P = zeros(n, n, n);
    
    Vcheck = zeros(n, n, n);
    V = zeros(n, n, n);
else
    load('MX')
end




for k = 1:Nl
    lambda = l(k);
    
    g1 = 6/5*lambda^2/(1+lambda)^3;
    g2 = 3/25*lambda*(1 + 17*lambda - 9*lambda^2)/(1+lambda)^3;
    g3 = 1/350*(5+272*lambda-831*lambda^2+1322*lambda^3-415*lambda^4)/(1+lambda)^3;
    g4 = 6/5*lambda^3/(1+lambda)^3;
    g5 = 3/25*lambda^2*( -4+ 17*lambda - 4*lambda^2)/(1+lambda)^3;
    g6 = 1/350*lambda*(-65+832*lambda-1041*lambda^2+832*lambda^3-65*lambda^4)/(1+lambda)^3;
    
    invlambda = 1/lambda;
    invg1 = 6/5*invlambda^2/(1+invlambda)^3;
    invg2 = 3/25*invlambda*(1 + 17*invlambda - 9*invlambda^2)/(1+invlambda)^3;
    invg3 = 1/350*(5+272*invlambda-831*invlambda^2+1322*invlambda^3-415*invlambda^4)/(1+invlambda)^3;
    invg4 = 6/5*invlambda^3/(1+invlambda)^3;
    invg5 = 3/25*invlambda^2*( -4+ 17*invlambda - 4*invlambda^2)/(1+invlambda)^3;
    invg6 = 1/350*invlambda*(-65+832*invlambda-1041*invlambda^2+832*invlambda^3-65*invlambda^4)/(1+invlambda)^3;
    
    
    
    M_11X = -1/4*g1 - g3 + 1;
    M_22X = -1/4*invg1 - invg3 + 1;
    for n = 2:2:m
        n
        fm  = 0;
        fm2 = 0;
        for q = 0:n
            [temp, Vcheck, V, Pcheck, P]  = P_MX(2, n-q, q, 0,  Vcheck, V, Pcheck, P);
            tempf=(1+lambda)^(q-n);
            invtempf=(1+invlambda)^(q-n);
            fm = fm + tempf*temp*(lambda/(1+lambda))^q;
            fm2 = fm2 + invtempf*temp*(invlambda/(1+invlambda))^q;
        end
        temp = fm - g1 - 2*g2/n + 4*g3/(n*(n+2));
        M_11X = M_11X + temp;
        temp = fm2 - invg1 - 2*invg2/n + 4*invg3/(n*(n+2));
        M_22X = M_22X + temp;
    end
    
    M_12X = 1/4*g4 + g5*log(4) - 2*g6;
    M_21X = 1/4*invg4 + invg5*log(4) - 2*invg6;
    
    for n = 1:2:m
        n
        fm = 0;
        fm2 = 0;
        for q = 0:n
            [temp, Vcheck, V, Pcheck, P]  = P_MX(2, n-q, q, 0,  Vcheck, V, Pcheck, P);
            tempf=(1+lambda)^(q+1-n);
            invtempf=(1+invlambda)^(q+1-n);
            fm = fm + tempf*temp*(lambda/(1+lambda))^(q+1);
            fm2 = fm2 + invtempf*temp*(invlambda/(1+invlambda))^(q+1);
        end
        temp = fm - g4 - 2*g5/n + 4*g6/(n*(n+2));
        M_12X = M_12X + temp;
        temp = fm2 - invg4 - 2*invg5/n + 4*invg6/(n*(n+2));
        M_21X = M_21X + temp;
        
    end
    M_12X = M_12X*8/(1+lambda)^3;
    M_21X = M_21X*8/(1+invlambda)^3;
    
    M(k, :) = [M_11X, M_12X, M_21X, M_22X];
    
end
% if i == 1
    save('MX','P','Pcheck','V','Vcheck');
% end

csvwrite('MX.dat',[l',M])

% end



function [val, Vcheck, V, Pcheck, P] = P_MX(n, p, q, val, Vcheck, V, Pcheck, P)
if Pcheck(n+200, p+200, q+200) == 1
    val = P(n+200, p+200, q+200);
else
    if p == 0 && q == 0
        if n == 2
            val = val + 1;
        end
    else
        for s = 1:q
            [temp1, Vcheck, V, Pcheck, P]  = P_MX(s, q-s, p-n+1, 0,  Vcheck, V, Pcheck, P);
            [temp2, Vcheck, V, Pcheck, P]  = P_MX(s, q-s, p-n-1, 0,  Vcheck, V, Pcheck, P);
            [temp3, Vcheck, V, Pcheck, P]  = V_MX(s, q-s-2, p-n+1, 0,  Vcheck, V, Pcheck, P);
            
            val = val+ nchoosek(n+s, n)/(n+1)*(...
                n*(n+0.5)*(2*n*s-n-s+2)/((2*s-1)*(n+s))*temp1 ...
                -n*(n-0.5)*temp2...
                -n*(4*n^2 -1)/(2*(2*s+1))*temp3);
        end
    end
    Pcheck(n+200, p+200, q+200) = 1;
    P(n+200, p+200, q+200) = val;
end
end


function [val, Vcheck, V, Pcheck, P] = V_MX(n, p, q, val, Vcheck, V, Pcheck, P)
if Vcheck(n+200, p+200, q+200) == 1
    val = V(n+200, p+200, q+200);
else
    if p == 0 && q == 0
        if n == 2
            val = val + 1;
        end
    else
        for s = 1:q
            [temp1, Vcheck, V, Pcheck, P]  = P_MX(s, q-s, p-n-1, 0,  Vcheck, V, Pcheck, P);
            val = val + nchoosek(n+s, n)/(n+1)*temp1;
        end
        val = -val*2*n/((2*n+3));
        [temp1, Vcheck, V, Pcheck, P]  = P_MX(n, p, q, 0,  Vcheck, V, Pcheck, P);
        val = val + temp1;
    end
    Vcheck(n+200, p+200, q+200) = 1;
    V(n+200, p+200, q+200) = val;
end
end
