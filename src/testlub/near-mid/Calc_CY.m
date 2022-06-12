%Zhengdong Yu
%Institute of Hydraulics, Department of Hydraulic Engineering, Tsinghua University
%April 2022
%yuzd17@mails.tsinghua.edu.cn

%This function calculates the order one corrections to the near-field form
%of the lubrication resistance functions given in:
%D. J. Jeffrey and Y. Onishi. Calculation of the resistance and mobility
%functions for two unequal rigid spheres in low-Reynolds-number flow.
%J. Fluid Mech., 139:261, Feb 1984. doi: 10.1017/S0022112084000355.

%The user must imput the particle size ratio lambda

%The output is of the form YC_11, YC_12, YC_22, YC_21

%This code relies on a recursive formula for finding the coefficients. The
%output from a perious run can be saved in a matrix CY.dat to be reused on
%future runs, greatly cutting the computational time. 
%%

clear;clc;
% function [] =   Calc_CY()
warning off;
l = [1.0,0.5,0.25,0.2,0.125,0.1,0.01,2,3,4,5,10,100]; %Lambda
m = 150; %Number of terms to include in the recursive sum
Nl = length(l);
A = zeros(Nl, 4);

i = 0; %Set i to 0 if you wish to load in a matrix CY.dat from a previous run. 
if i == 1
    n = 500;
    Pcheck = zeros(n, n, n);
    P = zeros(n, n, n);
    
    Vcheck = zeros(n, n, n);
    V = zeros(n, n, n);
    Qcheck = zeros(n, n, n);
    Q = zeros(n, n, n);
    
    
else
    load('CY')
end



for k = 1:Nl
    lambda = l(k);
    invlambda = 1/lambda;
    
    g2 = 2/5*lambda/(1+lambda);
    g3 = 1/125*(8+6*lambda+33*lambda^2)/(1+lambda);
    g4 = 4/5*lambda^2/(1+lambda)^4;
    g5 = 2/125*lambda*(43-24*lambda+43*lambda^2)/(1+lambda)^4;
    
    invg2 = 2/5*invlambda/(1+invlambda);
    invg3 = 1/125*(8+6*invlambda+33*invlambda^2)/(1+invlambda);
    invg4 = 4/5*invlambda^2/(1+invlambda)^4;
    invg5 = 2/125*invlambda*(43-24*invlambda+43*invlambda^2)/(1+invlambda)^4;
    
    
    C_11Y = 1-g3;
    C_22Y = 1-invg3;
    for n = 2:2:m
        
        fm  = 0;
        fm2 = 0;
        for q = 0:n
            [temp, Vcheck, V, Pcheck, P, Qcheck, Q] = Q_CY(1, n-q, q, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            fm = fm + temp*lambda^(q+mod(n,2));
            fm2 = fm2 + temp*invlambda^(q+mod(n,2));
        end
        if n == 2
            m1 = -2;
        else
            m1 = (n-2);
        end
        temp = (1+lambda)^(-n)*fm - 2/n*g2+4/n*g3/(n+2);
        C_11Y = C_11Y + temp;
        temp = (1+invlambda)^(-n)*fm2 - 2/n*invg2+4/n*invg3/(n+2);
        C_22Y = C_22Y + temp;
    end
    
    C_12Y = 2*g4*log(2)-2*g5;
    C_21Y = 2*invg4*log(2)-2*invg5;
    
    for n = 1:2:m
        m1 = (n+2);
        fm = 0;
        fm2 = 0;
        for q = 0:n
            [temp, Vcheck, V, Pcheck, P, Qcheck, Q] = Q_CY(1, n-q, q, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            fm = fm + temp*lambda^(q+1);
            fm2 = fm2+ temp*invlambda^(q+1);
        end
        
        temp = 8*(1+lambda)^(-n-3)*fm - 2/n*g4+4/n*g5/m1;
        C_12Y = C_12Y + temp;
        temp = 8*(1+invlambda)^(-n-3)*fm2 - 2/n*invg4+4/n*invg5/m1;
        C_21Y = C_21Y + temp;
        
    end
%     C_12Y = -C_12Y*8/(1+lambda)^3;
%     C_21Y = -C_21Y*8/(1+invlambda)^3;
    
    C(k, :) = [C_11Y, C_12Y, C_21Y, C_22Y];
    
end
save('CY','P','Pcheck','V','Vcheck', 'Q', 'Qcheck');

csvwrite('CY.dat',[l',C])

% end

function [val, Vcheck, V, Pcheck, P, Qcheck, Q] = P_CY(n, p, q, val, Vcheck, V, Pcheck, P, Qcheck, Q)
if Pcheck(n+200, p+200, q+200) == 1
    val = P(n+200, p+200, q+200);
else
    if p == 0 && q == 0
        val = val;
    else
        for s = 1:q
            [temp1, Vcheck, V, Pcheck, P, Qcheck, Q] = P_CY(s, q-s, p-n+1, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp2, Vcheck, V, Pcheck, P, Qcheck, Q] = P_CY(s, q-s, p-n-1, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp3, Vcheck, V, Pcheck, P, Qcheck, Q] = V_CY(s, q-s-2, p-n+1, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp4, Vcheck, V, Pcheck, P, Qcheck, Q] = Q_CY(s, q-s-1, p-n+1, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            
            val = val+ nchoosek(n+s, n+1)*...
                (((2*n+1)/(2*(n+1)) * (3*(n+s)-(n*s+1)*(2*n*s-s-n+2))/(s*(n+s)*(2*s-1)))*temp1 ...
                +n*(2*n-1)/(2*(n+1))*temp2...
                +n*(4*n^2 -1)/(2*(n+1)*(2*s+1))*temp3...
                -2*(4*n^2-1)/(3*(n+1))*temp4);
        end
    end
    Pcheck(n+200, p+200, q+200) = 1;
    P(n+200, p+200, q+200) = val;
end
end


function [val, Vcheck, V, Pcheck, P, Qcheck, Q] = V_CY(n, p, q, val, Vcheck, V, Pcheck, P, Qcheck, Q)
if Vcheck(n+200, p+200, q+200) == 1
    val = V(n+200, p+200, q+200);
else
    if p == 0 && q == 0
        val = val;
    else
        for s = 1:q
            [temp1, Vcheck, V, Pcheck, P, Qcheck, Q] = P_CY(s, q-s, p-n-1, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            val = val + nchoosek(n+s, n+1)*temp1;
        end
        val = val*2*n/((n+1)*(2*n+3));
        [temp1, Vcheck, V, Pcheck, P, Qcheck, Q] = P_CY(n, p, q, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
        val = val + temp1;
    end
    Vcheck(n+200, p+200, q+200) = 1;
    V(n+200, p+200, q+200) = val;
end
end


function [val, Vcheck, V, Pcheck, P, Qcheck, Q] = Q_CY(n, p, q, val, Vcheck, V, Pcheck, P, Qcheck, Q)
if Qcheck(n+200, p+200, q+200) == 1
    val = Q(n+200, p+200, q+200);
else
    if p == 0 && q == 0
        if n == 1
            val = val+1;
        end
    else
        for s = 1:q
            [temp1, Vcheck, V, Pcheck, P, Qcheck, Q] =Q_CY(s, q-s-1, p-n, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp2, Vcheck, V, Pcheck, P, Qcheck, Q] = P_CY(s, q-s, p-n, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            
            val = val + nchoosek(n+s, n+1)*...
                (s/(n+1)*temp1 ...
                -3/(2*n*s*(n+1))*temp2);
        end
    end
    Qcheck(n+200, p+200, q+200) = 1;
    Q(n+200, p+200, q+200) = val;
end
end
