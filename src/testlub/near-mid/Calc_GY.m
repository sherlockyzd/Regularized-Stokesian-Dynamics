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

%The output is of the form YG_11, YG_12, YG_22, YG_21

%This code relies on a recursive formula for finding the coefficients. The
%output from a perious run can be saved in a matrix GY.dat to be reused on
%future runs, greatly cutting the computational time. 
%%
clear;clc;
format long

% function [] =   Calc_GY()
warning off;
l = [1.0,0.5,0.25,0.2,0.125,0.1,0.01,2,3,4,5,10,100]; %Lambda
m = 150; %Number of terms to include in the recursive sum
Nl = length(l);
A = zeros(Nl, 4);

i = 1; %Set i to 0 if you wish to load in a matrix GY.dat from a previous run. 
if i == 0
    n = 500;
    Pcheck = zeros(n, n, n);
    P = zeros(n, n, n);
    
    Vcheck = zeros(n, n, n);
    V = zeros(n, n, n);
    
    Qcheck = zeros(n, n, n);
    Q = zeros(n, n, n);
else
    load('GY')
end





for k = 1:Nl
    lambda = l(k);
    
    g2 = 1/10*lambda*(4-lambda+7*lambda^2)/(1+lambda)^3;
    g3 = 1/500*(32 - 179*lambda+532*lambda^2-356*lambda^3+221*lambda^4)/(1+lambda)^3;
    invlambda = 1/lambda;
    invg2 = 1/10*invlambda*(4-invlambda+7*invlambda^2)/(1+invlambda)^3;
    invg3 = 1/500*(32 - 179*invlambda+532*invlambda^2-356*invlambda^3+221*invlambda^4)/(1+invlambda)^3;
    
    
    G_11Y = g2*log(4) - 2*g3;
    G_22Y = invg2*log(4) - 2*invg3;
    
    for n = 1:2:m
        n
        fm  = 0;
        fm2 = 0;
        for q = 0:n
            [temp, Vcheck, V, Pcheck, P, Qcheck, Q]  = P_AY(2, n-q, q, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            fm = fm + 3/4*temp*lambda^q;
            fm2 = fm2 + 3/4*temp*invlambda^q;
        end
        
        temp = (1+lambda)^(-n)*fm - 2/n*g2+4/n*g3/(n+2);
        G_11Y = G_11Y + temp;
        temp = (1+invlambda)^(-n)*fm2 - 2/n*invg2+4/n*invg3/(n+2);
        G_22Y = G_22Y + temp;
        
    end
    
    G_12Y =  g3;
    G_21Y =  invg3;
    for n = 2:2:m
        n
        fm = 0;
        fm2 = 0;
        for q = 0:n
            [temp, Vcheck, V, Pcheck, P, Qcheck, Q]  = P_AY(2, n-q, q, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            fm = fm + 3/4*temp*lambda^q;
            fm2 = fm2 + 3/4*temp*invlambda^q;
        end
        
        temp = (1+lambda)^(-n)*fm -2/n*g2+4/n*g3/(n+2);
        G_12Y = G_12Y - temp;
        temp = (1+invlambda)^(-n)*fm2 -2/n*invg2+4/n*invg3/(n+2);
        G_21Y = G_21Y - temp;
    end
    G_12Y = G_12Y*4/(1+lambda)^2;
    G_21Y = G_21Y*4/(1+invlambda)^2;
    
    G(k, :) = [G_11Y, G_12Y, G_21Y, G_22Y];
    
end
save('GY','P','Pcheck','V','Vcheck', 'Q', 'Qcheck');

csvwrite('GY.dat',[l',G])

% end

function [val, Vcheck, V, Pcheck, P, Qcheck, Q] = P_AY(n, p, q, val, Vcheck, V, Pcheck, P, Qcheck, Q)
if Pcheck(n+200, p+200, q+200) == 1
    val = P(n+200, p+200, q+200);
else
    if p == 0 && q == 0
        if n == 1
            val = val + 1;
        end
    else
        for s = 1:q
            [temp1, Vcheck, V, Pcheck, P, Qcheck, Q] = P_AY(s, q-s, p-n+1, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp2, Vcheck, V, Pcheck, P, Qcheck, Q] = P_AY(s, q-s, p-n-1, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp3, Vcheck, V, Pcheck, P, Qcheck, Q] = V_AY(s, q-s-2, p-n+1, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp4, Vcheck, V, Pcheck, P, Qcheck, Q] = Q_AY(s, q-s-1, p-n+1, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
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


function [val, Vcheck, V, Pcheck, P, Qcheck, Q] = V_AY(n, p, q, val, Vcheck, V, Pcheck, P, Qcheck, Q)
if Vcheck(n+200, p+200, q+200) == 1
    val = V(n+200, p+200, q+200);
else
    if p == 0 && q == 0
        if n == 1
            val = val + 1;
        end
    else
        for s = 1:q
            [temp, Vcheck, V, Pcheck, P, Qcheck, Q] = P_AY(s, q-s, p-n-1, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            val = val + nchoosek(n+s, n+1)*temp;
        end
        val = val*2*n/((n+1)*(2*n+3));
        [temp, Vcheck, V, Pcheck, P, Qcheck, Q] = P_AY(n, p, q, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
        val = val + temp;
    end
    Vcheck(n+200, p+200, q+200) = 1;
    V(n+200, p+200, q+200) = val;
end
end


function [val, Vcheck, V, Pcheck, P, Qcheck, Q] = Q_AY(n, p, q, val, Vcheck, V, Pcheck, P, Qcheck, Q)
if Qcheck(n+200, p+200, q+200) == 1
    val = Q(n+200, p+200, q+200);
else
    if p == 0 && q == 0
        val = val;
    else
        for s = 1:q
            [temp1, Vcheck, V, Pcheck, P, Qcheck, Q] = Q_AY(s, q-s-1, p-n, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp2, Vcheck, V, Pcheck, P, Qcheck, Q] = P_AY(s, q-s, p-n, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            val = val + nchoosek(n+s, n+1)*...
                (s/(n+1)*temp1 ...
                -3/(2*n*s*(n+1))*temp2);
        end
    end
    Qcheck(n+200, p+200, q+200) = 1;
    Q(n+200, p+200, q+200) = val;
end
end









