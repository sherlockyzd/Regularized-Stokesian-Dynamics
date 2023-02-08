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

%The output is of the form XA_11, XA_12, XA_22, XA_21

%This code relies on a recursive formula for finding the coefficients. The
%output from a perious run can be saved in a matrix AX.dat to be reused on
%future runs, greatly cutting the computational time. 
clear;clc;
% function [] = Calc_AX()
format long
warning off;

l = [1.0,0.5,0.25,0.2,0.125,0.1,0.01,2,3,4,5,10,100]; %Lambda
m = 150; %Number of terms to include in the recursive sum

Nl = length(l);
A = zeros(Nl, 4);

i = 0; %Set i to 0 if you wish to load in a matrix AX.dat from a previous run. 
if i == 1
    n = 500;
    Pcheck = zeros(n, n, n);
    P = zeros(n, n, n);
    
    Vcheck = zeros(n, n, n);
    V = zeros(n, n, n);
else
    load('AX.mat');
end
%%


for k = 1:Nl
    lambda = l(k);
    g1 = 2*lambda^2/(1+lambda)^3;
    g2 = 1/5*lambda*(1+7*lambda+lambda^2)/(1+lambda)^3;
    g3 = 1/42*(1+18*lambda-29*lambda^2+18*lambda^3+lambda^4)/(1+lambda)^3;
    A_11X = 1 - 1/4*g1;
    for n = 2:2:m
        n
        fm  = 0;
        tempf=(1+lambda)^(-n);
        for q = 0:n
            [temp, Vcheck, V, Pcheck, P]  = P_A(1, n-q, q, 0, Vcheck, V, Pcheck, P);
            fm = fm + tempf*temp*lambda^q;
        end
        f = fm;
        if n == 2
            m1 = -2;
        else
            m1 = (n-2);
        end
        temp = f - g1-2/n*g2+4/n*g3/m1;
        A_11X = A_11X + temp;
    end
    
    A_12X = 1/4*g1 + 2*g2*log(2)+2*g3;
    for n = 1:2:m
        n
        m1 = (n-2);
        fm = 0;
        tempf=(1+lambda)^(-n);
        for q = 0:n
            [temp, Vcheck, V, Pcheck, P]  = P_A(1, n-q, q, 0, Vcheck, V, Pcheck, P);
            fm = fm +tempf*temp*lambda^q;
        end
        f = fm;
        
        temp =f - g1-2/n*g2+4/n*g3/m1;
        A_12X = A_12X + temp;
        
    end
    A_12X = A_12X/(-.5*(1+lambda));
    
    
    lambda = 1/lambda;
    g1 = 2*lambda^2/(1+lambda)^3;
    g2 = 1/5*lambda*(1+7*lambda+lambda^2)/(1+lambda)^3;
    g3 = 1/42*(1+18*lambda-29*lambda^2+18*lambda^3+lambda^4)/(1+lambda)^3;
    A_22X = 1 - 1/4*g1;
    for n = 2:2:m
        n
        fm  = 0;
        tempf=(1+lambda)^(-n);
        for q = 0:n
            [temp, Vcheck, V, Pcheck, P]  = P_A(1, n-q, q, 0, Vcheck, V, Pcheck, P);
            
            fm = fm + tempf*temp*lambda^q;
        end
        f = fm;
        if n == 2
            m1 = -2;
        else
            m1 = (n-2);
        end
        temp = f - g1-2/n*g2+4/n*g3/m1;
        A_22X = A_22X + temp;
    end
    
    A_21X = 1/4*g1 + 2*g2*log(2)+2*g3;
    for n = 1:2:m
        n
        m1 = (n-2);
        fm = 0;
        tempf=(1+lambda)^(-n);
        for q = 0:n
            
            [temp, Vcheck, V, Pcheck, P]  = P_A(1, n-q, q, 0, Vcheck, V, Pcheck, P);
            
            fm = fm + tempf*temp*lambda^q;
        end
        f = fm;
        
        temp = f - g1-2/n*g2+4/n*g3/m1;
        A_21X = A_21X + temp;
        
    end
    A_21X = A_21X/(-.5*(1+lambda));
    A(k, :) = [A_11X, A_12X,A_21X,A_22X];
    
end
save('AX','P','Pcheck','V','Vcheck');

csvwrite('AX.dat',[l',A])


% end

%%
function [val, Vcheck, V, Pcheck, P] = P_A(n, p, q, val, Vcheck, V, Pcheck, P)
if Pcheck(n+200, p+200, q+200) == 1
    val = P(n+200, p+200, q+200);
else
    if p == 0 && q == 0
        if n == 1
            val = val + 1;
        else
            val = val;
        end
    elseif q <= 0
        val = val;
    else
        for s = 1:q
            [temp1, Vcheck, V, Pcheck, P] = P_A(s, q-s, p-n+1, 0, Vcheck, V, Pcheck, P);
            [temp2, Vcheck, V, Pcheck, P] = P_A(s, q-s, p-n-1, 0, Vcheck, V, Pcheck, P);
            [temp3, Vcheck, V, Pcheck, P] = V_A(s, q-s-2, p-n+1, 0, Vcheck, V, Pcheck, P);
            val = val+ nchoosek(n+s, n)*(n*(2*n+1)*(2*n*s-n-s+2)/(2*(n+1)*(2*s-1)*(n+s))*temp1 ...
                -n*(2*n-1)/(2*(n+1))*temp2...
                -n*(4*n^2 -1)/(2*(n+1)*(2*s+1))*temp3);
        end
    end
    Pcheck(n+200, p+200, q+200) = 1;
    P(n+200, p+200, q+200) = val;
end
end


function [val, Vcheck, V, Pcheck, P] = V_A(n, p, q, val, Vcheck, V, Pcheck, P)
if Vcheck(n+200, p+200, q+200) == 1
    val = V(n+200, p+200, q+200);
else
    
    if p == 0 && q == 0
        if n == 1
            val = val + 1;
        end
    else
        for s = 1:q
            [temp, Vcheck, V, Pcheck, P]  = P_A(s, q-s, p-n-1, 0, Vcheck, V, Pcheck, P);
            val = val + nchoosek(n+s, n)*temp;
        end
        val = -val*2*n/((n+1)*(2*n+3));
        [temp, Vcheck, V, Pcheck, P]  = P_A(n, p, q, 0, Vcheck, V, Pcheck, P);
        val = val + temp;
    end
    Vcheck(n+200, p+200, q+200) = 1;
    V(n+200, p+200, q+200) = val;
end
end
