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
m = 500; %Number of terms to include in the recursive sum
Nl = length(l);
A = zeros(Nl, 4);

for k = 1:Nl
    lambda = l(k);
    invlambda = 1/lambda;
    
    g1 = lambda^3*(1+lambda)^(-3);
    g2 = -lambda^3*(1+lambda)^(-3);
    invg1 = invlambda^3*(1+invlambda)^(-3);
    invg2 =  -invlambda^3*(1+invlambda)^(-3);

    C_11Y = g1*HZ(3,lambda/(1+lambda),m);
    C_22Y = invg1*HZ(3,invlambda/(1+invlambda),m);
    
    C_12Y = g2*HZ(3,1,m);
    C_21Y = invg2*HZ(3,1,m);
 
    C(k, :) = [C_11Y, C_12Y, C_21Y, C_22Y];
    
end
% save('CY','P','Pcheck','V','Vcheck', 'Q', 'Qcheck');

csvwrite('CX.dat',[l',C])
clear;
m = 150; %Number of terms to include in the recursive sum

i = 0; %Set i to 0 if you wish to load in a matrix CY.dat from a previous run. 
if i == 1
    n = 500;
    Qcheck = zeros(n, n, n);
    Q = zeros(n, n, n);
else
    load('CX')
end

for n=1:m
    n
    for q = 0:n
        [temp, Qcheck, Q] = Q_CX(1, n-q, q, 0,Qcheck, Q);
    end
end
save('CX', 'Q', 'Qcheck');



%%


function val =HZ(z,a,m)
val=0;
  for k = 0:m
       val=val+(k+a)^(-z);
  end
end


function [val,Qcheck, Q] = Q_CX(n, p, q, val,Qcheck, Q)
if Qcheck(n+200, p+200, q+200) == 1
    val = Q(n+200, p+200, q+200);
else
    if p == 0 && q == 0
        if n == 1
            val = val+1;
        end
    else
        for s = 0:q
            [temp1, Qcheck, Q] =Q_CX(s, q-s-1, p-n, 0,  Qcheck, Q);
%             [temp2, Vcheck, V, Pcheck, P, Qcheck, Q] = P_CY(s, q-s, p-n, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            
            val = val + nchoosek(n+s, n)*s/(n+1)*temp1 ;
        end
    end
    Qcheck(n+200, p+200, q+200) = 1;
    Q(n+200, p+200, q+200) = val;
end
end