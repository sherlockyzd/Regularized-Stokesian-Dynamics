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

%The output is of the form YH_11, YH_12, YH_22, YH_21

%This code relies on a recursive formula for finding the coefficients. The
%output from a perious run can be saved in a matrix HY.dat to be reused on
%future runs, greatly cutting the computational time. 
%%
clear;clc;
format long
% function [] =   Calc_HY()
warning off;
l = [1.0,0.5,0.25,0.2,0.125,0.1,0.01,2,3,4,5,10,20,100]; %Lambda
m = 180; %Number of terms to include in the recursive sum
Nl = length(l);
A = zeros(Nl, 4);

i = 0; %Set i to 0 if you wish to load in a matrix HY.dat from a previous run. 
if i == 1
    n = 500;
    Pcheck = zeros(n, n, n);
    P = zeros(n, n, n);
    
    Vcheck = zeros(n, n, n);
    V = zeros(n, n, n);
    Qcheck = zeros(n, n, n);
    Q = zeros(n, n, n);
    
    
else
    load('HY')
end




for k = 1:Nl
    lambda = l(k);
    
    g2 = 1/10*lambda*(2-lambda)/(1+lambda)^2;
    g3 = 1/500*(16-61*lambda + 180*lambda^2 + 2*lambda^3)/(1+lambda)^2;
    g5 = 1/20*(lambda^2 + 7*lambda^3) / (1+lambda)^2;
    g6 = 1/1000*(43*lambda + 147*lambda^2 - 185*lambda^3 +221*lambda^4)/(1+lambda)^2;
    
    
    invlambda = 1.0/lambda;
    invg2 = 1/10*invlambda*(2-invlambda)/(1+invlambda)^2;
    invg3 = 1/500*(16-61*invlambda + 180*invlambda^2 + 2*invlambda^3)/(1+invlambda)^2;
    invg5 = 1/20*invlambda^2*(1 + 7*invlambda) / (1+invlambda)^2;
    invg6 = invlambda/1000*(43 + 147*invlambda - 185*invlambda^2 +221*invlambda^3)/(1+invlambda)^2;
    
    H_11Y = -g3;
    H_22Y = -invg3;
    
    for n = 2:2:m
        n
        fm  = 0;
        fm2 = 0;
        
        for q = 0:n
            [temp, Vcheck, V, Pcheck, P, Qcheck, Q] = P_HY(2, n-q, q, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            tempf=(1+lambda)^(q-n);
            invtempf=(1+invlambda)^(q-n);
            fm = fm -3/8*tempf*temp*(lambda/(1+lambda))^q;
            fm2 = fm2 -3/8*invtempf*temp*(invlambda/(1+invlambda))^q;
        end
        
        temp = fm - 2/n*g2+4/n*g3/(n+2);
        H_11Y = H_11Y + temp;
        
        temp = fm2 - 2/n*invg2+4/n*invg3/(n+2);
        H_22Y = H_22Y + temp;
        
    end
    
    H_12Y =  g5*log(4) -2*g6;
    H_21Y =  invg5*log(4) -2*invg6;
    for n = 1:2:m
        n
        fm = 0;
        fm2 = 0;
  
        for q = 0:n
            [temp, Vcheck, V, Pcheck, P, Qcheck, Q] = P_HY(2, n-q, q, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            tempf=(1+lambda).^(q+1-n);
            invtempf=(1+invlambda)^(q+1-n);
            fm = fm -3/8*tempf*temp*(lambda/(1+lambda))^(q+1);
            fm2 = fm2-3/8*invtempf*temp*(invlambda/(1+invlambda))^(q+1);
        end
        
        temp = fm -2/n*g5+4/n*g6/(n+2);
        H_12Y = H_12Y + temp;
        temp = fm2 -2/n*invg5+4/n*invg6/(n+2);
        H_21Y = H_21Y + temp;
    end
    H_12Y = H_12Y*8.0/(1.0+lambda)^3;
    H_21Y = H_21Y*8.0/(1.0+invlambda)^3;
    H(k, :) = [H_11Y, H_12Y, H_21Y, H_22Y];
    
end
% if i == 1
save('HY','P','Pcheck','V','Vcheck', 'Q', 'Qcheck');
% end

csvwrite('HY.dat',[l',H])

% end

function [val, Vcheck, V, Pcheck, P, Qcheck, Q] = P_HY(n, p, q, val, Vcheck, V, Pcheck, P, Qcheck, Q)
if Pcheck(n+200, p+200, q+200) == 1
    val = P(n+200, p+200, q+200);
else
    if p == 0 && q == 0
            val = val;
    else
        for s = 1:q
            [temp1, Vcheck, V, Pcheck, P, Qcheck, Q] = P_HY(s, q-s, p-n+1, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp2, Vcheck, V, Pcheck, P, Qcheck, Q] = P_HY(s, q-s, p-n-1, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp3, Vcheck, V, Pcheck, P, Qcheck, Q] = V_HY(s, q-s-2, p-n+1, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp4, Vcheck, V, Pcheck, P, Qcheck, Q] = Q_HY(s, q-s-1, p-n+1, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            
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


function [val, Vcheck, V, Pcheck, P, Qcheck, Q] = V_HY(n, p, q, val, Vcheck, V, Pcheck, P, Qcheck, Q)
if Vcheck(n+200, p+200, q+200) == 1
    val = V(n+200, p+200, q+200);
else
    if p == 0 && q == 0
        val = val;
    else
        for s = 1:q
            [temp1, Vcheck, V, Pcheck, P, Qcheck, Q] = P_HY(s, q-s, p-n-1, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            val = val + nchoosek(n+s, n+1)*temp1;
        end
        val = val*2*n/((n+1)*(2*n+3));
        [temp1, Vcheck, V, Pcheck, P, Qcheck, Q] = P_HY(n, p, q, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
        val = val + temp1;
    end
    Vcheck(n+200, p+200, q+200) = 1;
    V(n+200, p+200, q+200) = val;
end
end


function [val, Vcheck, V, Pcheck, P, Qcheck, Q] = Q_HY(n, p, q, val, Vcheck, V, Pcheck, P, Qcheck, Q)
if Qcheck(n+200, p+200, q+200) == 1
    val = Q(n+200, p+200, q+200);
else
    if p == 0 && q == 0
        if n == 1
            val = val+1;
        end
    else
        for s = 1:q
            [temp1, Vcheck, V, Pcheck, P, Qcheck, Q] = Q_HY(s, q-s-1, p-n, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp2, Vcheck, V, Pcheck, P, Qcheck, Q] = P_HY(s, q-s, p-n, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            
            val = val + nchoosek(n+s, n+1)*...
                (s/(n+1)*temp1 ...
                -3/(2*n*s*(n+1))*temp2);
        end
    end
    Qcheck(n+200, p+200, q+200) = 1;
    Q(n+200, p+200, q+200) = val;
end
end
