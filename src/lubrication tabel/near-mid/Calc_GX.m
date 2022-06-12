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

%The output is of the form XG_11, XG_12, XG_22, XG_21

%This code relies on a recursive formula for finding the coefficients. The
%output from a perious run can be saved in a matrix GX.dat to be reused on
%future runs, greatly cutting the computational time. 
%%
clear;clc;
format long

% function [] =   Calc_GX()
warning off;
l = [1.0,0.5,0.25,0.2,0.125,0.1,0.01,2,3,4,5,10,100]; %Lambda

m = 180; %Number of terms to include in the recursive sum
Nl = length(l);
A = zeros(Nl, 4);

i = 0; %Set i to 0 if you wish to load in a matrix GX.dat from a previous run. 
if i == 1
    n = 500;
    Pcheck = zeros(n, n, n);
    P = zeros(n, n, n);
    
    Vcheck = zeros(n, n, n);
    V = zeros(n, n, n);
else
    load('GX')
end



for k = 1:Nl
    lambda = l(k);
    
    
    invlambda = 1/lambda;
    
    g1 = 3*lambda^2/(1+lambda)^3;
    g2 = 3/10*lambda*(1+12*lambda-4*lambda^2)/(1+lambda)^3;
    g3 = 1/140*(5+181*lambda-453*lambda^2+566*lambda^3-65*lambda^4)/(1+lambda)^3;
    G_11X = 1/4*g1 + g2*log(4) - 2*g3;
    
    invg1 = 3*invlambda^2/(1+invlambda)^3;
    invg2 = 3/10*invlambda*(1+12*invlambda-4*invlambda^2)/(1+invlambda)^3;
    invg3 = 1/140*(5+181*invlambda-453*invlambda^2+566*invlambda^3-65*invlambda^4)/(1+invlambda)^3;
    
    G_22X = 1/4*invg1 + invg2*log(4) - 2*invg3;
    
    
    for n = 1:2:m
        n
        fm  = 0;
        fm2 = 0;
        for q = 0:n
            [temp, Vcheck, V, Pcheck, P]  = P_A(2, n-q, q, 0,  Vcheck, V, Pcheck, P);
            tempf=(1+lambda)^(q-n);
            invtempf=(1+invlambda)^(q-n);
            fm = fm +3/4*tempf*temp*(lambda/(1+lambda))^q;
            fm2= fm2 +3/4*tempf*temp*(invlambda/(1+invlambda))^q;
        end
        %f = fm*2^(n);
        temp = fm - g1-2/n*g2+4*g3/(n*(n+2));
        G_11X = G_11X + temp;
        temp = fm2 - invg1-2/n*invg2+4*invg3/(n*(n+2));
        G_22X = G_22X + temp;
    end
    
    G_12X = 1/4*g1 + g3;
    G_21X = 1/4*invg1 + invg3;
    for n = 2:2:m
        n
        fm = 0;
        fm2 = 0;
        for q = 0:n
            [temp, Vcheck, V, Pcheck, P]  = P_A(2, n-q, q, 0,  Vcheck, V, Pcheck, P);
            tempf=(1+lambda)^(q-n);
            invtempf=(1+invlambda)^(q-n);
            fm = fm +3/4*tempf*temp*(lambda/(1+lambda))^q;
            fm2= fm2 +3/4*tempf*temp*(invlambda/(1+invlambda))^q;
        end
        
        temp =fm - g1-2/n*g2+4/n*g3/(n+2);
        G_12X = G_12X - temp;
        temp =fm2 - invg1-2/n*invg2+4/n*invg3/(n+2);
        G_21X = G_21X - temp;
    end
    G_12X = G_12X*4/(1+lambda)^2;
    G_21X = G_21X*4/(1+invlambda)^2;
    G(k, :) = [G_11X, G_12X, G_21X, G_22X];
    
end
save('GX','P','Pcheck','V','Vcheck');

csvwrite('GX.dat',[l',G])



% end


function [val, Vcheck, V, Pcheck, P] = P_A(n, p, q, val, Vcheck, V, Pcheck, P)
if Pcheck(n+200, p+200, q+200) == 1
    val = P(n+200, p+200, q+200);
else
    if p == 0 && q == 0
        if n == 1
            val = val + 1;
        end
    else
        for s = 1:q
            [temp1, Vcheck, V, Pcheck, P]  = P_A(s, q-s, p-n+1, 0,  Vcheck, V, Pcheck, P);
            [temp2, Vcheck, V, Pcheck, P]  =P_A(s, q-s, p-n-1, 0,  Vcheck, V, Pcheck, P);
            [temp3, Vcheck, V, Pcheck, P]  = V_A(s, q-s-2, p-n+1, 0,  Vcheck, V, Pcheck, P);
            
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
            [temp1, Vcheck, V, Pcheck, P]  =P_A(s, q-s, p-n-1, 0,  Vcheck, V, Pcheck, P);
            val = val + nchoosek(n+s, n)*temp1;
        end
        val = -val*2*n/((n+1)*(2*n+3));
        [temp1, Vcheck, V, Pcheck, P]  = P_A(n, p, q, 0,  Vcheck, V, Pcheck, P);
        val = val + temp1;
    end
    Vcheck(n+200, p+200, q+200) = 1;
    V(n+200, p+200, q+200) = val;
end
end

















