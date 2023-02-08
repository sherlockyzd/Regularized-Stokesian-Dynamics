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

%The output is of the form YM_11, YM_12, YM_22, YM_21

%This code relies on a recursive formula for finding the coefficients. The
%output from a perious run can be saved in a matrix MY.dat to be reused on
%future runs, greatly cutting the computational time.
%%
clear;clc;
format long

% function [] =     Calc_MY()
warning off;
l = [1.0,0.5,0.25,0.2,0.125,0.1,0.01,2,3,4,5,10,20,100]; %Lambda
m = 180; %Number of terms to include in the recursive sum
Nl = length(l);
A = zeros(Nl, 4);

i = 0; %Set i to 0 if you wish to load in a matrix MY.dat from a previous run. 
if i == 1
    n = 500;
    Pcheck = zeros(n, n, n);
    P = zeros(n, n, n);
    
    Qcheck = zeros(n, n, n);
    Q = zeros(n, n, n);
    
    
    Vcheck = zeros(n, n, n);
    V = zeros(n, n, n);
else
    load('MY')
end



for k = 1:Nl
    lambda = l(k);
    g2 = 6/25*lambda*(1-lambda+4*lambda^2)/(1+lambda)^3;
    g3 = 1/625*(24 - 201*lambda+882*lambda^2-1182*lambda^3+591*lambda^4)/(1+lambda)^3;
    g5 = 3/50*lambda^2*(7 -10*lambda +7*lambda^2)/(1+lambda)^3;
    g6 = 3/2500*lambda*(221-728*lambda+1902*lambda^2-728*lambda^3+221*lambda^4)/(1+lambda)^3;
    
    invlambda = 1/lambda;
    invg2 = 6/25*invlambda*(1-invlambda+4*invlambda^2)/(1+invlambda)^3;
    invg3 = 1/625*(24 - 201*invlambda+882*invlambda^2-1182*invlambda^3+591*invlambda^4)/(1+invlambda)^3;
    invg5 = 3/50*invlambda^2*(7 -10*invlambda +7*invlambda^2)/(1+invlambda)^3;
    invg6 = 3/2500*invlambda*(221-728*invlambda+1902*invlambda^2-728*invlambda^3+221*invlambda^4)/(1+invlambda)^3;
    
    M_11Y =  - g3 + 1;
    M_22Y = - invg3 + 1;
    for n = 2:2:m
        n
        fm  = 0;
        fm2 = 0;
        for q = 0:n
            [temp, Vcheck, V, Pcheck, P, Qcheck, Q] = P_MY(2, n-q, q, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            tempf=(1+lambda)^(q-n);
            invtempf=(1+invlambda)^(q-n);
            fm = fm + tempf*temp*(lambda/(1+lambda))^q;
            fm2 = fm2 + invtempf*temp*(invlambda/(1+invlambda))^q;
        end
        temp = fm - 2*g2/n + 4*g3/(n*(n+2));
        M_11Y = M_11Y + temp;
        temp = fm2- 2*invg2/n + 4*invg3/(n*(n+2));
        M_22Y = M_22Y + temp;
    end
    
    M_12Y =  g5*log(4) - 2*g6;
    M_21Y =  invg5*log(4) - 2*invg6;
    
    for n = 1:2:m
        n
        fm = 0;
        fm2 = 0;
        for q = 0:n
            [temp, Vcheck, V, Pcheck, P, Qcheck, Q] = P_MY(2, n-q, q, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            tempf=(1+lambda)^(q+1-n);
            invtempf=(1+invlambda)^(q+1-n);
            fm = fm + tempf*temp*(lambda/(1+lambda))^(q+1);
            fm2 = fm2 + invtempf*temp*(invlambda/(1+invlambda))^(q+1);
        end
        temp = fm - 2*g5/n + 4*g6/(n*(n+2));
        M_12Y = M_12Y + temp;
        temp = fm2- 2*invg5/n + 4*invg6/(n*(n+2));
        M_21Y = M_21Y + temp;
        
    end
    M_12Y = M_12Y*8/(1+lambda)^3;
    M_21Y = M_21Y*8/(1+invlambda)^3;
    M(k, :) = [M_11Y, M_12Y, M_21Y, M_22Y];
end

save('MY','P','Pcheck','V','Vcheck','Q','Qcheck')

csvwrite('MY.dat',[l',M])

% end

function [val, Vcheck, V, Pcheck, P, Qcheck, Q] = P_MY(n, p, q, val, Vcheck, V, Pcheck, P, Qcheck, Q)
if Pcheck(n+200, p+200, q+200) == 1
    val = P(n+200, p+200, q+200);
else
    
    if p == 0 && q == 0
        if n == 2
            val = val + 1;
        end
    else
        for s = 1:q
            [temp1, Vcheck, V, Pcheck, P, Qcheck, Q] = P_MY(s, q-s, p-n+1, 0,  Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp2, Vcheck, V, Pcheck, P, Qcheck, Q] = P_MY(s, q-s, p-n-1, 0,  Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp3, Vcheck, V, Pcheck, P, Qcheck, Q] = V_MY(s, q-s-2, p-n+1, 0,  Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp4, Vcheck, V, Pcheck, P, Qcheck, Q] = Q_MY(s, q-s-1, p-n+1, 0,  Vcheck, V, Pcheck, P, Qcheck, Q);
            val = val+ 1/(n+1)*nchoosek(n+s, n+1)*(...
                (n+0.5)*((n*s+4)*(n+s)-2*(n*s+1)^2)/(s*(2*s-1)*(n+s))*temp1...
                +n*(n-0.5)*temp2...
                +n*(4*n^2-1)/(2*(2*s+1))*temp3...
                -(4*n*n-1)*temp4);
        end
    end
    Pcheck(n+200, p+200, q+200) = 1;
    P(n+200, p+200, q+200) = val;
end
end


function [val, Vcheck, V, Pcheck, P, Qcheck, Q] = V_MY(n, p, q, val,  Vcheck, V, Pcheck, P, Qcheck, Q)
if Vcheck(n+200, p+200, q+200) == 1
    val = V(n+200, p+200, q+200);
else
    
    if p == 0 && q == 0
        if n == 2
            val = val + 1;
        end
    else
        for s = 1:q
            [temp, Vcheck, V, Pcheck, P, Qcheck, Q] = P_MY(s, q-s, p-n-1, 0,  Vcheck, V, Pcheck, P, Qcheck, Q);
            val = val + 1/(n+1)*nchoosek(n+s, n+1)*temp;
        end
        val = val*2*n/(2*n+3);
        [temp, Vcheck, V, Pcheck, P, Qcheck, Q] = P_MY(n, p, q, 0,  Vcheck, V, Pcheck, P, Qcheck, Q);
        val = val + temp;
    end
    Vcheck(n+200, p+200, q+200) = 1;
    V(n+200, p+200, q+200) = val;
end
end

function [val, Vcheck, V, Pcheck, P, Qcheck, Q] = Q_MY(n, p, q, val,  Vcheck, V, Pcheck, P, Qcheck, Q)
if Qcheck(n+200, p+200, q+200) == 1
    val = Q(n+200, p+200, q+200);
else
    if p == 0 && q == 0
        val = val;
    else
        for s = 1:q
            [temp, Vcheck, V, Pcheck, P, Qcheck, Q] = Q_MY(s, q-s-1, p-n, 0,  Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp1, Vcheck, V, Pcheck, P, Qcheck, Q] = P_MY(s, q-s, p-n, 0,  Vcheck, V, Pcheck, P, Qcheck, Q);
            val = val + 1/(n+1)*nchoosek(n+s, n+1)*( s*temp...
                - 1/(n*s)*temp1);
        end
    end
    Qcheck(n+200, p+200, q+200) = 1;
    Q(n+200, p+200, q+200) = val;
end
end


