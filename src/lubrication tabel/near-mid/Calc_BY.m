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

%The output is of the form XB_11, XB_12, XB_22, XB_21

%This code relies on a recursive formula for finding the coefficients. The
%output from a perious run can be saved in a matrix BY.dat to be reused on
%future runs, greatly cutting the computational time.
%%
clear;clc;

% function [] =   Calc_BY()
l = [1.0,0.5,0.25,0.2,0.125,0.1,0.01,2,3,4,5,10,100]; %Lambda
m = 150; %Number of terms to include in the recursive sum

Nl = length(l);
B = zeros(Nl, 4);
warning off;
i = 0; %Set i to 0 if you wish to load in a matrix BY.dat from a previous run. 
    if i == 1
        n = 500;
        Pcheck = zeros(n, n, n);
        P = zeros(n, n, n);

        Vcheck = zeros(n, n, n);
        V = zeros(n, n, n);

        Qcheck = zeros(n, n, n);
        Q = zeros(n, n, n);
    else
        load('BY.mat')
    end

for k = 1:Nl
    lambda = l(k);
    g2 = -1/5*lambda*(4+lambda)/(1+lambda)^2;
    g3 = -1/250*(32-33*lambda+83*lambda^2+43*lambda^3)/(1+lambda)^2;
    B_11Y = 0;
    for n = 1:2:m
        n
        fm = 0;
        tempf=(1+lambda)^(-n);
        for q = 0:n
            [temp, Vcheck, V, Pcheck, P, Qcheck, Q]  = Q_BY(1, n-q, q, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            fm = fm +tempf*temp*lambda^q;
        end
        f = 2*fm;
        temp = f - 2/n*g2+4/n*g3/(n+2);
        B_11Y = B_11Y + temp;
    end
    
    B_11Y = B_11Y + 2*g2*log(2)-2*g3;
    B_12Y = 0;
    for n = 2:2:m
        n
        fm = 0;
        tempf=(1+lambda)^(-n);
        for q = 0:n
            [temp, Vcheck, V, Pcheck, P, Qcheck, Q]  = Q_BY(1, n-q, q, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            fm = fm + tempf*temp*lambda^q;
        end
        f = 2*fm;
        
        temp = f - 2/n*g2+4/n*g3/(n+2);
        B_12Y = B_12Y + temp;
        
    end
    B_12Y = (B_12Y-g3)/(-.25*(1+lambda)^2);
    
    lambda = 1/lambda;
    g2 = -1/5*lambda*(4+lambda)/(1+lambda)^2;
    g3 = -1/250*(32-33*lambda+83*lambda^2+43*lambda^3)/(1+lambda)^2;
    
    B_22Y = 2*g2*log(2)-2*g3;
    for n = 1:2:m
        
        fm  = 0;
        tempf=(1+lambda)^(-n);
        for q = 0:n
            [temp, Vcheck, V, Pcheck, P, Qcheck, Q]  = Q_BY(1, n-q, q, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            fm = fm + tempf*temp*lambda^q;
        end
        f = 2*fm;
        temp = f - 2/n*g2+4/n*g3/(n+2);
        B_22Y = B_22Y + temp;
    end
    
    
    
    B_21Y = 0;
    for n = 2:2:m
        fm = 0;
        tempf=(1+lambda)^(-n);
        for q = 0:n
            [temp, Vcheck, V, Pcheck, P, Qcheck, Q]  = Q_BY(1, n-q, q, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            fm = fm + tempf*temp*lambda^q;
        end
        f = 2*fm;
        
        temp = f - 2/n*g2+4/n*g3/(n+2);
        B_21Y = B_21Y + temp;
        
    end
    B_21Y = (B_21Y-g3)/(-.25*(1+lambda)^2);
    
    
    B(k, :) =[B_11Y, B_12Y, B_21Y,B_22Y];
    
end

save('BY','P','Pcheck','V','Vcheck', 'Q', 'Qcheck');

csvwrite('BY.dat',[l',B])

% end

function [val, Vcheck, V, Pcheck, P, Qcheck, Q] = P_BY(n, p, q, val, Vcheck, V, Pcheck, P, Qcheck, Q)
if Pcheck(n+200, p+200, q+200) == 1
    val = P(n+200, p+200, q+200);
else
    if p == 0 && q == 0
        if n == 1
            val = val + 1;
        end
    else
        for s = 1:q
            [temp1, Vcheck, V, Pcheck, P, Qcheck, Q] = P_BY(s, q-s, p-n+1, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp2, Vcheck, V, Pcheck, P, Qcheck, Q] = P_BY(s, q-s, p-n-1, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp3, Vcheck, V, Pcheck, P, Qcheck, Q] = V_BY(s, q-s-2, p-n+1, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp4, Vcheck, V, Pcheck, P, Qcheck, Q] = Q_BY(s, q-s-1, p-n+1, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
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


function [val,  Vcheck, V, Pcheck, P, Qcheck, Q] = V_BY(n, p, q, val,  Vcheck, V, Pcheck, P, Qcheck, Q)

if Vcheck(n+200, p+200, q+200) == 1
    val = V(n+200, p+200, q+200);
else
    if p == 0 && q == 0
        if n == 1
            val = val + 1;
        end
    else
        for s = 1:q
            [temp1, Vcheck, V, Pcheck, P, Qcheck, Q] = P_BY(s, q-s, p-n-1, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            val = val + nchoosek(n+s, n+1)*temp1;
        end
        val = val*2*n/((n+1)*(2*n+3));
        [temp1, Vcheck, V, Pcheck, P, Qcheck, Q] = P_BY(n, p, q, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
        val = val + temp1;
    end
    Vcheck(n+200, p+200, q+200) = 1;
    V(n+200, p+200, q+200) = val;
end
end


function [val,  Vcheck, V, Pcheck, P, Qcheck, Q] = Q_BY(n, p, q, val,  Vcheck, V, Pcheck, P, Qcheck, Q)
if Qcheck(n+200, p+200, q+200) == 1
    val = Q(n+200, p+200, q+200);
else
    if p == 0 && q == 0
        val = val;
    else
        for s = 1:q
            [temp1, Vcheck, V, Pcheck, P, Qcheck, Q] = Q_BY(s, q-s-1, p-n, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            [temp2, Vcheck, V, Pcheck, P, Qcheck, Q] = P_BY(s, q-s, p-n, 0, Vcheck, V, Pcheck, P, Qcheck, Q);
            
            val = val + nchoosek(n+s, n+1)*...
                (s/(n+1)*temp1 ...
                -3/(2*n*s*(n+1))*temp2);
        end
    end
    Qcheck(n+200, p+200, q+200) = 1;
    Q(n+200, p+200, q+200) = val;
end

end









