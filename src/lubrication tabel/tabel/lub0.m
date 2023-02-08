clear;clc;
load FarfieldScalar.mat
xi=Lb0_list;
invxi=1.0./xi;
loginvxi=log(invxi);
xiloginvxi=xi.*loginvxi;
x11a=-1.23041+0.25*invxi+1.8918*xi+9.0*loginvxi/40.0+3.0*xiloginvxi/112.0;
x12a=-x11a+0.00312-0.0011*xi;
xm=invxi/6.0+3.0*loginvxi/20.0+47.0*xiloginvxi/280.0 ...
         -0.740815+0.706802*xi;