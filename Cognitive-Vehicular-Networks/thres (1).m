clc;
clear all;
close all;
pf=0.1;
N=1000;
k=0;
for i=1:length(weight)
    k=k+y(i)*((I(i)*SNR)+1);
%s=E^2;
n=noise^2;
%thres=(qfuninv(pf))*(sqrt((2*N)*(((y(1)^2)*(n^2))+(k^2))))+(N*(y(1)*(n^2)+k));
%sq_y0=y0^2;
%l=sq_y0*(n+sq)^2;
%thre=(qfuncinv(pd))*(sqrt(((2*N)*(l^2))+(sum^2)))+(N*(y0*n)+sum);

thres=(qfuninv(pf)*n)*(sqrt((2*N)*((y(1)^2)+(k^2))))+((n*N(y(1)+k)));

