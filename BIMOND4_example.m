%% BIMOND4_example.m
% M. Henry Linder (mhlinder@gmail.com)
% 
% This M-file reproduces the example given in Section 5 of "A
% Bivariate Interpolation Algorithm For Data That Are Monotone In One
% Variable" (Carlson and Fritsch, 1991).

clear
close all
addpath 'slatec_pchic/';

%% Input data
x = 1:4;
y = 1:4;
p = [0, 2, 19.998, 19.997; ...
     2.999, 3, 19.999, 19.998; ...
     3, 9, 20, 19.999; ...
     8, 10, 20.001, 20];

pp2d = BIMOND4(x, y, p);
% fnplt(pp2d);