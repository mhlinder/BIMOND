%% BIMOND3_example.m
% M. Henry Linder (mhlinder@gmail.com)
% 
% This M-file reproduces the example given in Appendix B of "An
% Algorithm For Monotone Piecewise Bicubic Interpolation" (Carlson and
% Fritsch, 1989).

clear
close all
addpath 'slatec_pchic/';

%% Input data
x = 1:4;
y = 1:4;
p = [0, 2.999, 3, 8; ...
     2, 3, 9, 10; ...
     19.998, 19.999, 20, 20.001; ...
     19.999, 20, 20.001, 20.002];

pp2d = BIMOND3(x, y, p);
% fnplt(pp2d);
