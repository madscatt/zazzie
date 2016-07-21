function [new]=recalc_D(A,B,C,D,E,dir_cos)

X=[A;B;C;D;E];

new=dir_cos*X;

return