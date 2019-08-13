function r = radius_variable_02(xy)

r1 = 0.01;  r2 = 0.05;  % specify the two constant r-values and the (bengt had 0.01 and 0.05)
d1 = 0.1;  d2 = 0.3;  % points of transition between them (bengt had 0.1 and 0.3, I had 1.5 and 1.8).

d = abs(xy(:,1)-xy(:,2))/sqrt(2);
r = 0.5*((r2-r1)/(d2-d1)*(abs(d-d1)-abs(d-d2))+(r2+r1));