function [V] = fwd_SP(k,x0,z,alpha,a,x)
for i = 1 : length(x)
    V(i) = k*(log10((((x(i)-x0)-a*cosd(alpha)).^2+(z - a*sind(alpha)).^2)/((((x(i)-x0)+a*cosd(alpha)).^2+(z + a*sind(alpha)).^2))));
end
end