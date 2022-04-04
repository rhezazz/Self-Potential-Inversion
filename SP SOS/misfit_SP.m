%Fungsi Objektif SP
function [misfit] = misfit_SP(V_obs,V_cal)
ls = length(V_cal);
for i = 1 : ls
    m(i) = abs(V_obs(i) - V_cal(i));
    n(i) = abs(V_obs(i) + V_cal(i));
end
misfit = (2.*sum(m))./(sum(m) + sum(n));
end