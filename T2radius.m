function [rt1, rt2] = T2radius(Simdata, tstep, K, h)

% K = 1 / (4/3 * pi * 1e-9);

c = max(Simdata.c(tstep,:), 0);
r = max(Simdata.r(tstep,:), 0);
p = max(Simdata.p(tstep,:), 0);


T = (c+r)/K;
% h = .1;

Area = length(find(T >= .16)) * h^2;
rt2 = sqrt(Area/pi);

Area = length(find(T >= .8)) * h^2;
rt1 = sqrt(Area/pi);

