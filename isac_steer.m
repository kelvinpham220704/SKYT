function A = isac_steer(thDeg, N)
% Steering vector ULA lambda/2: N x numel(thDeg)
A = exp(1j*pi*(0:N-1).'*sind(thDeg(:).'));
end
