function [cube, rawpw] = isac_gen_cube(Xd, tg, o)
% Mo hinh monostatic: Y[k,m,a] = X[k,m]*sum_p alpha_p e^{-j2pi k DF tau_p} e^{j2pi m Ts fd_p} a_a(theta_p) + SI + nhieu
% Sau do: H = Y./X (loc nguoc/ZF) -> cua so Hann -> MOT lan 2D-FFT moi anten -> cube [range x doppler x anten] (single)
%  tg: [R(m) v(m/s) theta(deg) SNR_RE(dB)];  o: DF, Tsym, NR, NRK, Nant, clutter, mti, siDb, analogDb, sic, nlDbc, adc, seed
c = 299792458; lambda = c/3.5e9; [K,M] = size(Xd); rng(o.seed);
P = [tg, 2*pi*rand(size(tg,1),1)]; C = zeros(0,5);
if o.clutter
    C1 = [15+365*rand(40,1), zeros(40,1), -70+140*rand(40,1), -33+11*rand(40,1)];
    C2 = [30+320*rand(25,1), 0.35*randn(25,1), -70+140*rand(25,1), -40+8*rand(25,1)];
    C = [C1; C2]; C(:,5) = 2*pi*rand(size(C,1),1);
end
Aall = [P; C]; R = Aall(:,1); v = Aall(:,2); th = Aall(:,3); amp = 10.^(Aall(:,4)/20); ph = Aall(:,5);
tau = 2*R/c; fd = 2*v/lambda; kk = (0:K-1).' - K/2;
Ek = exp(-1j*2*pi*(kk*o.DF)*tau.');                 % K x P
Em = exp(1j*2*pi*fd*((0:M-1)*o.Tsym));              % P x M
A  = isac_steer(th, o.Nant);                         % Nant x P
hn = @(n) 0.5*(1-cos(2*pi*(0:n-1).'/(n-1)));
win = hn(K)*hn(M).';
cube = zeros(o.NRK, M, o.Nant, 'single'); rawpw = zeros(o.NRK, M, 'single');
hasSI = ~isempty(o.siDb);
if hasSI
    hsi = (randn(o.Nant,3)+1j*randn(o.Nant,3))/sqrt(2);
    siTau = [0.5 12 30]*1e-9; siW = [1 0.2 0.05]/norm([1 0.2 0.05]);
    siAmp = 10^((o.siDb - o.analogDb)/20);
end
cn = @() (randn(K,M)+1j*randn(K,M))/sqrt(2);
for a = 1:o.Nant
    Hm = (Ek .* (amp.'.*A(a,:).*exp(1j*ph.'))) * Em;
    Y = Xd.*Hm + cn();                               % thu = X*H + nhieu nhiet (1 don vi / RE)
    if hasSI
        hk = exp(-1j*2*pi*(kk*o.DF)*siTau)*(hsi(a,:).'.*siW.');
        Y = Y + siAmp*hk.*Xd;                        % nhieu tu sinh (SI) tuyen tinh
        Y = Y + siAmp*10^(o.nlDbc/20)*cn();   % meo phi tuyen PA (khong khu duoc bang SIC tuyen tinh)
        if o.adc, Y = Y + siAmp/10^(61.9/20)*cn(); end                          % luong tu hoa ADC (SQNR ~62 dB)
    end
    Hh = Y./Xd;                                      % loc nguoc (ZF)
    if o.sic, Hh = Hh - mean(Hh(:,1:16),2); end      % SIC so: LS tren 16 symbol hieu chuan
    Rr = ifft(Hh.*win, o.NR, 1)*o.NR; Rr = Rr(1:o.NRK,:);
    rawpw = rawpw + single(abs(fftshift(fft(Rr,[],2),2)).^2);
    if o.mti, Hh = Hh - mean(Hh,2); end              % khu thanh phan tinh (slow-time)
    Rr = ifft(Hh.*win, o.NR, 1)*o.NR; Rr = Rr(1:o.NRK,:);
    cube(:,:,a) = single(fftshift(fft(Rr,[],2),2));  % 2D-FFT: range (IFFT theo k) x Doppler (FFT theo m)
end
end
