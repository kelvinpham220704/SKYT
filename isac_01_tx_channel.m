%Tx ISAC-OFDM 
clear; clc; close all; isac_params; rng(7);
nSlots = 2; nTx = 1; pdsch = isac_pdsch(1, NRB);
grid = []; typ = [];
for s = 1:nSlots
    carrier.NSlot = s-1;
    [g, pInd, ~, dInd] = isac_tx_slot(carrier, pdsch, nTx);
    t = zeros(size(g)); t(pInd) = 1; t(dInd) = 2;          % 1 = PDSCH data, 2 = DM-RS
    grid = [grid g]; typ = [typ t];                         %#ok<AGROW>
end
carrier.NSlot = 0;
wave = nrOFDMModulate(carrier, grid);  wave = wave/sqrt(mean(abs(wave).^2));
nS = size(wave,1); symLen = ofdmInfo.SymbolLengths(1:14*nSlots);
symStart = [0 cumsum(symLen(1:end-1))];
%resource grid
figure('Name','Resource grid','Position',[100 100 1100 420]);
subplot(1,2,1); imagesc(0:13, 0:NRB-1, typ(1:12:end,1:14)); axis xy; colormap(gca,[0.95 0.72 0.71; 0.12 0.23 0.58; 0.12 0.23 0.58]);
xlabel('OFDM symbol (1 slot)'); ylabel('Resource Block'); title('ISAC-OFDM resource grid (100 MHz, SCS 30 kHz)');
subplot(1,2,2); imagesc(0:13, 0:23, 10*log10(abs(grid(1:24,1:14)).^2)); axis xy; colorbar;
xlabel('OFDM symbol'); ylabel('Sub-carrier'); title('Zoom 2 RB: |X|^2 [dB] (DM-RS +3 dB neu co boost)');
saveas(gcf, fullfile(outDir,'tx_grid.png'));
%AWGN
snr = 10; sig = mean(abs(wave).^2);
n = sqrt(sig/10^(snr/10)/2)*(randn(nS,1)+1j*randn(nS,1)); xn = wave + n;
Yc = nrOFDMDemodulate(carrier, wave); Yn = nrOFDMDemodulate(carrier, xn);
msk = (typ==1); a0 = sum(conj(Yc(msk)).*grid(msk))/sum(abs(Yc(msk)).^2); Yn = a0*Yn;
evm = sqrt(mean(abs(Yn(msk)-grid(msk)).^2)/mean(abs(grid(msk)).^2));
figure('Name','AWGN','Position',[100 100 1200 380]);
k = 1:2*(Nfft+288); subplot(1,3,1); plot(k/fs*1e6, abs(wave(k)),'b', k/fs*1e6, abs(xn(k)),'r'); legend('Clean','+AWGN'); xlabel('us'); title('Mien thoi gian');
[f,p1] = welchPSD(wave,fs); [~,p2] = welchPSD(xn,fs); subplot(1,3,2); plot(f/1e6,10*log10(p1),'b',f/1e6,10*log10(p2),'r'); xlabel('MHz'); ylabel('dB/Hz'); title('Pho'); legend('Clean','+AWGN');
idx = find(msk); idx = idx(randperm(numel(idx),4000)); subplot(1,3,3); plot(real(Yn(idx)),imag(Yn(idx)),'r.'); axis equal; title(sprintf('QPSK sau FFT, EVM = %.1f%%',100*evm));
saveas(gcf, fullfile(outDir,'tx_awgn.png'));
%4 kenh TDL (A,B,C,D)
DS = struct('A',300e-9,'B',300e-9,'C',300e-9,'D',100e-9);
fgrid = linspace(-50e6,50e6,200); rxTDL = struct();
for m = 'ABCD'
    tdl = nrTDLChannel; tdl.DelayProfile = ['TDL-' m]; tdl.DelaySpread = DS.(m);
    tdl.MaximumDopplerShift = fD; tdl.SampleRate = fs;
    tdl.NumTransmitAntennas = 1; tdl.NumReceiveAntennas = 1; tdl.Seed = double(m);
    ti = info(tdl); pad = ti.MaximumChannelDelay;
    [rx, pg, ~] = tdl([xn; zeros(pad,1)]);               % pg: [Nsamples x Npaths x Ntx x Nrx]
    rxTDL.(m) = rx(1:nS);
    tau = ti.PathDelays(:).'; pw = ti.AveragePathGains(:).';
    ii = 1:round(fs*1e-3/80):nS; g = squeeze(pg(ii,:,1,1));
    Hft = g*exp(-1j*2*pi*tau.'*fgrid);
    figure('Name',['TDL-' m],'Position',[100 100 1300 360]);
    subplot(1,3,1); stem(tau*1e9, pw, 'r','filled'); ylim([-40 3]); xlabel('Tre [ns]'); ylabel('Cong suat [dB]'); title(sprintf('TDL-%s: PDP (%d tap)',m,numel(tau)));
    subplot(1,3,2); imagesc(ii/fs*1e3, fgrid/1e6, 20*log10(abs(Hft).'+1e-6), [-25 8]); axis xy; colorbar; xlabel('ms'); ylabel('MHz'); title(sprintf('|H(f,t)| [dB], fD = %.0f Hz',fD));
    kz = 3*(Nfft+288)+288+(1:160); subplot(1,3,3); plot(kz/fs*1e6,abs(xn(kz)),'b',kz/fs*1e6,abs(rx(kz)),'r'); legend('Tx+AWGN',['Sau TDL-' m]); xlabel('us'); title('Da duong + fading (zoom)');
    saveas(gcf, fullfile(outDir,['tx_tdl_' m '.png']));
end

%nhieu phi tuyen phan cung (TDL-D, SNR 25 dB)
tdl = nrTDLChannel; tdl.DelayProfile = 'TDL-D'; tdl.DelaySpread = 100e-9; tdl.MaximumDopplerShift = fD; tdl.SampleRate = fs;
tdl.NumTransmitAntennas = 1; tdl.NumReceiveAntennas = 1; tdl.Seed = 99;
ti = info(tdl); [rx, pg] = tdl([wave; zeros(ti.MaximumChannelDelay,1)]); y = rx(1:nS);
y = y + sqrt(mean(abs(y).^2)/10^(25/10)/2)*(randn(nS,1)+1j*randn(nS,1)); y = y/sqrt(mean(abs(y).^2));
Asat = 1.3; pA = 2; z = y./(1+(abs(y)/Asat).^(2*pA)).^(1/(2*pA)).*exp(1j*0.15*(abs(y)/Asat).^2);   % PA Rapp + AM/PM
z = 10^(0.5/20)*real(z) + 1j*(imag(z)*cosd(3)+real(z)*sind(3));                                    % IQ imbalance 0.5 dB / 3 do
z = z + (0.005+0.005j);                                                                             % DC offset
z = z.*exp(1j*cumsum(sqrt(2*pi*50/fs)*randn(nS,1)));                                                % phase noise (Wiener, 50 Hz)
q = 2*max(abs([real(z); imag(z)]))/2^10; z = round(real(z)/q)*q + 1j*round(imag(z)/q)*q;            % ADC 10 bit
[f,p1] = welchPSD(y,fs); [~,p2] = welchPSD(z,fs);
oob = @(p) 10*log10(mean(p(abs(f)>52e6 & abs(f)<60e6))/mean(p(abs(f)<40e6)));
evmCh = evmGenie(carrier,y,grid,typ,pg,ti,symStart,symLen,Nfft,fs,scs,Nsc);
evmHw = evmGenie(carrier,z,grid,typ,pg,ti,symStart,symLen,Nfft,fs,scs,Nsc);
figure('Name','HW impairments','Position',[100 100 1400 360]);
subplot(1,3,1); plot(f/1e6,10*log10(p1),'b',f/1e6,10*log10(p2),'r'); legend(sprintf('Truoc (OOB %.0f dBc)',oob(p1)),sprintf('Sau HW (OOB %.0f dBc)',oob(p2))); xlabel('MHz'); title('Pho: spectral regrowth');
aa = linspace(0,3,200); subplot(1,3,2); plot(aa,aa,'k--',aa,aa./(1+(aa/Asat).^(2*pA)).^(1/(2*pA)),'r'); xlabel('|x_{in}|'); ylabel('|x_{out}|'); title('AM/AM (Rapp)');
subplot(1,3,3); text(0.1,0.6,sprintf('EVM sau kenh: %.1f%%\nEVM sau HW: %.1f%%',100*evmCh,100*evmHw),'FontSize',14); axis off; title('EVM (TDL-D, SNR 25 dB)');
saveas(gcf, fullfile(outDir,'tx_hw.png'));
fprintf('EVM AWGN10 = %.1f%% | EVM kenh = %.1f%% | EVM HW = %.1f%%\n',100*evm,100*evmCh,100*evmHw);
%ham phu
function [f,p] = welchPSD(x,fs)
nseg = 8192; w = 0.5*(1-cos(2*pi*(0:nseg-1).'/nseg)); nb = floor(numel(x)/nseg); p = zeros(nseg,1);
for b = 1:nb, s = x((b-1)*nseg+(1:nseg)).*w; p = p + abs(fft(s)).^2; end
p = fftshift(p/(nb*sum(w.^2)*fs)); f = (-nseg/2:nseg/2-1).'*fs/nseg;
end
function e = evmGenie(carrier,y,grid,typ,pg,ti,symStart,symLen,Nfft,fs,scs,Nsc)
% EVM sau can bang ZF voi kenh biet truoc (tinh tu path gain tai tam moi symbol) + sua CPE theo symbol
Y = nrOFDMDemodulate(carrier, y); tau = ti.PathDelays(:).'; fsc = ((0:Nsc-1).'-Nsc/2)*scs;
err = []; ref = [];
for s = 1:size(grid,2)
    tc = min(round(symStart(s)+symLen(s)-Nfft/2)+1, size(pg,1));
    H = exp(-1j*2*pi*fsc*tau)*squeeze(pg(tc,:,1,1)).';
    Zs = Y(:,s)./H; m = typ(:,s)==1; cpe = sum(conj(grid(m,s)).*Zs(m))/sum(abs(grid(m,s)).^2);
    Zs = Zs/cpe; err = [err; Zs(m)-grid(m,s)]; ref = [ref; grid(m,s)];   %#ok<AGROW>
end
e = sqrt(mean(abs(err).^2)/mean(abs(ref).^2));
end
