%% ISAC 5G-A/6G - tham so he thong (monostatic, n78)
% ULA 32 anten, fc = 3.5 GHz, SCS = 30 kHz, BW = 100 MHz (273 RB)
c = 299792458;  fc = 3.5e9;  lambda = c/fc;
scsKHz = 30;  NRB = 273;  Nsc = NRB*12;  Nant = 32;
carrier = nrCarrierConfig;
carrier.SubcarrierSpacing = scsKHz;
carrier.NSizeGrid = NRB;
carrier.CyclicPrefix = 'normal';
ofdmInfo = nrOFDMInfo(carrier);
fs = ofdmInfo.SampleRate;           % 122.88 MHz
Nfft = ofdmInfo.Nfft;               % 4096
scs = scsKHz*1e3;
Tsym = mean(ofdmInfo.SymbolLengths)/fs;   % thoi gian 1 symbol (co CP) ~ 35.7 us
% --- tham so xu ly sensing (tiet kiem RAM)
M   = 512;        % so symbol / CPI (~18.3 ms)
DEC = 2;          % chi dung 1 sub-carrier trong 2 (BW van 98 MHz)
K   = Nsc/DEC;    % 1638
DF  = DEC*scs;    % 60 kHz
NR  = 2048;       % FFT range (zero-pad)
NRK = 330;        % giu 330 bin range (<= 400 m)
RBIN = c/(2*NR*DF);  VBIN = lambda/(2*M*Tsym);
fD  = (120/3.6)/lambda;   % Doppler cuc dai cua UE 120 km/h
outDir = 'out_matlab'; if ~exist(outDir,'dir'), mkdir(outDir); end
