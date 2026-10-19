function [ber, bler, sinrDb] = isac_link(profile, snrDb, eqType, nSlots, o)
% Lien ket 4x4 MIMO (4 lop), QPSK, LDPC 5G NR (R ~ 0.48), kenh nrTDLChannel, can bang ZF/MMSE.
% o.NRB (so RB), o.siCase: 'none' | 'nosic' | 'analog' | 'digital' (nhieu tu sinh full-duplex), o.seed
if ~isfield(o,'siCase'), o.siCase = 'none'; end
if ~isfield(o,'seed'), o.seed = 1; end
rng(o.seed); nTx = 4; nRx = 4; NRB = o.NRB; nVar = 10^(-snrDb/10);
carrier = nrCarrierConfig; carrier.SubcarrierSpacing = 30; carrier.NSizeGrid = NRB;
fs = nrOFDMInfo(carrier).SampleRate; pdsch = isac_pdsch(nTx, NRB);
DS = 300e-9; if strcmp(profile,'D'), DS = 100e-9; end
tdl = nrTDLChannel; tdl.DelayProfile = ['TDL-' profile]; tdl.DelaySpread = DS; tdl.MaximumDopplerShift = 389;
tdl.SampleRate = fs; tdl.NumTransmitAntennas = nTx; tdl.NumReceiveAntennas = nRx; tdl.Seed = o.seed;
pad = info(tdl).MaximumChannelDelay; tcr = 490/1024;
enc = nrDLSCH; enc.MultipleHARQProcesses = false; enc.TargetCodeRate = tcr;
dec = nrDLSCHDecoder; dec.MultipleHARQProcesses = false; dec.TargetCodeRate = tcr;
dec.LDPCDecodingAlgorithm = 'Normalized min-sum'; dec.MaximumLDPCIterationCount = 8;
nBitErr = 0; nBits = 0; nBlkErr = 0; sinrDb = [];
for slot = 0:nSlots-1
    carrier.NSlot = mod(slot,20);
    [pdschInd, pInfo] = nrPDSCHIndices(carrier, pdsch);
    trBlk = nrTBS(pdsch.Modulation, pdsch.NumLayers, numel(pdsch.PRBSet), pInfo.NREPerPRB, tcr, 0);
    dec.TransportBlockLength = trBlk;
    tb = randi([0 1], trBlk, 1); setTransportBlock(enc, tb);
    cw = enc(pdsch.Modulation, pdsch.NumLayers, pInfo.G, 0);
    txGrid = isac_tx_slot(carrier, pdsch, nTx, cw);
    txWave = nrOFDMModulate(carrier, txGrid);
    [rxWave, pg, st] = tdl([txWave; zeros(pad,nTx)]);
    pf = getPathFilters(tdl); off = nrPerfectTimingEstimate(pg, pf);
    rxWave = rxWave(1+off:end,:);
    Hest = nrPerfectChannelEstimate(carrier, pg, pf, off, st);
    rxGrid = nrOFDMDemodulate(carrier, rxWave);
    rxGrid = rxGrid + sqrt(nVar/2)*(randn(size(rxGrid)) + 1j*randn(size(rxGrid)));   % AWGN / RE
    [rxSym, Hre] = nrExtractResources(pdschInd, rxGrid, Hest);
    if ~strcmp(o.siCase,'none'), rxSym = rxSym + siTerm(o.siCase, size(rxSym), nVar); end
    [eq, csi] = isac_eq_mimo(rxSym, Hre, nVar, eqType);
    [llr, rxSymbols] = nrPDSCHDecode(carrier, pdsch, eq, 1);
    csiL = nrLayerDemap(csi);
    Qm = length(llr{1})/length(rxSymbols{1}); c1 = repmat(csiL{1}.', Qm, 1); llr{1} = llr{1}.*c1(:);
    resetSoftBuffer(dec, 0);
    [bits, blkErr] = dec(llr, pdsch.Modulation, pdsch.NumLayers, 0);
    nBitErr = nBitErr + sum(double(bits) ~= tb); nBits = nBits + trBlk; nBlkErr = nBlkErr + double(blkErr);
    if slot == 0, sinrDb = 10*log10(csi(:)); end
end
ber = nBitErr/nBits; bler = nBlkErr/nSlots;
end

function e = siTerm(sicCase, sz, nVar)
% Nhieu tu sinh (SI) tren moi RE: +80 dB (khong SIC), analog -45 dB, analog + LS so (16 pilot)
si0 = 80; ana = 0; if ~strcmp(sicCase,'nosic'), ana = 45; end
amp = 10^((si0-ana)/20); cn = @() (randn(sz)+1j*randn(sz))/sqrt(2);
nl = amp*10^(-50/20)*cn(); adc = amp/10^(61.9/20)*cn();
xsi = (sign(randn(sz(1),1)) + 1j*sign(randn(sz(1),1)))/sqrt(2);
if strcmp(sicCase,'digital')
    errv = (nVar + (amp*10^(-50/20))^2 + (amp/10^(61.9/20))^2)/16;   % sai so LS tren 16 pilot + thanh phan khong khu duoc
    e = -sqrt(errv)*cn().*xsi + nl + adc;
else
    e = amp*cn().*xsi + nl + adc;
end
end
