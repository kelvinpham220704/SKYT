%% 06 - Full-duplex: nhieu tu sinh (SI) va SIC.  (A) RDM truoc/sau SIC   (B) BER/BLER truoc/sau SIC
clear; clc; close all; isac_params;
%% (A) Sensing
nSlot = ceil(M/14); pdsch = isac_pdsch(1, NRB); X = zeros(Nsc, nSlot*14); rng(1);
for s = 1:nSlot, carrier.NSlot = mod(s-1,20); X(:,(s-1)*14+(1:14)) = isac_tx_slot(carrier, pdsch, 1); end
Xd = X(1:DEC:end,1:M);
TGT = [60 18 -35 -38; 120 -12 12 -41; 210 28 42 -45; 150 -30 -12 -43; 95 8 -40 -42; 175 6 22 -42];
base = struct('DF',DF,'Tsym',Tsym,'NR',NR,'NRK',NRK,'Nant',Nant,'clutter',false,'mti',true,'nlDbc',-50,'adc',true,'seed',5);
cases = {'before', 80, 0, false; 'after', 80, 45, true};     % SI +80 dB; sau: analog -45 dB + LS so
vAxis = ((1:M)-(M/2+1))*VBIN; rAxis = (0:NRK-1)*RBIN; figure('Position',[100 100 1200 420]);
for c_ = 1:2
    o = base; o.siDb = cases{c_,2}; o.analogDb = cases{c_,3}; o.sic = cases{c_,4};
    cube = isac_gen_cube(Xd, TGT, o); pw = sum(abs(cube).^2,3); det = isac_cfar(pw, Nant, 1e-7, 2, 8);
    subplot(1,2,c_); imagesc(vAxis,rAxis,10*log10(pw+1e-3),[60 95+10*(c_==1)]); axis xy; xlim([-60 60]); ylim([0 300]); colormap(jet); colorbar; hold on;
    [ri,di] = find(det); plot(vAxis(di),rAxis(ri),'wo'); title(sprintf('RDM %s SIC - san nhieu %.1f dB', cases{c_,1}, 10*log10(median(pw(:))))); xlabel('v [m/s]'); ylabel('R [m]');
end
saveas(gcf, fullfile(outDir,'sic_rdm.png'));
%% (B) Truyen thong (TDL-A, MMSE)
snrs = -4:2:12; nSlots = 60; names = {'none','nosic','analog','digital'}; BER = zeros(4,numel(snrs)); BLER = BER;
for k = 1:4, for s = 1:numel(snrs)
    [BER(k,s),BLER(k,s)] = isac_link('A', snrs(s), 'mmse', nSlots, struct('NRB',NRB,'siCase',names{k},'seed',1));
end, end
figure('Position',[100 100 1000 380]); lab = {'Khong SI (ly tuong)','Khong SIC (SI +80 dB)','Chi analog 45 dB','Analog + LS so'};
subplot(1,2,1); semilogy(snrs,max(BER',1e-6),'o-'); grid on; legend(lab); xlabel('SNR [dB]'); ylabel('BER'); title('Full-duplex: BER');
subplot(1,2,2); semilogy(snrs,max(BLER',1e-3),'o-'); grid on; xlabel('SNR [dB]'); ylabel('BLER'); title('Full-duplex: BLER');
saveas(gcf, fullfile(outDir,'sic_comm.png'));
