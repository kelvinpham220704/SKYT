%% 05 - Truyen thong: ZF vs MMSE tren 4 kenh TDL-A/B/C/D, BER (sau giai ma) va BLER
% Can 5G Toolbox (nrDLSCH, nrPDSCH, nrTDLChannel...). Chinh nSlots (do chinh xac) va NRB (toc do) neu chay lau.
clear; clc; close all; isac_params;
snrs = -4:2:12; nSlots = 60; o = struct('NRB',NRB,'siCase','none','seed',1);   % NRB = 273 (100 MHz)
eqs = {'zf','mmse'}; BER = zeros(4,numel(eqs),numel(snrs)); BLER = BER;
for pi_ = 1:4
    abc = 'ABCD'; prof = abc(pi_);
    for e = 1:2
        for s = 1:numel(snrs)
            [BER(pi_,e,s), BLER(pi_,e,s)] = isac_link(prof, snrs(s), eqs{e}, nSlots, o);
        end
    end
    figure('Name',['TDL-' prof],'Position',[100 100 1000 380]);
    subplot(1,2,1); semilogy(snrs,max(squeeze(BER(pi_,1,:)),1e-6),'s-','Color',[.9 .6 0]); hold on; semilogy(snrs,max(squeeze(BER(pi_,2,:)),1e-6),'ro-'); grid on; legend('ZF','MMSE'); xlabel('SNR [dB]'); ylabel('BER'); title(['TDL-' prof ': BER']);
    subplot(1,2,2); semilogy(snrs,max(squeeze(BLER(pi_,1,:)),1e-3),'s-','Color',[.9 .6 0]); hold on; semilogy(snrs,max(squeeze(BLER(pi_,2,:)),1e-3),'ro-'); grid on; legend('ZF','MMSE'); xlabel('SNR [dB]'); ylabel('BLER'); title(['TDL-' prof ': BLER']);
    saveas(gcf, fullfile(outDir,['comm_tdl_' prof '.png']));
end
save(fullfile(outDir,'comm_zf_mmse.mat'),'snrs','BER','BLER');
