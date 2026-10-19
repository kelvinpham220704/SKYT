function [txGrid, pdschInd, pdschInfo, dmrsInd] = isac_tx_slot(carrier, pdsch, nTx, cw)
% Tao resource grid ISAC-OFDM cho 1 slot: PDSCH (QPSK) + DM-RS. Moi RE deu duoc dung lam tham chieu sensing.
[pdschInd, pdschInfo] = nrPDSCHIndices(carrier, pdsch);
if nargin < 4 || isempty(cw), cw = randi([0 1], pdschInfo.G, 1); end
sym = nrPDSCH(carrier, pdsch, cw);
txGrid = nrResourceGrid(carrier, nTx);
txGrid(pdschInd) = sym;
dmrsInd = nrPDSCHDMRSIndices(carrier, pdsch);
txGrid(dmrsInd) = nrPDSCHDMRS(carrier, pdsch);
end
