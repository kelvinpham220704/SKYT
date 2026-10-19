function pdsch = isac_pdsch(nLayers, NRB)
% Cau hinh PDSCH + DM-RS (type 1, comb-2) cho 1 slot, toan bo RB
pdsch = nrPDSCHConfig;
pdsch.Modulation = 'QPSK';
pdsch.NumLayers = nLayers;
pdsch.PRBSet = 0:NRB-1;
pdsch.SymbolAllocation = [0 14];
pdsch.MappingType = 'A';
pdsch.DMRS.DMRSConfigurationType = 1;
pdsch.DMRS.DMRSTypeAPosition = 2;
pdsch.DMRS.DMRSAdditionalPosition = 2;
pdsch.DMRS.DMRSLength = 1;
pdsch.DMRS.NumCDMGroupsWithoutData = 1 + (nLayers > 2);
pdsch.DMRS.DMRSPortSet = 0:nLayers-1;
end
