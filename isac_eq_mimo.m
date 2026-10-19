function [x, csi] = isac_eq_mimo(Y, H, nVar, type)
% Can bang MIMO theo tung RE: 'zf' hoac 'mmse'.  Y: N x nRx, H: N x nRx x nTx
[N, nRx, nTx] = size(H); x = zeros(N,nTx); csi = zeros(N,nTx);
for n = 1:N
    Hn = reshape(H(n,:,:), nRx, nTx); y = Y(n,:).';
    G = Hn'*Hn;
    if strcmp(type,'zf')
        Gi = inv(G); W = Gi*Hn'; x(n,:) = (W*y).'; v = nVar*real(diag(Gi));
    else
        W = (G + nVar*eye(nTx))\Hn'; B = real(diag(W*Hn)); x(n,:) = ((W*y)./B).'; v = (1-B)./B;
    end
    csi(n,:) = (1./max(v,1e-9)).';
end
end
