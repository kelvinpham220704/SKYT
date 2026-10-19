function det = isac_cfar(pw, Lint, pfa, g, t)
% CA-CFAR 2D tren cong suat da tich phan khong dong pha Lint anten
w = 2*(g+t)+1; kern = ones(w); kern(t+1:t+2*g+1, t+1:t+2*g+1) = 0;
Z = conv2(pw, kern, 'same')/nnz(kern);
alpha = gammaincinv(pfa, Lint, 'upper')/Lint;
det = pw > alpha*Z*1.15;
h = g+t; det([1:h end-h+1:end],:) = false; det(:,[1:h end-h+1:end]) = false;
end
