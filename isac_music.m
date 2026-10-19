function [ps, D] = isac_music(X, L, fb, Dfix, grid)
% MUSIC / FBSS-MUSIC.  X: N x S snapshot; L: kich thuoc subarray; fb=true -> lam tron thuan-nghich
N = size(X,1);  S = size(X,2);
if ~fb, L = N; end
P = N-L+1; J = fliplr(eye(L)); R = zeros(L);
for p = 1:P
    Xs = X(p:p+L-1,:); Rp = (Xs*Xs')/S;
    if fb, R = R + Rp + J*conj(Rp)*J; else, R = R + Rp; end
end
if fb, R = R/(2*P); else, R = R/P; end
[V,E] = eig((R+R')/2); [ev,ix] = sort(real(diag(E)),'descend'); V = V(:,ix);
if isempty(Dfix), D = max(1, isac_mdl(ev, max(S,8), 3)); else, D = Dfix; end
En = V(:,D+1:L); A = isac_steer(grid, L);
ps = 1./sum(abs(En'*A).^2, 1);
end
