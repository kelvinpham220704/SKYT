%Do phan giai goc
clear; clc; close all; isac_params; rng(11);
grid = -80:0.05:80; Kn = 9; NT = 100; wk = [1 .6 .6 .6 .6 .3 .3 .3 .3];
dths = [0.5 0.75 1 1.25 1.5 2 2.5 3 3.5 4 5 6]; snr = 15; methods = {'bart','music','fbss'};
P = zeros(numel(methods), numel(dths));
for mi = 1:numel(methods)
    for di = 1:numel(dths)
        ok = 0; for t = 1:NT, ok = ok + trial(dths(di), snr, methods{mi}, Nant, grid, Kn, wk); end
        P(mi,di) = ok/NT;
    end
end
figure; plot(dths,P(1,:),'ko-',dths,P(2,:),'o-','Color',[0.9 0.6 0]); hold on; plot(dths,P(3,:),'ro-','LineWidth',1.5); yline(0.9,':');
xlabel('Khoang cach goc \Delta\theta [do]'); ylabel('Xac suat phan giai'); legend('Bartlett','MUSIC','FBSS-MUSIC'); title('SNR 15 dB / anten');
saveas(gcf, fullfile(outDir,'rx_resolution.png'));
thr = dths(find(P(3,:) >= 0.9, 1)); fprintf('FBSS-MUSIC phan giai duoc tu \\Delta\\theta = %.2f do (Rayleigh = %.2f do)\n', thr, rad2deg(2/Nant));

function ok = trial(dth, snr, method, Nant, grid, Kn, wk)
t1 = 23-dth/2; t2 = 23+dth/2; A = isac_steer([t1 t2], Nant);
s = exp(1j*2*pi*rand(2,1)); sig = (A*s)*wk;                       % 2 nguon ket hop
sg = 10^(-snr/20)/sqrt(2); X = sig + sg*(randn(Nant,Kn)+1j*randn(Nant,Kn));
switch method
    case 'bart', R = X*X'/Kn; Ab = isac_steer(grid,Nant); ps = real(sum(conj(Ab).*(R*Ab),1));
    case 'music', ps = isac_music(X, Nant, false, 2, grid);
    otherwise, ps = isac_music(X, 24, true, 2, grid);
end
i = find(ps(2:end-1) > ps(1:end-2) & ps(2:end-1) >= ps(3:end)) + 1; [~,o] = sort(ps(i),'descend'); i = i(o(1:min(2,numel(o)))); est = sort(grid(i));
ok = numel(est) == 2 && abs(est(1)-t1) < dth/2 && abs(est(2)-t2) < dth/2;
end
