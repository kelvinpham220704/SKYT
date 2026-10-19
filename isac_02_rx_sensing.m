%Rx monostatic:
clear; clc; close all; isac_params;
nSlot = ceil(M/14); pdsch = isac_pdsch(1, NRB); X = zeros(Nsc, nSlot*14);
rng(1);
for s = 1:nSlot
    carrier.NSlot = mod(s-1,20); g = isac_tx_slot(carrier, pdsch, 1);
    X(:,(s-1)*14+(1:14)) = g;
end
Xd = X(1:DEC:end, 1:M);                                   % K x M (da decimate)
TGT = [60 18 -35 -38; 120 -12 12 -41; 210 28 42 -45; 150 -30 -12 -43; 95 8 -40 -42; 175 6.0 22.0 -42; 175 6.3 24.5 -42];
o = struct('DF',DF,'Tsym',Tsym,'NR',NR,'NRK',NRK,'Nant',Nant,'clutter',true,'mti',true,'siDb',[],'analogDb',0,'sic',false,'nlDbc',-50,'adc',true,'seed',3);
tic; [cube, rawpw] = isac_gen_cube(Xd, TGT, o); fprintf('Tao cube: %.1f s, %.0f MB\n', toc, numel(cube)*8/1e6);
pw = sum(abs(cube).^2, 3);                                
det = isac_cfar(pw, Nant, 1e-7, 2, 8);
mx = pw; for s = -3:3, mx = max(mx, circshift(pw,[s 0])); end
mx2 = mx; for s = -3:3, mx2 = max(mx2, circshift(mx,[0 s])); end
[ri, di] = find(det & (pw==mx2) & abs((1:M)-(M/2+1)) > 1);
vAxis = ((1:M)-(M/2+1))*VBIN; rAxis = (0:NRK-1)*RBIN;
%RDM tho va RDM sau khu clutter + CFAR
figure('Position',[100 100 1250 450]);
subplot(1,2,1); imagesc(vAxis,rAxis,10*log10(double(rawpw)+1e-3),[62 102]); axis xy; xlim([-60 60]); ylim([0 300]); colormap(jet); colorbar; xlabel('Van toc [m/s]'); ylabel('Cu ly [m]'); title('RDM tho (co clutter)');
subplot(1,2,2); imagesc(vAxis,rAxis,10*log10(pw+1e-3),[62 90]); axis xy; xlim([-60 60]); ylim([0 300]); colorbar; xlabel('Van toc [m/s]'); ylabel('Cu ly [m]'); title('RDM sau khu clutter + CA-CFAR'); hold on;
plot(vAxis(di), rAxis(ri), 'wo', 'MarkerSize', 13, 'LineWidth', 1.5); plot(TGT(:,2), TGT(:,1), 'kx'); 
saveas(gcf, fullfile(outDir,'rx_rdm.png'));
figure; [VV,RR] = meshgrid(vAxis(abs(vAxis)<50), rAxis(rAxis<260)); sub = 10*log10(pw(rAxis<260, abs(vAxis)<50)+1e-3);
surf(VV, RR, sub, 'EdgeColor','none'); colormap(jet); caxis([62 90]); zlim([60 92]); xlabel('v [m/s]'); ylabel('R [m]'); zlabel('dB'); title('Range-Doppler Map 3D'); view(32,32);
saveas(gcf, fullfile(outDir,'rx_rdm3d.png'));
%Uoc luong R, v va AoA 
grid = -80:0.05:80; res = {}; figure('Position',[100 100 1250 400]); subplot(1,2,1); hold on;
for i = 1:numel(ri)
    r = ri(i); d = di(i); lp = log(pw+1e-30);
    dr = par3(lp(r-1,d),lp(r,d),lp(r+1,d)); dd = par3(lp(r,d-1),lp(r,d),lp(r,d+1));
    Rhat = (r-1+dr)*RBIN; vhat = (d-(M/2+1)+dd)*VBIN;
    Xs = []; for a = -1:1, for b = -1:1, Xs = [Xs squeeze(cube(r+a,d+b,:))]; end, end   %#ok<AGROW>
    [ps, D] = isac_music(Xs, 24, true, [], grid); ang = pickPeaks(ps, grid, D);
    plot(grid, 10*log10(ps/max(ps)));
    fprintf('R = %7.2f m  v = %6.2f m/s  theta = %s deg (D=%d)\n', Rhat, vhat, mat2str(ang,4), D);
end
xlim([-60 60]); ylim([-60 3]); xlabel('Goc [do]'); ylabel('MUSIC [dB]'); title('FBSS-MUSIC tung dinh');
% cap muc tieu cung o: Bartlett vs MUSIC vs FBSS-MUSIC
[~,k] = min(abs(rAxis(ri)-175)); Xs = []; for a = -1:1, for b = -1:1, Xs = [Xs squeeze(cube(ri(k)+a,di(k)+b,:))]; end, end
Rb = (Xs*Xs')/size(Xs,2); Ab = isac_steer(grid,Nant); bt = real(sum(conj(Ab).*(Rb*Ab),1));
pm = isac_music(Xs, Nant, false, 2, grid); pf = isac_music(Xs, 24, true, [], grid);
subplot(1,2,2); hold on;
plot(grid(:), 10*log10(bt(:)/max(bt)), 'k', 'LineWidth', 1.5);
plot(grid(:), 10*log10(pm(:)/max(pm)), 'Color', [0.9 0.6 0], 'LineWidth', 1.5);
plot(grid(:), 10*log10(pf(:)/max(pf)), 'r', 'LineWidth', 1.5);
xline(22,'--'); xline(24.5,'--'); xlim([12 34]); ylim([-45 3]); legend('Bartlett','MUSIC','FBSS-MUSIC'); xlabel('Goc [do]'); title('2 muc tieu cung o RD (22 va 24.5 do)');
saveas(gcf, fullfile(outDir,'rx_music_pair.png'));

function d = par3(y0,y1,y2), den = y0-2*y1+y2; d = 0; if den ~= 0, d = 0.5*(y0-y2)/den; end, end
function ang = pickPeaks(ps, grid, D)
i = find(ps(2:end-1) > ps(1:end-2) & ps(2:end-1) >= ps(3:end)) + 1; [~,o] = sort(ps(i),'descend'); i = i(o(1:min(D,numel(o)))); ang = sort(grid(i));
end
