%racking: 
clear; clc; close all; rng(21);
dt = 0.1; N = 300; T = (0:N-1)'*dt; SR = 0.4; STH = 0.3*pi/180; BLIND = [12 14.5]; PD = 0.95; GATE = 13.8;
x = zeros(N,2); vel = zeros(N,2); x(1,:) = [80 -70]; psi = 60*pi/180; sp = 15; ph = zeros(N,1);
for k = 1:N
    t = T(k);
    if t < 8, w = 0; a = 0; ph(k) = 0; elseif t < 16, w = 12*pi/180; a = 0; ph(k) = 1; elseif t < 24, w = 0; a = 2; ph(k) = 2; else, w = 0; a = 0; ph(k) = 3; end
    vel(k,:) = sp*[cos(psi) sin(psi)];
    if k < N, x(k+1,:) = x(k,:) + vel(k,:)*dt; psi = psi + w*dt; sp = sp + a*dt; end
end
Z = cell(N,1);
for k = 1:N
    m = zeros(0,2);
    if ~(T(k) >= BLIND(1) && T(k) < BLIND(2)) && rand < PD
        m = [hypot(x(k,1),x(k,2)) + SR*randn, atan2(x(k,2),x(k,1)) + STH*randn];
    end
    for i = 1:poisson1(1.5), pt = [-50+470*rand, -100+520*rand]; m(end+1,:) = [hypot(pt(1),pt(2)), atan2(pt(2),pt(1))]; end   %#ok<SAGROW>
    Z{k} = m;
end
[xk, Pk, sk] = runTracker('kf',  Z, x, vel, dt, SR, STH, GATE);
[xi, Pi, si] = runTracker('imm', Z, x, vel, dt, SR, STH, GATE);
ek = hypot(xk(:,1)-x(:,1), xk(:,2)-x(:,2)); ei = hypot(xi(:,1)-x(:,1), xi(:,2)-x(:,2));
rm = @(e,m) sqrt(mean(e(m).^2)); t1 = T > 1;
fprintf('RMSE [m]        KF-CV    IMM-KF\n');
names = {'Thang 0-8 s','Re 8-16 s','Tang toc 16-24 s','Thang 24-30 s'};
for p = 0:3, mk = (ph==p) & t1; fprintf('%-16s %7.2f %7.2f\n', names{p+1}, rm(ek,mk), rm(ei,mk)); end
bl = T >= BLIND(1) & T < BLIND(2); fprintf('%-16s %7.2f %7.2f\n','Vung mu', rm(ek,bl), rm(ei,bl));
fprintf('%-16s %7.2f %7.2f\n','Tong', rm(ek,t1), rm(ei,t1));
allm = []; for k = 1:N, if ~isempty(Z{k}), allm = [allm; Z{k}(:,1).*cos(Z{k}(:,2)) Z{k}(:,1).*sin(Z{k}(:,2))]; end, end   %#ok<AGROW>
outDir = 'out_matlab'; if ~exist(outDir,'dir'), mkdir(outDir); end
figure('Position',[100 100 1250 430]);
subplot(1,2,1); plot(allm(:,1),allm(:,2),'.','Color',[.8 .8 .8]); hold on; plot(x(:,1),x(:,2),'k','LineWidth',3); plot(xk(:,1),xk(:,2),'Color',[.9 .6 0],'LineWidth',1.4); plot(xi(:,1),xi(:,2),'r','LineWidth',1.6);
axis equal; xlim([40 460]); ylim([-110 330]); legend('Do + clutter','Quy dao goc','KF (CV)','IMM-KF'); xlabel('x [m]'); ylabel('y [m]'); title('Muc tieu phi tuyen: re roi tang toc');
subplot(1,2,2); semilogy(T,ek,'Color',[.9 .6 0]); hold on; semilogy(T,ei,'r','LineWidth',1.4); xline(BLIND(1),'--'); xline(BLIND(2),'--'); ylim([0.1 80]); legend('KF (CV)','IMM-KF'); xlabel('s'); ylabel('Sai so vi tri [m]'); title('Sai so theo thoi gian (vung mu 12-14.5 s)');
saveas(gcf, fullfile(outDir,'trk_maneuver.png'));

fig = figure('Position',[100 100 700 450]); fn = fullfile(outDir,'tracking_matlab.gif'); first = true;
for k = 11:2:N
    clf; hold on; plot(x(1:k,1),x(1:k,2),'k','LineWidth',2.5); plot(xk(1:k,1),xk(1:k,2),'Color',[.9 .6 0]); plot(xi(1:k,1),xi(1:k,2),'r','LineWidth',1.6);
    for j = max(1,k-25):k, if ~isempty(Z{j}), plot(Z{j}(:,1).*cos(Z{j}(:,2)), Z{j}(:,1).*sin(Z{j}(:,2)),'.','Color',[.6 .6 .6]); end, end
    plot(x(k,1),x(k,2),'ko','MarkerFaceColor','k'); axis equal; xlim([40 460]); ylim([-110 330]); title(sprintf('t = %.1f s',T(k)));
    if T(k) >= BLIND(1) && T(k) < BLIND(2), text(150,120,'VUNG MU - chi du doan','Color','r','FontSize',14,'FontWeight','bold'); end
    fr = getframe(fig); [ind,map] = rgb2ind(fr.cdata,256);
    if first, imwrite(ind,map,fn,'gif','LoopCount',Inf,'DelayTime',0.08); first = false; else, imwrite(ind,map,fn,'gif','WriteMode','append','DelayTime',0.08); end
end

function n = poisson1(lam), L = exp(-lam); n = 0; p = rand; while p > L, n = n+1; p = p*rand; end, end
function [z,Rc] = conv2c(m, SR, STH)
R = m(1); th = m(2); z = [R*cos(th); R*sin(th)]; J = [cos(th) -R*sin(th); sin(th) R*cos(th)]; Rc = J*diag([SR^2 STH^2])*J';
end
function Q = Qm(q, dt)
b = q*[dt^3/3 dt^2/2; dt^2/2 dt]; Q = zeros(4); for i = 1:2, Q([i i+2],[i i+2]) = b; end
end
function [xs, Ps, stat] = runTracker(kind, Z, x, vel, dt, SR, STH, GATE)
N = size(x,1); xs = zeros(N,4); Ps = zeros(4,4,N); stat = zeros(N,1); started = false;
H = [1 0 0 0; 0 1 0 0]; F = eye(4); F(1,3) = dt; F(2,4) = dt;
for k = 1:N
    zs = Z{k};
    if ~started
        if ~isempty(zs)   
            pts = [zs(:,1).*cos(zs(:,2)) zs(:,1).*sin(zs(:,2))]; [~,i] = min(hypot(pts(:,1)-x(1,1), pts(:,2)-x(1,2)));
            x0 = [pts(i,:) vel(1,:)].'; started = true; xs(k,:) = x0.'; Ps(:,:,k) = eye(4)*4;
            if strcmp(kind,'kf'), xe = x0; Pe = eye(4)*4; else, mu = [.5 .3 .2]; XS = {x0,x0,x0}; PS = {eye(4)*4,eye(4)*4,eye(4)*4}; end
        end
        continue
    end
    if strcmp(kind,'kf')
        xp = F*xe; Pp = F*Pe*F' + Qm(0.5,dt); best = []; bm = inf;
        for j = 1:size(zs,1)
            [z,Rc] = conv2c(zs(j,:),SR,STH); S = H*Pp*H' + Rc; r = z - H*xp; m = r'*(S\r);
            if m < GATE && m < bm, bm = m; best = {z,Rc,S,r}; end
        end
        if ~isempty(best), Kg = Pp*H'/best{3}; xe = xp + Kg*best{4}; Pe = (eye(4)-Kg*H)*Pp; stat(k) = 0; else, xe = xp; Pe = Pp; stat(k) = 1; end
        xs(k,:) = xe.'; Ps(:,:,k) = Pe;
    else
        qs = [0.3 6 40]; Pi = [.94 .04 .02; .05 .9 .05; .03 .12 .85];
        cbar = (Pi.'*mu.').'; mix = (Pi.*mu.')./cbar;           % mix(i,j)
        X0 = cell(1,3); P0 = cell(1,3);
        for j = 1:3
            X0{j} = zeros(4,1); for i = 1:3, X0{j} = X0{j} + mix(i,j)*XS{i}; end
            P0{j} = zeros(4); for i = 1:3, dX = XS{i}-X0{j}; P0{j} = P0{j} + mix(i,j)*(PS{i} + dX*dX'); end
        end
        XP = cell(1,3); PP = cell(1,3); xc = zeros(4,1);
        for j = 1:3, XP{j} = F*X0{j}; PP{j} = F*P0{j}*F' + Qm(qs(j),dt); xc = xc + cbar(j)*XP{j}; end
        Pc = zeros(4); for j = 1:3, dX = XP{j}-xc; Pc = Pc + cbar(j)*(PP{j} + dX*dX'); end
        best = []; bm = inf;
        for j = 1:size(zs,1)
            [z,Rc] = conv2c(zs(j,:),SR,STH); S = H*Pc*H' + Rc; r = z - H*xc; m = r'*(S\r);
            if m < GATE && m < bm, bm = m; best = {z,Rc}; end
        end
        if ~isempty(best)
            lik = zeros(1,3); XN = cell(1,3); PN = cell(1,3);
            for j = 1:3
                S = H*PP{j}*H' + best{2}; Kg = PP{j}*H'/S; r = best{1} - H*XP{j};
                XN{j} = XP{j} + Kg*r; PN{j} = (eye(4)-Kg*H)*PP{j}; lik(j) = exp(-.5*r'*(S\r))/sqrt(det(2*pi*S));
            end
            mu = cbar.*lik; mu = mu/sum(mu); mu = max(mu,1e-3); mu = mu/sum(mu); XS = XN; PS = PN; stat(k) = 0;
        else
            XS = XP; PS = PP; mu = cbar; stat(k) = 1;      % mat dau: chi du doan (coast)
        end
        xm = zeros(4,1); for j = 1:3, xm = xm + mu(j)*XS{j}; end
        Pm = zeros(4); for j = 1:3, dX = XS{j}-xm; Pm = Pm + mu(j)*(PS{j} + dX*dX'); end
        xs(k,:) = xm.'; Ps(:,:,k) = Pm;
    end
end
end
