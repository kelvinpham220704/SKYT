function D = isac_mdl(ev, Ns, Dmax)
ev = max(sort(ev,'descend'), eps); L = numel(ev); cost = zeros(Dmax+1,1);
for k = 0:Dmax
    g = exp(mean(log(ev(k+1:end)))); a = mean(ev(k+1:end));
    cost(k+1) = -Ns*(L-k)*log(g/a) + 0.5*k*(2*L-k)*log(Ns);
end
[~,i] = min(cost); D = i-1;
end
