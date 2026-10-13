function P = slip_period2(F, par)
%SLIP_PERIOD2  The period-2 gait born where family 2 stops being self-stable, followed in speed from the
%   first grid point past that speed, at the energy of the family-2 gait there.  Unknowns: one apex
%   [forward speed; leg angle; leg rate]; equations: the two-hop map returns to it.
%   P: .v (mean running speed of the gait) .vnom (speed of the family-2 gait of equal energy) .rho (error
%   growth per hop, the square root of the two-hop multiplier) .za .zb (the two apex states, faster one
%   first).  Empty if family 2 stays self-stable on the grid.
m = par(1);  g = par(4);  h = 1e-30;  vs = F.vs;
P = struct('v', [], 'vnom', [], 'rho', [], 'za', zeros(4,0), 'zb', zeros(4,0));
k = find(F.rho2 >= 1, 1);  if isempty(k) || k == 1, return; end
w = [];
for i = k:numel(vs)
    o = slip_orbit(par, vs(i), struct('w', F.W2(:,i), 'v', vs(i)));  z = o.z;  H = 0.5*m*z(2)^2 + m*g*z(1);
    if isempty(w)                                   % start: let the unstable family-2 gait fall onto the period-2 gait
        A = slip_lin(o, par);  [V, D] = eig(A);  [~, j] = max(abs(diag(D)));  e = real(V(:,j));  e = e/norm(e);
        Z = z;  Z(2:4) = z(2:4) + 1e-3*e;  Z(1) = (H - 0.5*m*Z(2)^2)/(m*g);
        for q = 1:600, [Z, info] = slip_hop(Z, 0, par);  if ~info.ok, return; end, end
        w = Z(2:4);
    end
    res = @(w) two_hop(w, H, m, g, par);  ok = false;
    for it = 1:40
        r = res(w);  if any(~isfinite(r)), break; end
        J = zeros(3);  for c = 1:3, ec = zeros(3,1);  ec(c) = 1i*h;  J(:,c) = imag(res(w + ec))/h; end
        if norm(r) < 1e-12, ok = true; break; end,  w = w - J\r;
    end
    if ~ok, return; end
    za = [(H - 0.5*m*w(1)^2)/(m*g); w];  zb = slip_hop(za, 0, par);
    if za(2) < zb(2), za = zb; end                                      % faster apex first
    [zb, qa] = slip_hop(za, 0, par);  [~, qb] = slip_hop(zb, 0, par);
    mu = max(abs(eig(slip_lin(struct('z', zb), par)*slip_lin(struct('z', za), par))));
    P.vnom(end+1) = vs(i);  P.v(end+1) = (qa.dist + qb.dist)/(qa.T + qb.T);  P.rho(end+1) = sqrt(mu);
    P.za(:,end+1) = za;  P.zb(:,end+1) = zb;
end
end

function r = two_hop(w, H, m, g, par)
z1 = [(H - 0.5*m*w(1)^2)/(m*g); w(:)];
z2 = slip_hop(z1, 0, par);  z3 = slip_hop(z2, 0, par);
r = z3(2:4) - z1(2:4);
end
