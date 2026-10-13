function R = slip_control(F, par, S, vlist, nMC)
%SLIP_CONTROL  The controllers and the two robustness tests on both families, at the speeds vlist.
%   Laws, all acting through the swing input u:
%     none     u = 0
%     lqr      LQR on the apex state, saturated at umax
%     lqr3x    the same, saturated at 3 umax
%     lqrinf   the same, unsaturated
%     nmpc     one-step optimal law: the u that minimises the LQR cost-to-go after the next hop, computed
%              on the nonlinear map with the next ground height known (terrain preview); saturated at umax
%     nmpc3x   the same at 3 umax
%   Tests, 100-hop runs from the nominal gait:
%     pothole  deepest pothole cleared in one landing [cm], with the peak leg force capped at 1.5x the
%              gait's own peak
%     dropout  hops survived by all nMC runs on 0.5 cm rough ground once the measurement is lost from hop 30
%   R: .laws .vlist .rho .F0 (speeds x families) .pothole .dropout (laws x speeds x families).  A value
%   stored negative means the largest size tested was survived (printed as >=).
um = S.umax;
LAW = { 'none',   0,    'none'
        'lqr',    um,   'lqr'
        'lqr3x',  3*um, 'lqr'
        'lqrinf', inf,  'lqr'
        'nmpc',   um,   'nmpc'
        'nmpc3x', 3*um, 'nmpc' };
nl = size(LAW,1);  nv = numel(vlist);
R = struct('laws', {LAW(:,1).'}, 'vlist', vlist, 'nMC', nMC, 'rho', nan(nv,2), 'F0', nan(nv,2), ...
           'pothole', nan(nl,nv,2), 'dropout', nan(nl,nv,2));
for b = 1:2
    W = F.(sprintf('W%d', b));
    for iv = 1:nv
        v = vlist(iv);  [~, i] = min(abs(F.vs - v));
        o = slip_orbit(par, v, struct('w', W(:,i), 'v', v));  assert(o.ok, 'family %d lost at v = %.2f', b, v);
        [K, P, R.rho(iv,b)] = slip_lqr(o, par, S.Q, S.R);
        gait = struct('z', o.z, 'K', K, 'P', P);
        [~, q] = slip_hop(o.z, 0, par, 0, 150);  R.F0(iv,b) = real(q.Fpk);
        for j = 1:nl
            R.pothole(j,iv,b) = scan_pothole(gait, LAW(j,:), par, S, 1.5*R.F0(iv,b));
            R.dropout(j,iv,b) = scan_dropout(gait, LAW(j,:), par, S, nMC);
            fprintf('.');
        end
    end
end
end

function u = apply(law, Z, gait, par, hg)
% the input for every column of Z; hg is the ground height of this hop (preview for nmpc)
umx = law{2};  N = size(Z,2);
switch law{3}
    case 'none', u = zeros(1,N);
    case 'lqr',  u = -gait.K*(Z(2:4,:) - gait.z(2:4));  u = max(-umx, min(umx, u));
    case 'nmpc'                                     % grid over u, one parabolic refinement
        ug = linspace(-umx, umx, 21);  C = nan(numel(ug), N);
        for j = 1:numel(ug)
            Zn = slip_hop(Z, ug(j)*ones(1,N), par, hg, 150);
            d = real(Zn(2:4,:)) - gait.z(2:4);  d(:, any(~isfinite(d),1)) = 1e6;
            C(j,:) = sum(d.*(gait.P*d), 1);
        end
        C(~isfinite(C)) = 1e12;
        [~, j0] = min(C, [], 1);  u = ug(j0);
        in = j0 > 1 & j0 < numel(ug);
        if any(in)
            k = j0(in);  ii = find(in);
            c0 = C(sub2ind(size(C), k, ii));  cm = C(sub2ind(size(C), k-1, ii));  cp = C(sub2ind(size(C), k+1, ii));
            den = cm - 2*c0 + cp;  st = 0.5*(cm - cp)./den;  st(~isfinite(st) | abs(st) > 1 | den <= 0) = 0;
            u(ii) = ug(k) + st.*(ug(2) - ug(1));
        end
        u = max(-umx, min(umx, u));
end
end

function alive = roll(gait, law, par, S, Z0, HG, DROP, Fcap)
% run every column of Z0 over the hops of HG (nH x N ground heights).  DROP marks hops without a
% measurement (u = 0).  A run ends on a failed hop, foot slip, the leg bottoming out, stopping, or a peak
% leg force above Fcap.
if nargin < 8 || isempty(Fcap), Fcap = inf; end
[nH, N] = size(HG);  Z = Z0;  alive = true(1,N);
for k = 1:nH
    ix = find(alive);  if isempty(ix), break; end
    h = HG(k, ix);  d = DROP(k, ix);  u = zeros(1, numel(ix));
    if ~all(d), u(~d) = apply(law, Z(:,ix(~d)), gait, par, h(~d)); end
    [Zn, q] = slip_hop(Z(:,ix), u, par, h, 150);
    fp = real(q.Fpk);  fp(~isfinite(fp)) = inf;
    bad = ~q.ok | abs(real(q.al_td)) > S.afric | real(q.lmin) < S.lfloor | ...
          real(Zn(2,:)) <= 0.05 | ~isfinite(real(Zn(1,:))) | fp > Fcap;
    alive(ix(bad)) = false;
    Z(:,ix(~bad)) = real(Zn(:,~bad));
end
end

function hm = scan_pothole(gait, law, par, S, Fcap)
% every depth as its own column; the result is the last depth before the first failure
hs = 0.005:0.005:3.00;  N = numel(hs);  nH = 100;
HG = zeros(nH, N);  HG(5,:) = -hs;
alive = roll(gait, law, par, S, repmat(gait.z,1,N), HG, false(nH,N), Fcap);
j = find(~alive, 1);
if isempty(j), k = N; elseif j == 1, k = 0; else, k = j-1; end
if k == 0, hm = 0; else, hm = 100*hs(k); end
if k == N, hm = -hm; end
end

function n = scan_dropout(gait, law, par, S, nMC)
% largest n such that all nMC runs survive n hops without a measurement, lost from hop 30
nH = 100;  nmax = 40;  rng(11);
HG = repmat(0.005*randn(nH, nMC), 1, nmax);
DROP = false(nH, nMC*nmax);
for nn = 1:nmax, DROP(30:min(30+nn-1, nH), (nn-1)*nMC + (1:nMC)) = true; end
alive = roll(gait, law, par, S, repmat(gait.z,1,nMC*nmax), HG, DROP);
allok = all(reshape(alive, nMC, nmax), 1);
j = find(~allok, 1);
if isempty(j), n = -nmax; elseif j == 1, n = 0; else, n = j-1; end
end
