function o = slip_orbit(par, v, warm)
%SLIP_ORBIT  A periodic gait at forward speed v: a fixed point of the apex map, found by Newton on
%   w = [height; leg angle; leg rate].  The model conserves energy, so the speed row of the map follows.
%   warm: a nearby solution to start from.  Without one, the gait is started where forward hopping
%   branches off hopping in place, and the speed is stepped up to v.
%   o: .w .v .z = [y; v; al; ald] .par .ok .info
if nargin < 3 || isempty(warm)
    yc = fzero(@(y) detAI(y, par), bracket(@(y) detAI(y, par), par(3)*linspace(1.05, 2.0, 60)));
    [~, A] = detAI(yc, par);  [V, D] = eig(A);  [~, i] = min(abs(diag(D) - 1));  nv = real(V(:,i))/real(V(1,i));
    warm = struct('w', [yc; 0; 0], 'v', 0, 'nv', nv);
    for a = [0.05 0.15 0.3 0.5 0.75 1:0.25:max(v,1)]
        if a > v, break; end
        warm = slip_orbit(par, a, warm);  assert(warm.ok, 'branch start failed at v = %.2f', a);
    end
end
w = warm.w(:);  if isfield(warm, 'nv'), w(2:3) = w(2:3) + (v - warm.v)*warm.nv(2:3); end
h = 1e-30;
for it = 1:30
    [r, info] = res(w);  J = zeros(3);
    for i = 1:3, e = zeros(3,1);  e(i) = 1i*h;  J(:,i) = imag(res(w + e))/h; end
    if ~info.ok || norm(r) < 1e-12, break; end
    w = w - J\r;
end
o = struct('w', w, 'v', v, 'z', [w(1); v; w(2); w(3)], 'par', par, 'ok', info.ok && norm(r) < 1e-9, 'info', info);
    function [r, info] = res(w)
        z = [w(1); v; w(2); w(3)];  [zn, info] = slip_hop(z, 0, par);  r = zn(2:4) - z(2:4);
    end
end

function [d, A] = detAI(y, par)
% hopping in place at apex height y; det(A - I) = 0 where forward hopping branches off
h = 1e-30;  A = zeros(3);  z = [y; 0; 0; 0];  H = par(1)*par(4)*y;
for i = 1:3
    zc = z;  zc(i+1) = zc(i+1) + 1i*h;  zc(1) = (H - 0.5*par(1)*zc(2)^2)/(par(1)*par(4));
    zn = slip_hop(zc, 0, par);  A(:,i) = imag(zn(2:4))/h;
end
d = det(A - eye(3));
end

function b = bracket(f, ys)
fv = arrayfun(f, ys);  i = find(fv(1:end-1).*fv(2:end) < 0, 1);
assert(~isempty(i), 'no bifurcation found');  b = ys(i:i+1);
end
