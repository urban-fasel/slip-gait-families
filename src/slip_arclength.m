function B = slip_arclength(o0, ds, nmax, zlim)
%SLIP_ARCLENGTH  Pseudo-arclength continuation of a gait in the swing damping zeta_s, through folds.
%   Unknowns [height; leg angle; leg rate; zeta_s/0.1].  Starts at o0 towards increasing zeta_s and stops
%   after nmax points or when zeta_s leaves zlim.  B(k): .w .zeta .turn (true once a fold has been passed)
h = 1e-30;  v = o0.v;  par = o0.par;  x = [o0.w(1:3); par(6)/0.1];  t = [0;0;0;1];
B = pack(x, false);  turned = false;
for k = 1:nmax
    [~, Jx] = resid(x, par, v, h);  tn = [Jx; t.'] \ [0;0;0;1];  tn = tn/norm(tn);
    if tn(4)*t(4) < 0, turned = true; end
    t = tn;  xp = x + ds*t;  xc = xp;  ok = false;
    for it = 1:12
        [r, Jx, good] = resid(xc, par, v, h);  if ~good, break; end
        q = [r; t.'*(xc - xp)];  if norm(q) < 1e-11, ok = true; break; end
        xc = xc - [Jx; t.'] \ q;
    end
    if ~ok, ds = ds/2;  if ds < 1e-5, break; end,  continue; end
    x = xc;  B(end+1) = pack(x, turned); %#ok<AGROW>
    if x(4)*0.1 < zlim(1) || x(4)*0.1 > zlim(2), break; end
end
end

function [r, Jx, good] = resid(x, par, v, h)
Jx = zeros(3,4);  p = par;  p(6) = x(4)*0.1;  z = [x(1); v; x(2); x(3)];
[zn, q] = slip_hop(z, 0, p);  r = zn(2:4) - z(2:4);  good = q.ok && real(q.clear) > 0;
if nargout < 2 || ~good, return; end
for i = 1:4
    xc = x;  xc(i) = xc(i) + 1i*h;  p = par;  p(6) = xc(4)*0.1;  z = [xc(1); v; xc(2); xc(3)];
    zn = slip_hop(z, 0, p);  Jx(:,i) = imag(zn(2:4) - z(2:4))/h;
end
end

function b = pack(x, turned)
b.w = x(1:3);  b.zeta = x(4)*0.1;  b.turn = turned;
end
