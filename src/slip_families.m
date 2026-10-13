function F = slip_families(par, vs)
%SLIP_FAMILIES  The two families of natural gaits on the same hardware, over the speeds vs.
%   Family 1 is followed in speed from where forward hopping branches off hopping in place.  Family 2 is
%   reached from family 1 by continuation in the swing damping: through a fold and back to zero damping,
%   where it is a different gait; it is then followed in speed.
%   F: .vs .par .W1 .W2 (gait unknowns [height; leg angle; leg rate]) .rho1 .rho2 .Fpk1 .Fpk2
n = numel(vs);

% family 1
O1 = cell(1,n);  warm = [];
for i = 1:n
    o = slip_orbit(par, vs(i), warm);  assert(o.ok, 'family 1 lost at v = %.2f', vs(i));
    O1{i} = o;  warm = struct('w', o.w, 'v', o.v);
end

% family 2: through the fold once, then follow in speed with a secant predictor
O2 = cell(1,n);  O2{1} = second_family(par, vs(1));  wprev = [];
for i = 2:n
    wg = O2{i-1}.w;  if ~isempty(wprev), wg = 2*O2{i-1}.w - wprev; end
    o = slip_orbit(par, vs(i), struct('w', wg, 'v', vs(i)));
    if ~o.ok || abs(o.w(1) - O2{i-1}.w(1)) > 0.05
        o = slip_orbit(par, vs(i), struct('w', O2{i-1}.w, 'v', vs(i)));
    end
    assert(o.ok && abs(o.w(1) - O2{i-1}.w(1)) < 0.05, 'family 2 lost at v = %.2f', vs(i));
    wprev = O2{i-1}.w;  O2{i} = o;
end

F = struct('vs', vs, 'par', par);
[F.W1, F.rho1, F.Fpk1] = quantities(O1, par);
[F.W2, F.rho2, F.Fpk2] = quantities(O2, par);
end

function o = second_family(par, v)
B = slip_arclength(slip_orbit(par, v, []), 0.02, 700, [-1e-6 0.6]);
k = find([B.turn], 1, 'last');
assert(~isempty(k) && abs(B(k).zeta) < 5e-3, 'fold not passed at v = %.2f', v);
o = slip_orbit(par, v, struct('w', B(k).w, 'v', v));
assert(o.ok && abs(o.w(1) - B(k).w(1)) < 0.05 && o.par(6) == 0, 'family 2 lost in the pin at v = %.2f', v);
end

function [W, rho, Fpk] = quantities(O, par)
n = numel(O);  W = zeros(3,n);  rho = zeros(1,n);  Fpk = zeros(1,n);
for i = 1:n
    W(:,i) = O{i}.w;  Fpk(i) = O{i}.info.Fpk;  rho(i) = max(abs(eig(slip_lin(O{i}, par))));
end
end
