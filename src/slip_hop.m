function [Zn, info] = slip_hop(Z, U, par, hg, ns)
%SLIP_HOP  One hop of the spring-mass hopper, from one apex to the next.
%   Z: apex states [height; forward speed; leg angle; leg rate], one column per hopper.  U: swing input
%   during the descent [rad/s^2], one per column; it only places the foot and does no work on the body.
%   par = [m k l0 g omega zeta_s].  hg: ground height at this touchdown (default 0).  ns: integration
%   steps in stance (default 300).  Works on many columns at once and with complex-step derivatives.
%   info: ok, al_td (leg angle at touchdown), lmin (shortest leg length), Fpk (peak leg force), yd_lo
%   (vertical speed at lift-off), clear (foot clearance at the apex), Ts, T (stance time, hop period),
%   dist (hop length).
if nargin < 4 || isempty(hg), hg = 0; end,  if nargin < 5, ns = 300; end
m = par(1);  k = par(2);  l0 = par(3);  g = par(4);  w = par(5);  zs = par(6);
N = size(Z,2);  y = Z(1,:) - hg;  xd = Z(2,:);  al0 = Z(3,:);  ald0 = Z(4,:);  c = U./w^2;

% descent: the body falls freely; the leg swings as a damped oscillator driven by u,  al'' = -w^2 al - 2 zs w al' + u
wd = w*sqrt(1 - zs^2);  Bs = (ald0 + zs*w*(al0 - c))/wd;
alf  = @(t) c + exp(-zs*w*t).*((al0 - c).*cos(wd*t) + Bs.*sin(wd*t));
aldf = @(t) exp(-zs*w*t).*((wd*Bs - zs*w*(al0 - c)).*cos(wd*t) - (zs*w*Bs + wd*(al0 - c)).*sin(wd*t));
gtd  = @(t) y - 0.5*g*t.^2 - l0*cos(alf(t));
info.clear = y - l0*cos(al0);  ok = real(info.clear) > 0;
tmax = sqrt(2*max(real(y),0)/g);  tg = linspace(0,1,60).' * tmax;
Gneg = real(y) - 0.5*g*tg.^2 - l0*cos(real(alf(tg))) < 0;
[hit, i1] = max(Gneg, [], 1);  ok = ok & hit & i1 > 1;  i1 = max(i1, 2);
ta = tg(sub2ind(size(tg), i1-1, 1:N));  tb = tg(sub2ind(size(tg), i1, 1:N));
for it = 1:30, tm = (ta+tb)/2;  neg = real(gtd(tm)) < 0;  tb(neg) = tm(neg);  ta(~neg) = tm(~neg); end
t1 = (ta+tb)/2;
for it = 1:3, t1 = t1 - gtd(t1)./(-g*t1 + l0*sin(alf(t1)).*aldf(t1)); end
al = alf(t1);  yd = -g*t1;
S = [al; l0 + 0*al; (xd.*cos(al) - yd.*sin(al))/l0; xd.*sin(al) + yd.*cos(al)];   % stance state [al; l; ald; ld]
info.al_td = al;  ok = ok & real(S(4,:)) < 0 & abs(real(al)) < pi/2;

% stance: integrate until the leg is back at its rest length and extending; Newton refines the end time
f  = @(S) stance(S, m, k, l0, g);
dt = 2.5*pi*sqrt(m/k)/ns;  done = false(1,N);  S0 = S;  jc = ones(1,N);  l1 = S(2,:);
lmin = l0 + 0*al;  gotmin = false(1,N);  lprev = S(2,:);  lcur = S(2,:);
for j = 1:ns
    Sn = rk4(f, S, dt);
    ismin = ~gotmin & j > 1 & real(lcur) <= real(lprev) & real(lcur) < real(Sn(2,:));
    lmin(ismin) = lcur(ismin) - (Sn(2,ismin) - lprev(ismin)).^2 ./ (8*(Sn(2,ismin) - 2*lcur(ismin) + lprev(ismin)));
    gotmin = gotmin | ismin;  lprev = lcur;  lcur = Sn(2,:);
    cross = ~done & real(Sn(2,:)) >= l0 & real(Sn(4,:)) > 0;
    S0(:,cross) = S(:,cross);  l1(cross) = Sn(2,cross);  jc(cross) = j;  done = done | cross;
    S = Sn;  S(:, ~isfinite(real(S(2,:))) | real(S(2,:)) < 0.05*l0) = NaN;
    if all(done | ~ok | isnan(real(S(2,:)))), break; end
end
ok = ok & done & gotmin;
tau = (l0 - S0(2,:))./(l1 - S0(2,:));  tau(~done) = 0.5;
for it = 1:4, Sl = rk4(f, S0, tau*dt);  tau = tau - (Sl(2,:) - l0)./(Sl(4,:)*dt); end
Sl = rk4(f, S0, tau*dt);  info.Ts = (jc - 1 + tau)*dt;  info.lmin = lmin;  info.Fpk = k*(l0 - lmin);

% lift-off and flight to the next apex
a2 = Sl(1,:);  l = Sl(2,:);  ald = Sl(3,:);  ld = Sl(4,:);
yl = l.*cos(a2);  xd1 = ld.*sin(a2) + l.*ald.*cos(a2);  yd1 = ld.*cos(a2) - l.*ald.*sin(a2);
info.yd_lo = yd1;  ok = ok & real(yd1) > 0 & real(yl) > 0;  t3 = yd1/g;
B3 = (ald + zs*w*a2)/wd;  e3 = exp(-zs*w*t3);
Zn = [hg + yl + yd1.^2/(2*g);  xd1;  e3.*(a2.*cos(wd*t3) + B3.*sin(wd*t3)); ...
      e3.*((wd*B3 - zs*w*a2).*cos(wd*t3) - (zs*w*B3 + wd*a2).*sin(wd*t3))];
Zn(:,~ok) = NaN;  info.ok = ok;  info.T = t1 + info.Ts + t3;
info.dist = xd.*t1 + (l.*sin(a2) - l0*sin(al)) + xd1.*t3;
end

function dS = stance(S, m, k, l0, g)
% stance dynamics: a spring-loaded inverted pendulum about the foot, S = [al; l; ald; ld]
al = S(1,:);  l = S(2,:);  ald = S(3,:);  ld = S(4,:);
dS = [ald;  ld;  -2*ld.*ald./l + (g./l).*sin(al);  l.*ald.^2 - g*cos(al) - (k/m)*(l - l0)];
end

function Sn = rk4(f, S, h)
k1 = f(S);  k2 = f(S + 0.5*h.*k1);  k3 = f(S + 0.5*h.*k2);  k4 = f(S + h.*k3);
Sn = S + (h/6).*(k1 + 2*k2 + 2*k3 + k4);
end
