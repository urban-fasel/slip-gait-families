function slip_animate(F, P, par, v, file)
%SLIP_ANIMATE  Animate the natural motions at speed v: family 1, family 2 and, if it exists there, the
%   period-2 gait at the energy of family 2.  Five hops each, played at half speed, apices marked.
%   file: GIF to write; without one the animation plays on screen.
m = par(1);  g = par(4);  grey = [60 60 60]/255;  blue = [0 78 140]/255;
[~, i] = min(abs(F.vs - v));
o1 = slip_orbit(par, v, struct('w', F.W1(:,i), 'v', v));  rho1 = max(abs(eig(slip_lin(o1, par))));
o2 = slip_orbit(par, v, struct('w', F.W2(:,i), 'v', v));  rho2 = max(abs(eig(slip_lin(o2, par))));
held = {'not self-stable', 'self-stable'};
G = { sprintf('family 1, leg swinging forward at touchdown.   v = %.1f m/s,   \\rho = %.2f:  %s', v, rho1, held{1 + (rho1 < 1)}), o1.z, grey
      sprintf('family 2, leg retracting at touchdown.   v = %.1f m/s,   \\rho = %.2f:  %s', v, rho2, held{1 + (rho2 < 1)}), o2.z, blue };
j = find(abs(P.vnom - v) < 1e-9, 1);
if ~isempty(j)
    G(3,:) = { sprintf('period-2 gait at the energy of family 2.   mean speed %.2f m/s,   \\rho = %.2f per hop:  %s', ...
                       P.v(j), P.rho(j), held{1 + (P.rho(j) < 1)}), P.za(:,j), blue };
end
ng = size(G,1);  nhop = 5;  dt = 0.004;  T = cell(1,ng);
for k = 1:ng, T{k} = trajectory(G{k,2}, par, nhop, dt); end
tend = min(cellfun(@(T) T.t(end), T));
xmax = max(cellfun(@(T) T.xb(find(T.t <= tend, 1, 'last')), T)) + 0.5;
fps = 25;  tf = 0:1/(2*fps):tend;                                    % played at half speed

vis = 'on';  if nargin > 4, vis = 'off'; end
fig = figure('Color', 'w', 'Position', [100 100 1000 240*ng], 'Visible', vis);
ax = gobjects(1,ng);
for k = 1:ng
    ax(k) = subplot(ng,1,k);  hold(ax(k), 'on');  box(ax(k), 'on')
    plot(ax(k), [-1 xmax+1], [0 0], 'k', 'LineWidth', 1.2)
    set(ax(k), 'DataAspectRatio', [1 1 1], 'XLim', [-0.4 xmax], 'YLim', [-0.08 2.1], 'XTick', 0:xmax, 'YTick', 0:0.5:2, 'FontSize', 9)
    text(ax(k), -0.25, 1.92, G{k,1}, 'FontSize', 10)
    ylabel(ax(k), 'y [m]');  if k == ng, xlabel(ax(k), 'x [m]'); end
end
for f = 1:numel(tf)
    for k = 1:ng
        Tk = T{k};  q = find(Tk.t <= tf(f), 1, 'last');  c = G{k,3};
        delete(findobj(ax(k), 'Tag', 'moving'))
        plot(ax(k), Tk.xb(1:q), Tk.yb(1:q), ':', 'Color', c, 'LineWidth', 0.8, 'Tag', 'moving')
        plot(ax(k), Tk.xa(Tk.ta <= tf(f)), Tk.ya(Tk.ta <= tf(f)), 'v', 'Color', c, 'MarkerSize', 4, 'Tag', 'moving')
        plot(ax(k), [Tk.xf(q) Tk.xb(q)], [Tk.yf(q) Tk.yb(q)], '-', 'Color', c, 'LineWidth', 2.5, 'Tag', 'moving')
        plot(ax(k), Tk.xb(q), Tk.yb(q), 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'MarkerSize', 11, 'Tag', 'moving')
        plot(ax(k), Tk.xf(q), Tk.yf(q), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', c, 'MarkerSize', 4, 'Tag', 'moving')
    end
    text(ax(1), xmax - 0.1, 1.92, sprintf('t = %.2f s', tf(f)), 'HorizontalAlignment', 'right', 'FontSize', 9, 'Tag', 'moving');
    if nargin > 4
        [im, map] = rgb2ind(print(fig, '-RGBImage', '-r80'), 128);
        if f == 1, imwrite(im, map, file, 'gif', 'LoopCount', inf, 'DelayTime', 1/fps);
        else,      imwrite(im, map, file, 'gif', 'WriteMode', 'append', 'DelayTime', 1/fps); end
    else
        drawnow;  pause(1/fps)
    end
end
if nargin > 4, close(fig);  fprintf('wrote %s\n', file); end
end

function T = trajectory(z0, par, nhop, dt)
% body (xb, yb), foot (xf, yf) and apex points over nhop hops from the apex state z0
m = par(1);  k = par(2);  l0 = par(3);  g = par(4);  w = par(5);  hs = 2.5*pi*sqrt(m/k)/300;
T = struct('t', [], 'xb', [], 'yb', [], 'xf', [], 'yf', [], 'ta', [], 'xa', [], 'ya', []);
tnow = 0;  xnow = 0;  Z = z0;
for hop = 1:nhop
    y0 = Z(1);  xd = Z(2);  al0 = Z(3);  ald0 = Z(4);
    [Zn, info] = slip_hop(Z, 0, par);  assert(info.ok, 'hop failed')
    t3 = info.yd_lo/g;  t1 = info.T - info.Ts - t3;
    T.ta(end+1) = tnow;  T.xa(end+1) = xnow;  T.ya(end+1) = y0;
    % descent: the body falls freely, the leg swings about the vertical
    tt = (0:dt:t1).';  al = al0*cos(w*tt) + (ald0/w)*sin(w*tt);  X = xnow + xd*tt;  Y = y0 - 0.5*g*tt.^2;
    T = add(T, tnow + tt, X, Y, X - l0*sin(al), Y - l0*cos(al));
    % stance about the fixed foot, S = [al; l; ald; ld]
    altd = al0*cos(w*t1) + (ald0/w)*sin(w*t1);  yd = -g*t1;
    S = [altd; l0; (xd*cos(altd) - yd*sin(altd))/l0; xd*sin(altd) + yd*cos(altd)];
    xfoot = xnow + xd*t1 - l0*sin(altd);  ts = 0;
    while true
        Sn = rk4(S, hs, m, k, l0, g);
        if Sn(2) >= l0 && Sn(4) > 0, break; end
        S = Sn;  ts = ts + hs;  assert(ts < 5, 'stance did not end')
        T = add(T, tnow + t1 + ts, xfoot + S(2)*sin(S(1)), S(2)*cos(S(1)), xfoot, 0);
    end
    tau = (l0 - S(2))/(Sn(2) - S(2));
    for it = 1:4, Sl = rk4(S, tau*hs, m, k, l0, g);  tau = tau - (Sl(2) - l0)/(Sl(4)*hs); end
    Sl = rk4(S, tau*hs, m, k, l0, g);  ts = ts + tau*hs;
    % flight to the next apex
    a2 = Sl(1);  l = Sl(2);  ald = Sl(3);  ld = Sl(4);
    xl = xfoot + l*sin(a2);  yl = l*cos(a2);  xd1 = ld*sin(a2) + l*ald*cos(a2);  yd1 = ld*cos(a2) - l*ald*sin(a2);
    tt = (0:dt:yd1/g).';  al = a2*cos(w*tt) + (ald/w)*sin(w*tt);  X = xl + xd1*tt;  Y = yl + yd1*tt - 0.5*g*tt.^2;
    T = add(T, tnow + t1 + ts + tt, X, Y, X - l0*sin(al), Y - l0*cos(al));
    tnow = tnow + info.T;  xnow = xnow + info.dist;  Z = Zn;
end
end

function T = add(T, t, xb, yb, xf, yf)
T.t = [T.t; t(:)];  T.xb = [T.xb; xb(:)];  T.yb = [T.yb; yb(:)];  T.xf = [T.xf; xf(:)];  T.yf = [T.yf; yf(:)];
end

function Sn = rk4(S, h, m, k, l0, g)
f = @(S) [S(3); S(4); -2*S(4)*S(3)/S(2) + (g/S(2))*sin(S(1)); S(2)*S(3)^2 - g*cos(S(1)) - (k/m)*(S(2) - l0)];
k1 = f(S);  k2 = f(S + 0.5*h*k1);  k3 = f(S + 0.5*h*k2);  k4 = f(S + h*k3);
Sn = S + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
end
