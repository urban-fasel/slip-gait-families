%% Families of periodic motions of a spring-mass hopper, and what a controller buys on each
%
%   A spring-mass hopper with a freely swinging leg (massless-leg SLIP with swing-leg dynamics; parameter
%   values of Raff & Remy, Nonlinear Dynamics 113, 8221-8247, 2025) carries three natural motions on one
%   design: two families of gaits, and the period-2 gait the self-stable family hands over to.  Section 2
%   computes the three, plots their error growth per cycle against speed, and animates them at one speed.
%   Section 3 runs the controllers and the two robustness tests on each family and prints the results.
%
%   Urban Fasel, 2026.
clear; close all; clc
addpath(fullfile(fileparts(mfilename('fullpath')), 'src'))

speeds = 2.0;      % speeds of the controller comparison [m/s]
nMC    = 60;       % runs in the dropout test; 6 for a quick look (the pothole column is deterministic)
vanim  = 3.0;      % animate the three motions at this speed and write gaits.gif; [] to skip

%% 1. The hopper
m = 1;  l0 = 1;  g = 9.81;
par = [m, 40*m*g/l0, l0, g, sqrt(5*g/l0), 0];        % [m k l0 g omega zeta_s]
S.umax   = 15;                                       % swing input limit [rad/s^2]
S.afric  = 0.5;                                      % friction cone: |leg angle| at touchdown [rad]
S.lfloor = 0.6*l0;                                   % leg stroke limit [m]
S.Q = diag([5 1 1]);  S.R = 0.05;                    % LQR weights on [xd; al; ald] and on u

%% 2. The three motions, and their error growth against speed
vs = 1.0:0.05:3.6;
F = slip_families(par, vs);
P = slip_period2(F, par);
show_families(F, P)
plot_backbone(F, P)
if ~isempty(vanim), slip_animate(F, P, par, vanim, fullfile(fileparts(mfilename('fullpath')), 'gaits.gif')); end

%% 3. What a controller can buy, and what it cannot
fprintf('controller comparison at v = %s m/s, nMC = %d ', mat2str(speeds), nMC);  t0 = tic;
R = slip_control(F, par, S, speeds, nMC);
fprintf(' (%.0f s)\n', toc(t0));
show_control(R)


%% ---------------------------------------------------------------------------------------------------
function show_families(F, P)
fprintf('\n%6s | %-24s | %-24s\n', 'v', 'family 1, leg swinging', 'family 2, leg retracting');
fprintf('%6s | %11s %11s | %11s %11s\n', 'm/s', 'Fpk [N]', 'rho', 'Fpk [N]', 'rho');
for i = 1:numel(F.vs)
    if abs(mod(F.vs(i)*10 + 1e-9, 2)) > 1e-6, continue; end                     % every 0.2 m/s
    fprintf('%6.2f | %11.2f %11.3f | %11.2f %11.3f\n', F.vs(i), F.Fpk1(i), F.rho1(i), F.Fpk2(i), F.rho2(i));
end
fprintf('\nrho < 1: the body holds the motion with nothing switched on.  Family 2 does so up to v = %.2f m/s\n', max(F.vs(F.rho2 < 1)));
if ~isempty(P.v)
    j = P.rho < 1;
    fprintf(['and then hands over to a period-2 gait, self-stable up to v = %.2f m/s (rho per hop %.2f to %.2f); v is the\n' ...
             'speed of the family-2 gait of equal energy, the mean running speed of the period-2 gait is %.2f to %.2f m/s.\n'], ...
            max(P.vnom(j)), min(P.rho(j)), max(P.rho(j)), min(P.v(j)), max(P.v(j)));
end
fprintf('\n');
end

function plot_backbone(F, P)
grey = [60 60 60]/255;  blue = [0 78 140]/255;
figure('Color', 'w', 'Position', [100 100 620 380]);  hold on
k = find(F.rho2 < 1, 1, 'last');
h1 = plot(F.vs, F.rho1, '--', 'Color', grey, 'LineWidth', 1.2);
h2 = plot(F.vs(1:k), F.rho2(1:k), '-', 'Color', blue, 'LineWidth', 2);
plot(F.vs(k:end), F.rho2(k:end), '--', 'Color', blue, 'LineWidth', 2);
hs = [h1 h2];  names = {'family 1', 'family 2'};
if ~isempty(P.v)                                 % plotted at the speed of the family-2 gait of equal energy
    j = find(P.rho < 1, 1, 'last');  if isempty(j), j = 1; end
    h3 = plot([F.vs(k) P.vnom(1:j)], [1 P.rho(1:j)], '-', 'Color', blue, 'LineWidth', 1.2);
    plot(P.vnom(j:end), P.rho(j:end), '--', 'Color', blue, 'LineWidth', 1.2);
    plot(F.vs(k), 1, 'o', 'Color', blue, 'MarkerFaceColor', 'w', 'MarkerSize', 5);
    hs(end+1) = h3;  names{end+1} = 'period-2 gait';
end
yline(1, ':k');
set(gca, 'YScale', 'log', 'YTick', [0.5 1 2 5 10], 'YTickLabel', {'0.5' '1' '2' '5' '10'});  grid on;  box on
xlim([F.vs(1) F.vs(end)]);  ylim([0.8 10])
xlabel('running speed [m/s]');  ylabel('error growth per cycle, |\rho|')
title('solid: self-stable;  dashed: not self-stable')
legend(hs, names, 'Location', 'northwest')
end

function show_control(R)
for iv = 1:numel(R.vlist)
    fprintf('\n=== v = %.2f m/s:  rho = %.2f (family 1) / %.2f (family 2);  nominal peak leg force %.1f / %.1f N\n', ...
            R.vlist(iv), R.rho(iv,1), R.rho(iv,2), R.F0(iv,1), R.F0(iv,2));
    fprintf('    %-8s %17s   %17s\n', '', 'family 1', 'family 2');
    fprintf('    %-8s %8s %8s   %8s %8s\n', 'law', 'pothole', 'dropout', 'pothole', 'dropout');
    for j = 1:numel(R.laws)
        fprintf('    %-8s %s   %s\n', R.laws{j}, cells([R.pothole(j,iv,1) R.dropout(j,iv,1)]), ...
                                                 cells([R.pothole(j,iv,2) R.dropout(j,iv,2)]));
    end
end
fprintf(['\n    pothole  deepest pothole cleared in one landing [cm], peak leg force capped at 1.5x nominal\n' ...
         '    dropout  hops survived by every run once the measurement is removed;  >=40: every run was still\n' ...
         '             going when the test stopped at 40 hops\n\n' ...
         '    none = controller off; lqr = the simple controller, the same regulator on both families; nmpc3x = the\n' ...
         '    high-authority controller.  The pothole column moves with the law; the dropout column moves only with\n' ...
         '    the family.\n\n']);
end

function s = cells(v)
s = '';
for x = v
    if isnan(x), s = [s sprintf('%8s', '--')]; elseif x < 0, s = [s sprintf('%8s', sprintf('>=%.0f', -x))];
    else, s = [s sprintf('%8.2f', x)]; end %#ok<AGROW>
end
end
