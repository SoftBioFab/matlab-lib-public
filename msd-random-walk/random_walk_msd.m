%% 2D random walk and mean squared displacement (MSD)
% Each step: dx, dy drawn uniformly from [-1, 1].
% Expected MSD for this walk: MSD(lag) = 2 * (1/3) * lag = (2/3) * lag

clear; close all; clc;

% Folder for saved figures (next to this script)
scriptDir = fileparts(mfilename('fullpath'));
repoDir   = fileparts(scriptDir);                      % main repo folder (one level up)
fprintf('Repo folder (for git):\ncd ''%s''\n\n', repoDir);   % ready to copy
figDir    = fullfile(scriptDir, 'figures');
if ~exist(figDir, 'dir'); mkdir(figDir); end
rng(1);                       % fixed seed for reproducibility

%% Parameters
nTraj  = 100;                 % number of trajectories
nSteps = 1000;                % number of points per trajectory
maxLag = 250;                 % largest lag time for the MSD (in steps)

%% Generate random walks
% Steps: rows = time, columns = trajectories
dx = 2*rand(nSteps-1, nTraj) - 1;
dy = 2*rand(nSteps-1, nTraj) - 1;

% Positions, starting at the origin
x = [zeros(1, nTraj); cumsum(dx, 1)];
y = [zeros(1, nTraj); cumsum(dy, 1)];

%% Plot a few trajectories
nShow = 5;
figTraj = figure;
hold on;
for k = 1:nShow
    plot(x(:,k), y(:,k), 'LineWidth', 1, ...
         'DisplayName', sprintf('Traj %d', k));
end
plot(0, 0, 'k.', 'MarkerSize', 20, ...
     'HandleVisibility', 'off');      % common starting point (not in legend)
hold off;
legend('Location', 'best');
axis equal; box on;
xlabel('x'); ylabel('y');
title(sprintf('%d example trajectories', nShow));
exportgraphics(figTraj, fullfile(figDir, 'trajectories.png'), 'Resolution', 300);

%% Calculate MSD
% 1) Time-averaged MSD for each trajectory (average over time origins)
% 2) Mean over trajectories
% 3) Error = standard error of the mean over trajectories (std/sqrt(nTraj)).
%    Trajectories are independent, so this is a valid error estimate.
%    (The points within one trajectory are NOT independent: overlapping
%    windows are correlated, so 1/sqrt(number of points) underestimates it.)
lags    = (1:maxLag)';
msdTraj = zeros(maxLag, nTraj);      % rows = lag, columns = trajectory
for i = 1:maxLag
    tau = lags(i);
    ddx = x(1+tau:end, :) - x(1:end-tau, :);
    ddy = y(1+tau:end, :) - y(1:end-tau, :);
    msdTraj(i, :) = mean(ddx.^2 + ddy.^2, 1);
end
msd    = mean(msdTraj, 2);                 % ensemble-averaged MSD
msdStd = std(msdTraj, 0, 2);               % spread between trajectories
msdSem = msdStd / sqrt(nTraj);             % error on the mean MSD (standard uncertainty)
k      = 2;                                % coverage factor (~95%)
msdU   = k * msdSem;                       % expanded uncertainty, used for display

%% Fit and plot MSD
% Weighted linear fit: MSD = slope*lag + offset, weights = 1/SEM^2
% MSD values at different lags are correlated, so a single fit to the mean
% curve underestimates the error. We fit each trajectory with the same
% weights and take the mean and SEM of the parameters over trajectories,
% which are independent.
A  = [lags, ones(maxLag, 1)];
sw = 1 ./ msdSem;                          % sqrt of weights 1/SEM^2
P  = (A .* sw) \ (msdTraj .* sw);          % 2 x nTraj: [slope; offset] per trajectory
p    = mean(P, 2);                         % mean slope and offset
pErr = k * std(P, 0, 2) / sqrt(nTraj);     % expanded uncertainty (k*SEM) of slope and offset

figMsd = figure;
idx = 1:10:maxLag;                         % show every 10th lag to avoid clutter
errorbar(lags(idx), msd(idx), msdU(idx), 'o', 'MarkerSize', 4); hold on;
plot(lags, (2/3)*lags, 'k--', 'LineWidth', 1.5);
plot(lags, A*p, 'r-', 'LineWidth', 1);
hold off; box on;
xlabel('Lag (steps)'); ylabel('MSD');
legend(sprintf('Simulation (mean \\pm %d SEM)', k), 'Theory: (2/3) \tau', ...
       sprintf('Weighted fit: slope = %s (k=%d)', fmtErr(p(1), pErr(1)), k), ...
       'Location', 'northwest');
title('Mean squared displacement');
exportgraphics(figMsd, fullfile(figDir, 'msd.png'), 'Resolution', 300);

fprintf('Errors are expanded uncertainties, k = %d (~95%%)\n', k);
fprintf('Fitted slope:  %s (theory: 0.667)\n', fmtErr(p(1), pErr(1)));
fprintf('Fitted offset: %s (theory: 0)\n', fmtErr(p(2), pErr(2)));

% Diffusion coefficient in 2D: MSD = 4*D*lag
D = p(1) / 4;  DErr = pErr(1) / 4;
fprintf('Diffusion coefficient D: %s (theory: 0.167)\n', fmtErr(D, DErr));

%% Helper: round error to 1 significant digit, value to the same decimal
function s = fmtErr(val, err)
    n = -floor(log10(err));             % decimal position of first significant digit
    err = round(err, n);
    n = -floor(log10(err));             % recompute in case rounding gave e.g. 0.01
    nd = max(n, 0);                     % number of decimals to show
    s = sprintf('%.*f %s %.*f', nd, round(val, n), char(177), nd, err);   % char(177) = ±
end

