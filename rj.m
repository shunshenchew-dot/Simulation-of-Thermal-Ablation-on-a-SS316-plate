clear; clc; close all;

% Material Properties of SS316 
mat.rho = 7900;                 % kg/m^3
mat.cp  = 470;                  % J/(kg*K)
mat.k   = 15;                   % W/(m*K)
mat.a_th = mat.k/(mat.rho*mat.cp);

mat.T0 = 293.15;                % K
mat.Tm = 1663.15;               % K
mat.Tv = 3073.15;               % K

mat.Lm = 290e3;                 % J/kg
mat.Lv = 6090e3;                % J/kg

% Minimum entalphy for material removal
mat.h_remove_mass = mat.cp*(mat.Tm-mat.T0) + mat.Lm + mat.cp*(mat.Tv-mat.Tm) + mat.Lv; % J/kg
mat.H_remove      = mat.rho * mat.h_remove_mass; % J/m^3

fprintf('\n=== MINIMUM ENTHALPY FOR MATERIAL REMOVAL ===\n');
fprintf('h_remove (J/kg)  = %.6e\n', mat.h_remove_mass);
fprintf('H_remove (J/m^3) = %.6e\n\n', mat.H_remove);

%% ---------------- DAMAGE / MICROSTRUCTURE PROXIES (VISUAL MODEL) ----------------
mat.T_HAZ      = 0.60 * mat.Tm;     % K (DamageMap HAZ)
mat.T_damage   = 0.85 * mat.Tm;     % K
mat.GradT_crit = 2.0e9;             % K/m (tune 1e9–5e9)

% HAZ outline threshold for temperature plots:
mat.T_HAZ_outline = 773.15;         % K

% Fixed Tmin for temperature plots:
mat.Tmin_plot = 293.15;             % K

%Laser Settings
laser.P_av   = 30;              % W
laser.f      = 102e3;           % Hz
laser.tau    = 140e-9;          % s
laser.d_spot = 40e-6;           % m  (1/e^2 diameter)
laser.t_peak = 0.5*laser.tau;   % peak at 50% of tau

laser.eta_abs   = 0.0170;
laser.alpha_abs = 5.3e7;        % 1/m

% Geometry definition
dom.W  = 400e-6;                % m
dom.H  = 80e-6;                 % m
dom.x0 = dom.W/2;

%% ---------------- SNAPSHOTS ----------------
snap_fracs = 0:0.1:1.0;

% ---------------- ZOOM WINDOW (snapshot close-ups) ----------------
zoom.x_halfwidth_um = 50;   % ±50 µm around crater centre
zoom.z_max_um       = 15;   % top 15 µm for temp/phase zoom

%% ---------------- BASELINE DISCRETISATION ----------------
base.dx = 2e-6;
base.dz = 0.4e-6;
base.cfl_safety = 0.2;
base.forceNtDivisibleBy10 = true;

%% ---------------- PROFILOMETER BASELINE ----------------
surface_offset = 688.2; % um (arbitrary baseline)

%% =====================================================================
% (0) BASELINE SINGLE CASE
%% =====================================================================
opts0 = struct();
opts0.doSnapshots = true;
opts0.doAllPlots  = true;
opts0.caseName    = 'BASELINE_SINGLE_FULL';

res0 = run_single_case_STABLE_FULL(mat, laser, dom, base.dx, base.dz, base.cfl_safety, ...
    base.forceNtDivisibleBy10, snap_fracs, opts0, surface_offset, zoom);

fprintf('\nBASELINE SINGLE | Max depth = %.3f um\n', res0.maxDepth_um);
fprintf('H_remove per cell (J) = %.6e (cellVol=dx*dz*1m)\n', res0.Hremove_cell_J);

fprintf('\n=== SNAPSHOT ENTHALPY (TOTAL ABOVE T0, incl latent history) ===\n');
for k = 1:numel(res0.snap_pct)
    fprintf('Snapshot %3d%% : H = %.6e J (per 1 um width)\n', res0.snap_pct(k), res0.Hsnap_J(k));
end

% Enthalpy vs snapshot + ablation threshold line (SAFE)
figure('Name','Enthalpy vs Snapshot','Color','w','Position',[120 120 980 520]);
plot(res0.snap_pct, res0.Hsnap_J, '-o', 'LineWidth', 2);
grid on; hold on;

N_ab_cols = sum(res0.surf_k > 1);
E_threshold_total = res0.Hremove_cell_J * N_ab_cols;

safe_hline(E_threshold_total, ':', [1 0 0], 2, 'Ablation threshold');
xlabel('Snapshot (%)');
ylabel('Total enthalpy above T_0 (J per 1 um width)');
title('Total domain enthalpy vs snapshot with ablation threshold');
legend({'Total domain enthalpy','Ablation threshold'}, 'Location','best');
hold off;

%% =====================================================================
% (0b) MESH VISUALISATION (BASELINE)
%% =====================================================================
plot_simulation_mesh(dom, res0.dx, res0.dz, zoom); drawnow;

%% =====================================================================
% PAPER-STYLE FIGURES
%   - Top-view reconstruction uses instantaneous surface T (Tfinal)
%   - Distance decay uses Tmax
%   - Depth plot uses INSTANTANEOUS (Tfinal) at final time
%% =====================================================================
plot_topview_Tinst_with_crater(res0, dom, mat);
plot_Tmax_vs_distance_from_crater(res0, dom);
plot_centerline_Tinst_vs_depth_FINAL(res0, dom, mat);

%% =====================================================================
% RELATIVE TEMPORAL INTENSITY (peaks at 0.5*tau)
%% =====================================================================
plot_relative_temporal_intensity(laser);

%% =====================================================================
% CRATER DEPTH (0–2.0 um)
%% =====================================================================
figure('Name','Simulated Crater Profile','Color','w','Position',[100,100,1200,520]);
plot(res0.x_um, res0.depth_um, 'b-', 'LineWidth', 2);
grid on; grid minor;
xlabel('X (\mum)'); ylabel('Depth (\mum)');
title('Simulated Crater Profile');
xlim([0 400]); ylim([0 2.0]);
set(gca,'YDir','reverse');
pbaspect([400 (2.0*6) 1]);

%% =====================================================================
% PROFILOMETER STYLE PLOT (FINAL) + RECAST THICKNESS OVERLAY
%% =====================================================================
z_profile = surface_offset - res0.depth_um;
x_profile = res0.x_um;

figure('Name', 'Simulated Profilometer Scan (Final)', 'Color', 'w', 'Position', [100, 100, 1200, 340]);
plot(x_profile, z_profile, 'b-', 'LineWidth', 2);
hold on; grid on; grid minor;
ax = gca; ax.GridAlpha = 0.4; ax.MinorGridAlpha = 0.2; ax.FontSize = 10; ax.LineWidth = 1;
axis tight;
ylim([surface_offset-2.0, surface_offset+0.2]);
xlabel('X (\mum)', 'FontWeight', 'bold'); ylabel('Z (\mum)', 'FontWeight', 'bold');
title('Simulated Crater Profile (Profilometer View, Final)');

x_range = max(x_profile) - min(x_profile);
y_range = max(z_profile) - min(z_profile);
if y_range == 0, y_range = 1; end
pbaspect([x_range (y_range*6) 1]);

ar_val = x_range / y_range;
text(max(x_profile), min(z_profile)+0.02, sprintf('simulated aspect ratio: 1:%.2f', ar_val), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
    'Color', 'r', 'FontWeight', 'bold', 'FontSize', 10);
hold off;

% Recast thickness along X from RecastMask
dz_um = res0.dz * 1e6;
Nx = res0.Nx;
recast_thick_um = zeros(1, Nx);
surf_k_from_depth = round(res0.depth_um / dz_um) + 1;

for j = 1:Nx
    ks = max(1, min(res0.Nz, surf_k_from_depth(j)));
    k = ks; cnt = 0;
    while (k <= res0.Nz) && res0.RecastMask(k,j)
        cnt = cnt + 1; k = k + 1;
    end
    recast_thick_um(j) = cnt * dz_um;
end

figure('Name','Recast thickness along X (overlay)','Color','w','Position',[100,100,1200,380]);
yyaxis left
plot(res0.x_um, surface_offset - res0.depth_um, 'LineWidth', 2); grid on; hold on;
ylabel('Surface height Z (\mum)');
xlabel('X (\mum)');
title('Profilometer surface + Recast thickness (RecastMask = Melted & ~Removed)');
ylim([surface_offset-2.0, surface_offset+0.2]);
yyaxis right
plot(res0.x_um, recast_thick_um, 'LineWidth', 2);
ylabel('Recast thickness (\mum)');
legend({'Surface profile','Recast thickness'}, 'Location','best');

fprintf('Any recast? %d | Max Tmax = %.0f K\n', any(res0.RecastMask(:)), max(res0.Tmax(:)));

%% =====================================================================
% NEW: CUMULATIVE MATERIAL REMOVAL RATE (V/tau ONLY)
%% =====================================================================
plot_cumulative_MRR_VoverTau(res0, laser);

%% =====================================================================
% (1) BATCH RUN: Settings 1–10 + MAE/MAPE + MRR (V/tau ONLY)
%% =====================================================================
expCases = [
    1  220  70   2.28
    2  200  78   2.42
    3  190  88   2.04
    4  160  96   1.88
    5  140  102  2.04
    6  120  110  1.18
    7  100  120  0.88
    8   85  136  0.67
    9   75  148  0.78
    10  55  180  0.44
];
Ncases = size(expCases,1);

SimMaxDepth_um      = zeros(Ncases,1);
AbsErr_um           = zeros(Ncases,1);

Vsim_m3_per_mwidth  = zeros(Ncases,1);
MRR_VoverTau        = zeros(Ncases,1);

optsBatch = struct('doSnapshots',false,'doAllPlots',false,'caseName','BATCH_BASE_FULL');

fprintf('\n===================== BATCH RUN (BASELINE PARAMS) =====================\n');
for i = 1:Ncases
    setting = expCases(i,1);
    tau_ns  = expCases(i,2);
    f_kHz   = expCases(i,3);
    exp_um  = expCases(i,4);

    laser_i = laser;
    laser_i.tau    = tau_ns*1e-9;
    laser_i.f      = f_kHz*1e3;
    laser_i.t_peak = 0.5*laser_i.tau;

    optsBatch.caseName = sprintf('SET_%02d_BASE', setting);

    res = run_single_case_STABLE_FULL(mat, laser_i, dom, base.dx, base.dz, base.cfl_safety, ...
        base.forceNtDivisibleBy10, snap_fracs, optsBatch, surface_offset, zoom);

    SimMaxDepth_um(i) = res.maxDepth_um;
    AbsErr_um(i)      = abs(SimMaxDepth_um(i) - exp_um);

    [Vsim_m3_per_mwidth(i), MRR_VoverTau(i)] = ...
    compute_volume_and_MRR(res.x_um, res.depth_um, laser_i.f);

    fprintf('Setting %02d | Exp=%.3f | Sim=%.3f | AbsErr=%.3f | V_pulse=%.3e m^3 | MRR(V*f)=%.3e m^3/s\n', ...
        setting, exp_um, SimMaxDepth_um(i), AbsErr_um(i), Vsim_m3_per_mwidth(i), MRR_VoverTau(i));
end

OverallMAE_um = mean(AbsErr_um);
BaselineMAPE  = mean(AbsErr_um ./ expCases(:,4))*100;
fprintf('\nBaseline MAE  = %.4f um\n', OverallMAE_um);
fprintf('Baseline MAPE = %.2f %%\n\n', BaselineMAPE);

figure('Name','Experimental vs Simulated (Baseline)','NumberTitle','off'); set(gcf,'Color','w');
plot(expCases(:,1), expCases(:,4), '-o', 'LineWidth', 2); hold on; grid on;
plot(expCases(:,1), SimMaxDepth_um, '-o', 'LineWidth', 2);
xlabel('Setting #'); ylabel('Max depth [\mum]');
title(sprintf('Exp vs Sim (baseline) | MAPE = %.2f%%', BaselineMAPE));
legend({'Experimental','Sim (baseline)'}, 'Location','best');

figure('Name','Abs Error per Setting (Baseline)','NumberTitle','off'); set(gcf,'Color','w');
bar(expCases(:,1), AbsErr_um); grid on;
xlabel('Setting #'); ylabel('Abs error [\mum]');
title(sprintf('Baseline abs error (MAE = %.4f \\mum)', OverallMAE_um));

figure('Name','MRR (V*f) per Setting','NumberTitle','off'); set(gcf,'Color','w');
plot(expCases(:,1), MRR_VoverTau, '-o', 'LineWidth', 2); grid on;
xlabel('Setting #');
ylabel('MRR_{equiv} = V_{pulse} \cdot f (m^3/s per 1 um width)');
title('Equivalent MRR assuming repetition rate f');

%% =====================================================================
% (2) CONVERGENCE STUDY — TRUE 2D MESH REFINEMENT (dx AND dz)
%% =====================================================================
dx_list_um = [10 8 6 4 3 2];
dx_list = dx_list_um * 1e-6;
cfl_fixed = base.cfl_safety;

optsSweep = struct('doSnapshots',false,'doAllPlots',false,'caseName','MESH_SWEEP_FULL');

dxCell = cell(1,numel(dx_list));
for i = 1:numel(dx_list)
    dx_i = dx_list(i);
    dz_i = dx_i/5; % keep ratio
    optsSweep.caseName = sprintf('MESH_dx_%dum', round(dx_i*1e6));

    dxCell{i} = run_single_case_STABLE_FULL(mat, laser, dom, dx_i, dz_i, cfl_fixed, ...
        false, snap_fracs, optsSweep, surface_offset, zoom);
end
dxResults = [dxCell{:}];

[~, idxMinDx] = min([dxResults.dx]);
refDX = dxResults(idxMinDx);

for i = 1:numel(dxResults)
    dxResults(i).L2err_vs_ref = L2_on_common_support(dxResults(i).x_um, dxResults(i).depth_um, ...
                                                    refDX.x_um, refDX.depth_um);
end

dx_um_list  = [dxResults.dx]'*1e6;
dz_um_list  = [dxResults.dz]'*1e6;
Nz_list     = [dxResults.Nz]';
Nx_list     = [dxResults.Nx]';
Nodes_list  = [dxResults.Nodes]';
Depth_list  = [dxResults.maxDepth_um]';
CPU_list    = [dxResults.cpu_s]';
L2_list     = [dxResults.L2err_vs_ref]';

Tmesh_dx = table(dx_um_list, dz_um_list, Nx_list, Nz_list, Nodes_list, Depth_list, L2_list, CPU_list, ...
    'VariableNames', {'dx_um','dz_um','Nx','Nz','Nodes','MaxDepth_um','L2err_um','CPU_s'});
disp('--- TRUE 2D Mesh Sweep Convergence Table ---');
disp(Tmesh_dx);

plot_mesh_family_one_window(dom, dx_list, zoom);

% DT sweep via CFL list
cfl_list = [0.45, 0.30, 0.20, 0.15, 0.10];
optsSweep.caseName = 'DT_SWEEP_FULL';

dtCell = cell(1,numel(cfl_list));
for i = 1:numel(cfl_list)
    dtCell{i} = run_single_case_STABLE_FULL(mat, laser, dom, base.dx, base.dz, cfl_list(i), false, snap_fracs, optsSweep, surface_offset, zoom);
end
dtResults = [dtCell{:}];

[~, idxMinDt] = min([dtResults.dt]);
refDT = dtResults(idxMinDt);
for i = 1:numel(dtResults)
    dtResults(i).L2err_vs_ref = L2_on_common_support(dtResults(i).x_um, dtResults(i).depth_um, ...
                                                    refDT.x_um, refDT.depth_um);
end

% Solve-time vs nodes (from dx sweep)
Nodes = [dxResults.Nodes];
Depth = [dxResults.maxDepth_um];
TimeS = [dxResults.cpu_s];
[NodesS, idxN] = sort(Nodes,'ascend');

figure('Name','Mesh Convergence Analysis','NumberTitle','off'); set(gcf,'Color','w');
yyaxis left
plot(NodesS, Depth(idxN), '-o', 'LineWidth', 2); grid on; hold on;
ylabel('Max Depth (\mum)');
yyaxis right
plot(NodesS, TimeS(idxN), '-o', 'LineWidth', 2);
ylabel('Solve Time (s)');
set(gca,'XScale','log');
xlabel('Number of Nodes (Nx*Nz)');
title('Mesh Convergence Analysis');
legend({'Depth','Solve time'}, 'Location','best');

%% =====================================================================
% (3) HYBRID RUN (OPT only for settings 6–10)
%% =====================================================================
opt.eta_abs   = 0.130912;
opt.alpha_abs = 1.082956e7;
opt.kH        = 0.830951;

SimHybrid_um   = zeros(Ncases,1);
AbsErr_hyb_um  = zeros(Ncases,1);

optsHyb = struct('doSnapshots',false,'doAllPlots',false,'caseName','HYBRID_FULL');

fprintf('\n===================== HYBRID RUN (OPT APPLIED TO SETTINGS 6–10 ONLY) =====================\n');
% ==== PREALLOCATE (avoid growth + avoid old values leaking in) ====
SimSimulation_um = nan(Ncases,1);
AbsErr_sim_um    = nan(Ncases,1);

% ==== RUN CASES ====
for i = 1:Ncases
    setting = expCases(i,1);
    tau_ns  = expCases(i,2);
    f_kHz   = expCases(i,3);
    exp_um  = expCases(i,4);

    laser_i = laser;
    mat_i   = mat;

    % Tag for naming only (BASE vs OPT)
    tag = 'BASE';
    if setting >= 6 && setting <= 10
        laser_i.eta_abs   = opt.eta_abs * (4*log(2)/8);
        laser_i.alpha_abs = opt.alpha_abs;
        mat_i.H_remove    = opt.kH * mat.H_remove;
        tag = 'OPT';
    end

    % Set pulse params for this case
    laser_i.tau    = tau_ns * 1e-9;
    laser_i.f      = f_kHz  * 1e3;
    laser_i.t_peak = 0.5 * laser_i.tau;

    % Case name (pure label)
    optsHyb.caseName = sprintf('SET_%02d_%s', setting, tag);

    % Run
    res = run_single_case_STABLE_FULL(mat_i, laser_i, dom, ...
        base.dx, base.dz, base.cfl_safety, ...
        base.forceNtDivisibleBy10, snap_fracs, ...
        optsHyb, surface_offset, zoom);

    % Store results
    SimSimulation_um(i) = res.maxDepth_um;
    AbsErr_sim_um(i)    = abs(SimSimulation_um(i) - exp_um);

    fprintf('Setting %02d [%s] | Exp=%.3f | Sim=%.3f | AbsErr=%.3f\n', ...
        setting, tag, exp_um, SimSimulation_um(i), AbsErr_sim_um(i));
end

% ==== ERROR METRICS (ONCE ONLY) ====
exp_all_um = expCases(:,4);

MAE_sim_um = mean(abs(SimSimulation_um - exp_all_um), 'omitnan');
MAPE_sim   = mean(abs(SimSimulation_um - exp_all_um) ./ exp_all_um, 'omitnan') * 100;

fprintf('\nSimulation MAE  = %.4f um\n', MAE_sim_um);
fprintf('Simulation MAPE = %.2f %%\n\n', MAPE_sim);

% ==== PLOT (ONCE ONLY) ====
figure('Name','Exp vs Simulation','NumberTitle','off');
set(gcf,'Color','w');

plot(expCases(:,1), exp_all_um, '-o', 'LineWidth', 2); hold on; grid on;
plot(expCases(:,1), SimSimulation_um, '-o', 'LineWidth', 2);

xlabel('Setting #');
ylabel('Max depth [\mum]');
title(sprintf('Experimental vs Simulation | MAE = %.3f \\mum | MAPE = %.2f%%', ...
    MAE_sim_um, MAPE_sim));

legend({'Experimental','Simulation'}, 'Location','best');



%% =====================================================================
% (4) BEAM PLOTS (baseline laser) — 2D + 3D
%% =====================================================================
plot_beam_2d_3d(laser);

%% =====================================================================
% ===================== LOCAL FUNCTIONS ================================
%% =====================================================================

function res = run_single_case_STABLE_FULL(mat, laser, dom, dx, dz, cfl_safety, forceNtDivisibleBy10, snap_fracs, opts, surface_offset, zoom)
    tStart = tic;

    % Grid
    W  = dom.W;  H  = dom.H;
    Nx = round(W/dx) + 1;
    Nz = round(H/dz) + 1;
    x = linspace(0, W, Nx);
    z = linspace(0, H, Nz);

    % Snap x0 to nearest node
    [~, j0] = min(abs(x - dom.x0));
    x0 = x(j0);

    % Laser normalisation (1/e^2)
    E_pulse  = laser.P_av / laser.f;
    time_int = laser.tau * sqrt(pi/(4*log(2)));
    P_peak   = E_pulse / time_int;
    d = laser.d_spot;
    I0 = (8 / (pi*d^2)) * P_peak;

    % x profile
    px_col = gaussian_cell_average_x(x, x0, dx, d);
    r_cut = 40e-6;           % 40 µm radius
    r = abs(x - x0);
    taper = 0.5*(1 - tanh((r - r_cut)/(5e-6)));  % 5 µm smooth edge
    px_col = px_col .* taper;


    if max(px_col) > 0
        px_col = px_col / max(px_col);
    end
    px = repmat(px_col, Nz, 1);

    % Time step (explicit stability)
    dt_max = 0.5*(mat.rho*mat.cp*dx^2*dz^2)/(mat.k*(dx^2 + dz^2));
    dt = cfl_safety*dt_max;

    Nt = floor(laser.tau/dt);
    if Nt < 5, Nt = 5; end
    if forceNtDivisibleBy10
        Nt = ceil(Nt/10)*10;
    end
    dt = laser.tau/Nt;
    tvec = linspace(0, laser.tau, Nt+1);

    % Snapshots
    snap_times = snap_fracs * laser.tau;
    snap_idx = zeros(size(snap_times));
    for k = 1:numel(snap_times)
        [~, snap_idx(k)] = min(abs(tvec - snap_times(k)));
    end
    snap_pct = round(100*snap_fracs);

    % Store evolving crater profiles at snapshots (for cumulative V/tau)
    depth_um_snap = nan(numel(snap_idx), Nx);
    t_ns_snap     = nan(1, numel(snap_idx));

    % Initialise
    T = mat.T0*ones(Nz,Nx);
    Tnew = T;

    surf_k = ones(1,Nx);
    Ecell  = zeros(1,Nx);
    cellVol = dx*dz*1e-6;

    Hremove_cell_J = mat.H_remove * cellVol;

    Melted    = false(size(T));
    Vaporized = false(size(T));

    Tmax     = T;
    MaxGradT = zeros(size(T));

    Hsnap_J = nan(1, numel(snap_idx));

    % Snapshot figure handles
    tileAxesT = []; tileAxesD = []; tileAxesP = [];
    tileAxesTz = []; tileAxesDz = []; tileAxesPz = [];

    if opts.doSnapshots
        % NOTE: temperature snapshots are now INSTANTANEOUS T (not Tmax)
        figT = figure('Name',[opts.caseName ' TEMP snapshots (ABS T, instantaneous)'], 'Color','w', 'Position',[60 80 1550 820]);
        tiledlayout(figT, 3, 4, 'Padding','compact', 'TileSpacing','compact');
        tileAxesT = gobjects(1, numel(snap_idx));
        for k = 1:numel(snap_idx), tileAxesT(k) = nexttile; end

        figD = figure('Name',[opts.caseName ' DEPTH snapshots (MIRRORED, 0–2.0um)'], 'Color','w', 'Position',[60 40 1550 820]);
        tiledlayout(figD, 3, 4, 'Padding','compact', 'TileSpacing','compact');
        tileAxesD = gobjects(1, numel(snap_idx));
        for k = 1:numel(snap_idx), tileAxesD(k) = nexttile; end

        figP = figure('Name',[opts.caseName ' PHASE snapshots'], 'Color','w', 'Position',[80 60 1550 820]);
        tiledlayout(figP, 3, 4, 'Padding','compact', 'TileSpacing','compact');
        tileAxesP = gobjects(1, numel(snap_idx));
        for k = 1:numel(snap_idx), tileAxesP(k) = nexttile; end

        figTz = figure('Name',[opts.caseName ' TEMP snapshots — ZOOMED (instantaneous)'], 'Color','w', 'Position',[1650 80 1550 820]);
        tiledlayout(figTz, 3, 4, 'Padding','compact', 'TileSpacing','compact');
        tileAxesTz = gobjects(1, numel(snap_idx));
        for k = 1:numel(snap_idx), tileAxesTz(k) = nexttile; end

        figDz = figure('Name',[opts.caseName ' DEPTH snapshots — ZOOMED'], 'Color','w', 'Position',[1650 40 1550 820]);
        tiledlayout(figDz, 3, 4, 'Padding','compact', 'TileSpacing','compact');
        tileAxesDz = gobjects(1, numel(snap_idx));
        for k = 1:numel(snap_idx), tileAxesDz(k) = nexttile; end

        figPz = figure('Name',[opts.caseName ' PHASE snapshots — ZOOMED'], 'Color','w', 'Position',[1670 60 1550 820]);
        tiledlayout(figPz, 3, 4, 'Padding','compact', 'TileSpacing','compact');
        tileAxesPz = gobjects(1, numel(snap_idx));
        for k = 1:numel(snap_idx), tileAxesPz(k) = nexttile; end

        % 0% renders (use INSTANTANEOUS T)
        t0 = 0;
        render_temp_snapshot_ABS(tileAxesT(1), 0, t0, T, x, z, surf_k, mat);
        render_depth_snapshot_mirrored(tileAxesD(1), 0, t0, x, surf_k, dz);
        [~, Phase3tmp0] = build_removed_and_phase(Melted, Vaporized, surf_k, Nz);
        render_phase_snapshot(tileAxesP(1), 0, t0, x, z, Phase3tmp0);

        render_temp_snapshot_ABS_ZOOM(tileAxesTz(1), 0, t0, T, x, z, surf_k, mat, zoom);
        render_depth_snapshot_mirrored_ZOOM(tileAxesDz(1), 0, t0, x, surf_k, dz, zoom);
        render_phase_snapshot_ZOOM(tileAxesPz(1), 0, t0, x, z, Phase3tmp0, zoom);

        Hsnap_J(1) = compute_domain_enthalpy_J(T, Melted, Vaporized, surf_k, mat, cellVol);

        depth_um_snap(1,:) = (surf_k-1)*dz*1e6;
        t_ns_snap(1) = 0.0;
    end

    % MAIN LOOP
    for n = 1:Nt
        t = tvec(n);

        % Temporal Gaussian (FWHM tau), peak at 0.5*tau
        pt = exp(-4*log(2)*((t - laser.t_peak)/laser.tau)^2);

        % Beer–Lambert from current surface
        pz = zeros(Nz,Nx);
        for j = 1:Nx
            ks = surf_k(j);
            if ks <= Nz
                zrel = z(ks:end) - z(ks);
                pz(ks:end,j) = exp(-laser.alpha_abs*zrel(:));
            end
        end

        A = laser.eta_abs*(laser.alpha_abs*I0) .* pt .* px .* pz; % W/m^3

        % Explicit diffusion
        Txx = (T(:,3:end)-2*T(:,2:end-1)+T(:,1:end-2))/dx^2;
        Tzz = (T(3:end,:)-2*T(2:end-1,:)+T(1:end-2,:))/dz^2;

        Tnew(2:end-1,2:end-1) = T(2:end-1,2:end-1) ...
            + mat.a_th*dt*(Txx(2:end-1,:) + Tzz(:,2:end-1)) ...
            + (dt/(mat.rho*mat.cp))*A(2:end-1,2:end-1);

        % Insulated (zero-flux) boundaries instead of forced ambient
        Tnew(:,1)   = Tnew(:,2);
        Tnew(:,end) = Tnew(:,end-1);
        Tnew(end,:) = Tnew(end-1,:);


        % Material mask now
        inMatNow = false(size(Tnew));
        for j = 1:Nx
            if surf_k(j) <= Nz
                inMatNow(surf_k(j):Nz, j) = true;
            end
        end

        % Phase history
        Melted    = Melted    | (inMatNow & (Tnew >= mat.Tm));
        Vaporized = Vaporized | (inMatNow & (Tnew >= mat.Tv));

        % Latent subtraction (visual)
        melt = (T < mat.Tm) & (Tnew >= mat.Tm);
        vap  = (T < mat.Tv) & (Tnew >= mat.Tv);
        Tnew(melt) = Tnew(melt) - (mat.Lm/mat.cp);
        Tnew(vap)  = Tnew(vap)  - (mat.Lv/mat.cp);

        % Air reset above surface; special handling at ks==1
        for j = 1:Nx
            ks = surf_k(j);
            if ks > 1
                Tnew(1:ks-1, j) = mat.T0;
            end
            if ks <= Nz && ks == 1
                term_source = (dt/(mat.rho*mat.cp)) * A(ks,j);
                if ks+1 <= Nz
                    term_cond   = (mat.a_th * dt / dz^2) * (T(ks+1,j) - T(ks,j));
                else
                    term_cond = 0;
                end
                Tnew(ks,j)  = T(ks,j) + term_source + term_cond;
            end
        end

        % ABLATION ENERGY BUCKET
        for j = 1:Nx
            ks = surf_k(j);
            if ks <= Nz
                if Tnew(ks,j) >= mat.Tv
                    Ecell(j) = Ecell(j) + A(ks,j)*dt*cellVol;

                    while (Ecell(j) >= mat.H_remove*cellVol) && (surf_k(j) <= Nz)
                        Ecell(j)  = Ecell(j) - mat.H_remove*cellVol;
                        surf_k(j) = surf_k(j) + 1;

                        if surf_k(j) <= Nz
                            Tnew(surf_k(j), j) = max(Tnew(surf_k(j), j), mat.Tv);
                        end
                    end
                end
            end
        end

        % History updates
        Tmax = max(Tmax, Tnew);
        [Gz, Gx] = gradient(Tnew, dz, dx);
        GradMag  = sqrt(Gx.^2 + Gz.^2);
        MaxGradT = max(MaxGradT, GradMag);

        % CommFit
        T = Tnew;

        % Snapshots
        if opts.doSnapshots
            hit = find((n+1) == snap_idx, 1);
            if ~isempty(hit)
                [~, Phase3tmp] = build_removed_and_phase(Melted, Vaporized, surf_k, Nz);

                t_snap = tvec(snap_idx(hit)); % actual snapshot time (s)

                % NORMAL (TEMP uses INSTANTANEOUS T)
                render_temp_snapshot_ABS(tileAxesT(hit), snap_pct(hit), t_snap, T, x, z, surf_k, mat);
                render_depth_snapshot_mirrored(tileAxesD(hit), snap_pct(hit), t_snap, x, surf_k, dz);
                render_phase_snapshot(tileAxesP(hit), snap_pct(hit), t_snap, x, z, Phase3tmp);

                % ZOOMED (TEMP uses INSTANTANEOUS T)
                render_temp_snapshot_ABS_ZOOM(tileAxesTz(hit), snap_pct(hit), t_snap, T, x, z, surf_k, mat, zoom);
                render_depth_snapshot_mirrored_ZOOM(tileAxesDz(hit), snap_pct(hit), t_snap, x, surf_k, dz, zoom);
                render_phase_snapshot_ZOOM(tileAxesPz(hit), snap_pct(hit), t_snap, x, z, Phase3tmp, zoom);

                % Enthalpy
                Hsnap_J(hit) = compute_domain_enthalpy_J(T, Melted, Vaporized, surf_k, mat, cellVol);

                % store crater profile for cumulative V/tau
                depth_um_snap(hit,:) = (surf_k-1)*dz*1e6;
                t_ns_snap(hit) = t_snap*1e9;
            end
        end
    end

    % Final depth
    crater_depth = (surf_k - 1)*dz;
    x_um = x*1e6;
    depth_um = crater_depth*1e6;

    % Removed mask
    Removed = false(Nz,Nx);
    for j = 1:Nx
        if surf_k(j) > 1
            Removed(1:surf_k(j)-1, j) = true;
        end
    end

    % Masks
    RecastMask = Melted & ~Removed;
    HAZMask    = (Tmax >= mat.T_HAZ) & ~Melted & ~Removed;
    CrackMask  = (MaxGradT >= mat.GradT_crit) & (Tmax >= mat.T_damage) & ~Removed;

    DamageMap = ones(Nz,Nx);
    DamageMap(Removed)     = 0;
    DamageMap(HAZMask)     = 2;
    DamageMap(RecastMask)  = 3;
    DamageMap(CrackMask)   = 4;

    cpu_s = toc(tStart);

    if opts.doAllPlots
        figure('Name',[opts.caseName ' Damage Map'], 'Color','w');
        imagesc(x_um, z*1e6, DamageMap);
        set(gca,'YDir','reverse'); axis tight;
        xlabel('X (\mum)'); ylabel('Z (\mum)');
        title('Damage Map: 0 Removed/Air, 1 Virgin, 2 HAZ, 3 Recast, 4 Crack-risk proxy');

        cmap = [
            1.00 1.00 1.00;  % 0
            0.00 0.10 0.40;  % 1
            1.00 0.90 0.10;  % 2
            1.00 0.55 0.00;  % 3
            0.90 0.00 0.00   % 4
        ];
        colormap(cmap); caxis([0 4]);
        cb = colorbar; cb.Ticks = 0:4;
        cb.TickLabels = {'Removed/Air','Virgin','HAZ','Recast','Crack-risk'};
    end

    res = struct();
    res.dx = dx; res.dz = dz; res.cfl_safety = cfl_safety;
    res.Nx = Nx; res.Nz = Nz;
    res.Nt = Nt; res.dt = dt;
    res.Nodes = Nx*Nz;
    res.cpu_s = cpu_s;

    res.x0 = x0;
    res.x_um = x_um;
    res.depth_um = depth_um;
    res.maxDepth_um = max(depth_um);

    res.Tmax   = Tmax;
    res.Tfinal = T;
    res.surf_k = surf_k;

    res.MaxGradT = MaxGradT;
    res.Removed = Removed;
    res.RecastMask = RecastMask;
    res.HAZMask = HAZMask;
    res.CrackMask = CrackMask;
    res.DamageMap = DamageMap;

    res.I0 = I0;
    res.px_col = px_col;

    res.Hremove_cell_J = Hremove_cell_J;
    res.Hsnap_J = Hsnap_J;
    res.snap_pct = snap_pct;

    res.depth_um_snap = depth_um_snap;
    res.t_ns_snap = t_ns_snap;
end

function H_J = compute_domain_enthalpy_J(T, Melted, Vaporized, surf_k, mat, cellVol)
    [Nz, Nx] = size(T);

    MaterialMask = false(Nz,Nx);
    for j = 1:Nx
        ks = surf_k(j);
        if ks <= Nz
            MaterialMask(ks:Nz, j) = true;
        end
    end

    dT = max(0, T - mat.T0);
    Hsensible = mat.rho * mat.cp .* dT;

    Hlatent = mat.rho * mat.Lm .* (Melted & MaterialMask) + ...
              mat.rho * mat.Lv .* (Vaporized & MaterialMask);

    Htot = (Hsensible + Hlatent) .* MaterialMask;
    H_J = sum(Htot(:)) * cellVol;
end

function px_col = gaussian_cell_average_x(x, x0, dx, d)
    a = 8/(d^2);
    xl = x - dx/2;
    xr = x + dx/2;
    C = sqrt(pi)/(2*sqrt(a));
    px_col = C * (erf(sqrt(a)*(xr - x0)) - erf(sqrt(a)*(xl - x0))) / dx;
end

function render_temp_snapshot_ABS(ax, pct, t_s, Tfield, x, z, surf_k, mat)
    axes(ax); cla(ax);

    Tmask = Tfield;
    Nx = length(surf_k);
    for j = 1:Nx
        ks = surf_k(j);
        if ks > 1
            Tmask(1:ks-1, j) = NaN;
        end
    end

    Tmin_plot = mat.Tmin_plot;
    Tmax_plot = mat.Tv;

    levels = linspace(Tmin_plot, Tmax_plot, 140);
    contourf(ax, x*1e6, z*1e6, Tmask, levels, 'LineColor', 'none');
    axis(ax, 'tight');
    set(ax, 'YDir', 'reverse');
    xlim(ax, [0 400]); ylim(ax, [0 80]);

    colormap(ax, jet(256));
    caxis(ax, [Tmin_plot Tmax_plot]);

    xlabel(ax, 'X (\mum)'); ylabel(ax, 'Z (\mum)');

    t_ns = t_s*1e9;
    title(ax, sprintf('Snapshot %d%%  (t = %.1f ns)', pct, t_ns));

    cb = colorbar(ax);
    cb.Label.String = 'T (K)';
    cb.Ticks = unique([Tmin_plot, round(mat.Tm), round(mat.Tv)]);

    hold(ax, 'on');
    Thaz = mat.T_HAZ_outline;
    try
        contour(ax, x*1e6, z*1e6, Tmask, [Thaz Thaz], ...
            'LineStyle', ':', 'LineWidth', 1.2, 'LineColor', [1 0.2 0.2]);
    catch
    end
    hold(ax, 'off');

    set(ax, 'Color', 'w');
end

function render_temp_snapshot_ABS_ZOOM(ax, pct, t_s, Tfield, x, z, surf_k, mat, zoom)
    axes(ax); cla(ax);

    Tmask = Tfield;
    Nx = length(surf_k);
    for j = 1:Nx
        ks = surf_k(j);
        if ks > 1
            Tmask(1:ks-1, j) = NaN;
        end
    end

    Tmin_plot = mat.Tmin_plot;
    Tmax_plot = mat.Tv;

    levels = linspace(Tmin_plot, Tmax_plot, 140);
    contourf(ax, x*1e6, z*1e6, Tmask, levels, 'LineColor', 'none');
    set(ax, 'YDir', 'reverse');

    x0_um = mean(x)*1e6;
    xlim(ax, [x0_um - zoom.x_halfwidth_um, x0_um + zoom.x_halfwidth_um]);
    ylim(ax, [0 zoom.z_max_um]);

    colormap(ax, jet(256));
    caxis(ax, [Tmin_plot Tmax_plot]);

    xlabel(ax, 'X (\mum)'); ylabel(ax, 'Z (\mum)');

    t_ns = t_s*1e9;
    title(ax, sprintf('Snapshot %d%% (Zoom)  (t = %.1f ns)', pct, t_ns));

    cb = colorbar(ax);
    cb.Label.String = 'T (K)';
    cb.Ticks = unique([Tmin_plot, round(mat.Tm), round(mat.Tv)]);

    hold(ax, 'on');
    Thaz = mat.T_HAZ_outline;
    try
        contour(ax, x*1e6, z*1e6, Tmask, [Thaz Thaz], ...
            'LineStyle', ':', 'LineWidth', 1.2, 'LineColor', [1 0.2 0.2]);
    catch
    end
    hold(ax, 'off');

    set(ax,'Color','w');
end

function render_depth_snapshot_mirrored(ax, pct, t_s, x, surf_k, dz)
    axes(ax); cla(ax);

    x_um = x*1e6;
    depth_um = (surf_k-1)*dz*1e6;

    plot(ax, x_um, depth_um, 'k-', 'LineWidth', 2);
    grid(ax,'on'); grid(ax,'minor');

    xlabel(ax,'X (\mum)'); ylabel(ax,'Depth (\mum)');

    t_ns = t_s*1e9;
    title(ax, sprintf('Snapshot %d%%  (t = %.1f ns)', pct, t_ns));

    xlim(ax, [0 400]);
    ylim(ax, [0 2.0]);
    set(ax, 'YDir','reverse');
    pbaspect(ax, [400 (2.0*6) 1]);
end

function render_depth_snapshot_mirrored_ZOOM(ax, pct, t_s, x, surf_k, dz, zoom)
    axes(ax); cla(ax);

    x_um = x*1e6;
    depth_um = (surf_k-1)*dz*1e6;

    plot(ax, x_um, depth_um, 'k-', 'LineWidth', 2);
    grid(ax,'on'); grid(ax,'minor');

    x0_um = mean(x_um);
    xlim(ax, [x0_um - zoom.x_halfwidth_um, x0_um + zoom.x_halfwidth_um]);

    ylim(ax, [0 2.0]);
    set(ax,'YDir','reverse');

    xlabel(ax,'X (\mum)'); ylabel(ax,'Depth (\mum)');

    t_ns = t_s*1e9;
    title(ax, sprintf('Snapshot %d%% (Zoom)  (t = %.1f ns)', pct, t_ns));

    pbaspect(ax, [2*zoom.x_halfwidth_um (2.0*6) 1]);
end

function render_phase_snapshot(ax, pct, t_s, x, z, Phase3)
    axes(ax); cla(ax);

    imagesc(ax, x*1e6, z*1e6, Phase3);
    set(ax,'YDir','reverse'); axis(ax,'tight');

    xlabel(ax,'X (\mum)'); ylabel(ax,'Z (\mum)');

    t_ns = t_s*1e9;
    title(ax, sprintf('Snapshot %d%%  (t = %.1f ns)', pct, t_ns));

    colormap(ax, [0 0 1; 1 0.5 0; 1 0 0]);
    caxis(ax, [1 3]);

    cb = colorbar(ax);
    cb.Ticks = [1 2 3];
    cb.TickLabels = {'Solid','Molten','Vaporised'};
    cb.Label.String = 'Phase';
end

function render_phase_snapshot_ZOOM(ax, pct, t_s, x, z, Phase3, zoom)
    axes(ax); cla(ax);

    imagesc(ax, x*1e6, z*1e6, Phase3);
    set(ax,'YDir','reverse');

    x_um = x*1e6;
    x0_um = mean(x_um);
    xlim(ax, [x0_um - zoom.x_halfwidth_um, x0_um + zoom.x_halfwidth_um]);
    ylim(ax, [0 zoom.z_max_um]);

    xlabel(ax,'X (\mum)'); ylabel(ax,'Z (\mum)');

    t_ns = t_s*1e9;
    title(ax, sprintf('Snapshot %d%% (Zoom)  (t = %.1f ns)', pct, t_ns));

    colormap(ax, [0 0 1; 1 0.5 0; 1 0 0]);
    caxis(ax,[1 3]);

    cb = colorbar(ax);
    cb.Ticks = [1 2 3];
    cb.TickLabels = {'Solid','Molten','Vaporised'};
    cb.Label.String = 'Phase';
end

function [Removed, Phase3] = build_removed_and_phase(Melted, Vaporized, surf_k, Nz)
    Nx = numel(surf_k);
    Removed = false(Nz, Nx);
    for j = 1:Nx
        if surf_k(j) > 1
            Removed(1:surf_k(j)-1, j) = true;
        end
    end
    Phase3 = ones(size(Removed));
    Phase3(Melted & ~Vaporized & ~Removed) = 2;
    Phase3(Vaporized | Removed) = 3;
end

function [V_pulse, MRR] = compute_volume_and_MRR(x_um, depth_um, f_Hz)
    % Volume removed per pulse (per 1 um width, 2D model)
    V_pulse = trapz(x_um*1e-6, depth_um*1e-6) * 1e-6;  
    % Standard MRR (pulse train equivalent)
    MRR = V_pulse * f_Hz;                              % m^3/s
end



function L2 = L2_on_common_support(x1, y1, x2, y2)
    xmin = max(min(x1), min(x2));
    xmax = min(max(x1), max(x2));
    if xmax <= xmin
        L2 = NaN; return;
    end
    xg = linspace(xmin, xmax, 400);
    y1g = interp1(x1, y1, xg, 'linear', 'extrap');
    y2g = interp1(x2, y2, xg, 'linear', 'extrap');
    L2 = sqrt(mean((y1g - y2g).^2));
end

function plot_relative_temporal_intensity(laser)
    t = linspace(0, laser.tau, 800);
    pt = exp(-4*log(2)*((t - 0.5*laser.tau)/laser.tau).^2);
    pt = pt / max(pt);

    figure('Name','Relative Temporal Intensity (Peak at 0.5\tau)','Color','w','Position',[120 120 900 480]);
    plot(t*1e9, pt, 'LineWidth', 2);
    grid on;
    xlabel('Time (ns)');
    ylabel('Relative Intensity (I/I_{peak})');
    title('Temporal Gaussian Pulse (peak at 0.5\tau)');

    safe_vline(0.5*laser.tau*1e9, '--', [0 0 0], 1.5, '0.5\tau');
    ylim([0 1.05]);
end

function plot_beam_2d_3d(laser)
    grid_size = 60e-6; N_pts = 400;
    xv = linspace(-grid_size/2, grid_size/2, N_pts);
    yv = linspace(-grid_size/2, grid_size/2, N_pts);
    [XX, YY] = meshgrid(xv, yv);
    RR = sqrt(XX.^2 + YY.^2);

    E_pulse = laser.P_av / laser.f;
    P_peak  = E_pulse / (laser.tau * sqrt(pi/(4*log(2))));
    I0      = (8 / (pi*laser.d_spot^2)) * P_peak;

    I = I0 * exp( -8 * (RR.^2) / (laser.d_spot^2) );
    I(RR > (laser.d_spot/2)) = 0;

    figure('Name', 'Laser Beam Intensity Profile (2D)', 'Color', 'w', 'Position', [100, 100, 650, 540]);
    surf(XX*1e6, YY*1e6, I, 'EdgeColor', 'none');
    view(2); axis square; axis tight;
    colormap(jet(256));
    cb = colorbar; cb.Label.String = 'Intensity (W/m^2)';
    xlabel('X (\mum)'); ylabel('Y (\mum)');
    title('2D Laser Beam Profile');

    figure('Name', 'Laser Beam Intensity Profile (3D)', 'Color', 'w', 'Position', [800, 100, 760, 600]);
    surf(XX*1e6, YY*1e6, I, 'EdgeColor', 'none');
    shading interp; colormap(jet(256));
    axis square; grid on; view(3);
    xlabel('X (\mum)'); ylabel('Y (\mum)'); zlabel('Intensity (W/m^2)');
    title('3D Beam Profile (Gaussian, hard-truncated)');
    cb = colorbar; cb.Label.String = 'Intensity (W/m^2)';
end

function plot_simulation_mesh(dom, dx, dz, zoom)
    x = 0:dx:dom.W;
    z = 0:dz:dom.H;
    [XX, ZZ] = meshgrid(x, z);

    Xum = XX * 1e6;
    Zum = ZZ * 1e6;

    Nmax = 250000;
    Npts = numel(Xum);
    if Npts > Nmax
        step = ceil(sqrt(Npts / Nmax));
        Xum_p = Xum(1:step:end, 1:step:end);
        Zum_p = Zum(1:step:end, 1:step:end);
        tag = sprintf(' (decimated 1/%d)', step);
    else
        Xum_p = Xum;
        Zum_p = Zum;
        tag = '';
    end

    figure('Name','Mesh used in simulation (full domain)','Color','w','Position',[120 120 980 520]);
    scatter(Xum_p(:), Zum_p(:), 14, 'r', 'filled');
    axis tight; axis ij; axis equal;
    xlabel('X (\mum)'); ylabel('Z (\mum)');
    title(sprintf('Simulation Mesh Nodes (dx=%.3f \\mum, dz=%.3f \\mum)%s', dx*1e6, dz*1e6, tag));
    grid on; grid minor;

    x0_um = dom.x0 * 1e6;
    xlimZ = [x0_um - zoom.x_halfwidth_um, x0_um + zoom.x_halfwidth_um];
    ylimZ = [0, zoom.z_max_um];

    figure('Name','Mesh used in simulation (zoom)','Color','w','Position',[120 680 980 520]);
    scatter(Xum_p(:), Zum_p(:), 14, 'r', 'filled');
    axis ij; axis equal;
    xlim(xlimZ); ylim(ylimZ);
    xlabel('X (\mum)'); ylabel('Z (\mum)');
    title(sprintf('Mesh Zoom (dx=%.3f \\mum, dz=%.3f \\mum)%s', dx*1e6, dz*1e6, tag));
    grid on; grid minor;
end

function plot_mesh_family_one_window(dom, dx_list, zoom)
    figure('Name','Mesh family (ALL mesh sizes in one window)','Color','w','Position',[80 80 1400 900]);
    tiledlayout(numel(dx_list),1,'Padding','compact','TileSpacing','compact');

    for i = 1:numel(dx_list)
        dx = dx_list(i);
        dz = dx/5;

        x = 0:dx:dom.W;
        z = 0:dz:dom.H;
        [XX, ZZ] = meshgrid(x, z);

        Xum = XX*1e6;
        Zum = ZZ*1e6;

        ax = nexttile;
        scatter(ax, Xum(:), Zum(:), 10, 'r', 'filled');
        axis(ax,'ij'); axis(ax,'equal');
        grid(ax,'on'); grid(ax,'minor');

        x0_um = dom.x0*1e6;
        xlim(ax, [x0_um-zoom.x_halfwidth_um, x0_um+zoom.x_halfwidth_um]);
        ylim(ax, [0, zoom.z_max_um]);

        title(ax, sprintf('Mesh Zoom (dx=%.1f \\mum, dz=%.3f \\mum)', dx*1e6, dz*1e6));
        xlabel(ax,'X (\mum)'); ylabel(ax,'Z (\mum)');
    end
end

function plot_topview_Tinst_with_crater(res, dom, mat)
    x  = res.x_um * 1e-6;
    x0 = dom.x0;
    Nx = res.Nx;

    Tsurf = nan(1, Nx);
    for j = 1:Nx
        ks = res.surf_k(j);
        ks = max(1, min(res.Nz, ks));
        Tsurf(j) = res.Tfinal(ks, j);
    end

    dz_um = res.dz * 1e6;
    surf_k_est = round(res.depth_um / dz_um) + 1;
    crater_cols = find(surf_k_est > 1);
    if isempty(crater_cols)
        r_crater = 0;
    else
        x_left  = x(min(crater_cols));
        x_right = x(max(crater_cols));
        r_crater = 0.5 * abs(x_right - x_left);
    end

    r_line = abs(x - x0);
    [r_sorted, idx] = sort(r_line, 'ascend');
    T_sorted = Tsurf(idx);

    [r_u, T_u] = make_unique_profile(r_sorted, T_sorted);

    W = dom.W;  Nxy = 420;
    xv = linspace(0, W, Nxy);
    yv = linspace(0, W, Nxy);
    [XX, YY] = meshgrid(xv, yv);
    RR = sqrt((XX - x0).^2 + (YY - x0).^2);

    Tgrid = interp1(r_u, T_u, RR, 'linear', mat.T0);
    if r_crater > 0
        Tgrid(RR <= r_crater) = NaN;
    end

    figure('Name','Top-view T (instantaneous) + crater','Color','w','Position',[80 80 1100 560]);
    imagesc(xv*1e6, yv*1e6, Tgrid);
    axis image;
    set(gca,'YDir','normal');
    xlabel('X (\mum)'); ylabel('Y (\mum)');
    title('Hypothetical Laser Top-view T');

    colormap(jet(256));
    set(gca,'Color',[1 1 1]);
    cb = colorbar; cb.Label.String = 'T (K)';
    caxis([mat.Tmin_plot mat.Tv]);
end

function [xu, yu] = make_unique_profile(x, y)
    [xu, ~, ic] = unique(x, 'stable');
    yu = accumarray(ic(:), y(:), [], @mean);
    [xu, ord] = sort(xu, 'ascend');
    yu = yu(ord);
end

function plot_Tmax_vs_distance_from_crater(res, dom)
    x = res.x_um * 1e-6;
    Nx = res.Nx;
    x0 = dom.x0;

    dz_um = res.dz * 1e6;
    surf_k_est = round(res.depth_um / dz_um) + 1;

    crater_cols = find(surf_k_est > 1);
    if isempty(crater_cols)
        r_crater = 0;
    else
        x_left  = x(min(crater_cols));
        x_right = x(max(crater_cols));
        r_crater = 0.5 * abs(x_right - x_left);
    end

    Tmax_surf = zeros(1, Nx);
    for j = 1:Nx
        ks = max(1, min(res.Nz, surf_k_est(j)));
        Tmax_surf(j) = res.Tmax(ks, j);
    end

    r = abs(x - x0);
    d = max(0, r - r_crater);
    d_um = d*1e6;

    nbins = 120;
    edges = linspace(0, max(d_um), nbins+1);
    centers = 0.5*(edges(1:end-1)+edges(2:end));
    Tbin = nan(1, nbins);

    for b = 1:nbins
        mask = (d_um >= edges(b)) & (d_um < edges(b+1));
        if any(mask)
            Tbin(b) = mean(Tmax_surf(mask));
        end
    end

    figure('Name','Paper-style: Tmax vs distance from crater','Color','w','Position',[1200 80 900 560]);
    plot(centers, Tbin, 'LineWidth', 2);
    grid on;
    xlabel('Distance from ablation crater (\mum)');
    ylabel('T_{max} (K)');
    title('T_{max} decay away from crater edge');
    xlim([0, max(centers)]);
end

function plot_centerline_Tinst_vs_depth_FINAL(res, dom, mat)
    Nz = res.Nz;
    Nx = res.Nx;

    z  = linspace(0, dom.H, Nz);
    j0 = round(Nx/2);

    % --- Surface index at beam centre (based on crater depth) ---
    dz_um      = res.dz * 1e6;
    surf_k_est = round(res.depth_um / dz_um) + 1;
    ks         = max(1, min(Nz, surf_k_est(j0)));

    % Start from TOP SURFACE node so local depth starts at 0
    kStart   = ks;
    z_rel_um = (z(kStart:end) - z(ks)) * 1e6;     % local depth below current surface (µm)

    % Temperatures along centreline (final time)
    T_raw  = res.Tfinal(kStart:end, j0);
    T_phys = min(T_raw, mat.Tv);                  % physical interpretation cap at Tv

    % Use Tmax to define HAZ (peak exposure)
    Tmax_line = res.Tmax(kStart:end, j0);

    % --- User presentation choices ---
    shift_um      = 2.0;                          % shift so surface appears at Z = 2 µm
    zmax_plot_um  = 11.0;                         % show up to 11 µm
    Thaz          = mat.T_HAZ_outline;            % HAZ threshold line (e.g., 773.15 K)

    % Keep points within x-range after shifting
    keep = (z_rel_um + shift_um <= zmax_plot_um);

    % Shifted depth axis for material profile
    Zp    = z_rel_um(keep) + shift_um;
    Tp    = T_phys(keep);
    Tmaxp = Tmax_line(keep);

    % --- Force ambient in removed/air region Z < 2, but KEEP a surface point at Z=2 ---
    Tsurf_cap = min(res.Tfinal(ks, j0), mat.Tv);  % capped surface/interface point at Z=2

    % Build final curve:
    %   [0, 2-) : air at T0
    %   Z=2     : surface/interface at capped T (≈Tv if boiling)
    %   >2      : material profile
    Z_full = [0, shift_um-1e-6, shift_um, Zp(:)'];
    T_full = [mat.T0, mat.T0, Tsurf_cap, Tp(:)'];

    % Ensure uniqueness / monotonicity
    [Z_full, uidx] = unique(Z_full, 'stable');
    T_full = T_full(uidx);

    % ---- Plot ----
    figure('Name','Temperature distribution across depth at beam centre (Snapshot 100%)', ...
           'Color','w','Position',[1200 650 980 540]);
    hold on; grid on;

    % Top limit for shading and axes
    yTop = max([max(T_full)*1.10, mat.Tv*1.10, mat.Tm*1.10, mat.T0*1.10]);

    % ---- Shade removed crater region (air): Z in [0,2] ----
    patch([0 shift_um shift_um 0], ...
          [mat.T0 mat.T0 yTop yTop], ...
          [0.85 0.85 0.85], 'EdgeColor','none', 'FaceAlpha',0.25);
    text(shift_um/2, mat.T0 + 40, 'Removed material / air', ...
         'HorizontalAlignment','center', 'FontSize', 9, 'Color',[0.4 0.4 0.4]);

    % ---- Shade MELTING and VAPORISATION temperature bands ----
    patch([0 zmax_plot_um zmax_plot_um 0], ...
          [mat.Tm mat.Tm mat.Tv mat.Tv], ...
          [1 0.90 0.60], 'EdgeColor','none', 'FaceAlpha',0.25);  % melting band

    patch([0 zmax_plot_um zmax_plot_um 0], ...
          [mat.Tv mat.Tv yTop yTop], ...
          [1 0.60 0.60], 'EdgeColor','none', 'FaceAlpha',0.18);  % vapor band

    % ---- Plot temperature curve ----
    plot(Z_full, T_full, '-o', 'LineWidth', 2.2, 'MarkerSize', 5);

    % ---- Reference lines ----
    yline(mat.Tm, '--k', 'T_m', 'LabelHorizontalAlignment','left');
    yline(mat.Tv, '--k', 'T_v', 'LabelHorizontalAlignment','left');
    yline(Thaz,  ':k', 'HAZ threshold', 'LabelHorizontalAlignment','left');

    % ---- HAZ depth (defined by Tmax >= Thaz, on shifted axis) ----
    haz_mask = (Tmaxp >= Thaz);
    if any(haz_mask)
        z_haz_end = max(Zp(haz_mask));   % already shifted
        xline(z_haz_end, 'r--', 'HAZ depth', ...
            'LabelVerticalAlignment','bottom', ...
            'LabelHorizontalAlignment','left', ...
            'LineWidth', 1.6);

        text(z_haz_end, Thaz, sprintf('  HAZ depth = %.2f \\mum', z_haz_end), ...
            'Color','r', 'FontWeight','bold', ...
            'VerticalAlignment','bottom', 'HorizontalAlignment','left');
    else
        text(shift_um, Thaz, 'HAZ not reached', ...
            'Color','r', 'FontWeight','bold', ...
            'VerticalAlignment','bottom', 'HorizontalAlignment','left');
    end

    % ---- Labels & legend ----
    xlabel('Depth below surface Z (\mum)');
    ylabel('Temperature (K)');
    title('Temperature distribution across depth at beam centre (100%)');

    xlim([0 zmax_plot_um]);
    ylim([mat.T0 yTop]);

    legend({'Removed material / air (Z<2)', ...
            'Melting region (T_m to T_v)', ...
            'Vaporisation region (T \geq T_v)', ...
            'T (physical, capped at T_v)'}, ...
            'Location','northeast');

    hold off;
end





function plot_cumulative_MRR_VoverTau(res, laser)
    if ~isfield(res,'depth_um_snap') || isempty(res.depth_um_snap)
        warning('No snapshot crater profiles stored. Enable opts.doSnapshots=true for cumulative V/tau plot.');
        return;
    end

    x_m = res.x_um * 1e-6;
    Ns = size(res.depth_um_snap, 1);

    V_m3 = nan(1, Ns);
    MRR = nan(1, Ns);

    for k = 1:Ns
        d_um = res.depth_um_snap(k,:);
        if any(isnan(d_um))
            continue;
        end
        d_m = d_um * 1e-6;
        V_m3(k) = trapz(x_m, d_m) * 1e-6;
        MRR(k)  = V_m3(k) * laser.f;   
    end

    pct = res.snap_pct(:)';

    figure('Name','Cumulative material removal rate','Color','w','Position',[120 120 1100 520]);
    yyaxis left
    plot(pct, V_m3, '-o', 'LineWidth', 2); grid on; hold on;
    xlabel('Snapshot (%)');
    ylabel('Cumulative removed volume V (m^3 per 1 um width)');
    title('Cumulative material removal rate');

    yyaxis right
    plot(pct, MRR, '-o', 'LineWidth', 2);
    ylabel('MRR_{equiv} = V \cdot f (m^3/s per 1 um width)');

    legend({'V (cumulative)','MRR = V·f (equiv)'}, 'Location','best');
end

%% ===================== SAFE LINE HELPERS =====================

function h = safe_hline(y, style, rgb, lw, labelStr)
    ax = gca;
    xL = xlim(ax);
    h = plot(ax, xL, [y y], style, 'LineWidth', lw);
    try, h.Color = rgb; catch, end
    if nargin >= 5 && ~isempty(labelStr)
        text(xL(1) + 0.02*(xL(2)-xL(1)), y, labelStr, ...
            'Color', rgb, 'FontWeight','bold', ...
            'VerticalAlignment','bottom', 'HorizontalAlignment','left');
    end
end

function h = safe_vline(x, style, rgb, lw, labelStr)
    ax = gca;
    yL = ylim(ax);
    h = plot(ax, [x x], yL, style, 'LineWidth', lw);
    try, h.Color = rgb; catch, end
    if nargin >= 5 && ~isempty(labelStr)
        text(x, yL(1) + 0.02*(yL(2)-yL(1)), labelStr, ...
            'Color', rgb, 'FontWeight','bold', ...
            'VerticalAlignment','bottom', 'HorizontalAlignment','left');
    end
end

function tune_eta_for_alpha()
% =====================================================================
% TUNE eta_abs TO KEEP CRATER DEPTH THE SAME WHEN CHANGING alpha_abs
%
% How to use:
%   1) Put this file in the SAME folder as your big model script (so
%      run_single_case_STABLE_FULL is on the MATLAB path).
%   2) Edit the "USER INPUTS" section below:
%        - alpha_baseline, eta_baseline (your current values)
%        - alpha_new (the new alpha you want, e.g. 4e6)
%        - eta_lo, eta_hi (search bounds)
%   3) Run:  tune_eta_for_alpha
%
% Output:
%   - prints eta_new that makes maxDepth_um match baseline
% =====================================================================

clc;

%% ========================= USER INPUTS =========================
alpha_baseline = 5.3e7;   % 1/m  (your current alpha)
eta_baseline   = 0.0170;  % your current eta_abs

alpha_new      = 4.0e6;   % 1/m  (target alpha)

% Binary search bounds for eta_abs (increase eta_hi if needed)
eta_lo = 1e-4;
eta_hi = 1.0;

tol_um  = 0.01;     % crater depth match tolerance in micrometers
maxIter = 30;       % binary search iterations
% ===============================================================

%% =================== DEFINE YOUR MODEL INPUTS ==================
% Paste your definitions here (same as your main script).
% Keep them identical to avoid mismatch.

% ---------------- MATERIAL: SS316 ----------------
mat.rho = 7900; mat.cp = 470; mat.k = 15;
mat.a_th = mat.k/(mat.rho*mat.cp);
mat.T0 = 293.15; mat.Tm = 1663.15; mat.Tv = 3073.15;
mat.Lm = 290e3; mat.Lv = 6090e3;
mat.h_remove_mass = mat.cp*(mat.Tm-mat.T0) + mat.Lm + mat.cp*(mat.Tv-mat.Tm) + mat.Lv;
mat.H_remove = mat.rho * mat.h_remove_mass;

mat.T_HAZ      = 0.60 * mat.Tm;
mat.T_damage   = 0.85 * mat.Tm;
mat.GradT_crit = 2.0e9;
mat.T_HAZ_outline = 773.15;
mat.Tmin_plot = 293.15;

% ---------------- LASER ----------------
laser.P_av   = 30;
laser.f      = 102e3;
laser.tau    = 140e-9;
laser.d_spot = 40e-6;
laser.t_peak = 0.5*laser.tau;

laser.eta_abs   = eta_baseline;
laser.alpha_abs = alpha_baseline;

% ---------------- DOMAIN ----------------
dom.W  = 400e-6;
dom.H  = 80e-6;
dom.x0 = dom.W/2;

% ---------------- DISCRETISATION ----------------
base.dx = 2e-6;
base.dz = 0.4e-6;
base.cfl_safety = 0.2;
base.forceNtDivisibleBy10 = true;

% ---------------- SNAPSHOTS (not used in calibration run) ----------------
snap_fracs = 0:0.1:1.0;

% ---------------- PROFILOMETER + ZOOM (required by your function signature) ----------------
surface_offset = 688.2;
zoom.x_halfwidth_um = 50;
zoom.z_max_um       = 15;

optsCal = struct('doSnapshots',false,'doAllPlots',false,'caseName','ETA_TUNER');
%% ===============================================================

%% =================== 1) RUN BASELINE (TARGET) ==================
res_base = run_single_case_STABLE_FULL(mat, laser, dom, base.dx, base.dz, base.cfl_safety, ...
    base.forceNtDivisibleBy10, snap_fracs, optsCal, surface_offset, zoom);

targetDepth_um = res_base.maxDepth_um;

fprintf('\nBASELINE:\n');
fprintf('  alpha = %.3e 1/m\n', alpha_baseline);
fprintf('  eta   = %.6f\n', eta_baseline);
fprintf('  depth = %.4f um\n\n', targetDepth_um);

%% =================== 2) SWITCH alpha, TUNE eta =================
laser.alpha_abs = alpha_new;

% Evaluate at eta_lo
laser.eta_abs = eta_lo;
res_lo = run_single_case_STABLE_FULL(mat, laser, dom, base.dx, base.dz, base.cfl_safety, ...
    base.forceNtDivisibleBy10, snap_fracs, optsCal, surface_offset, zoom);
d_lo = res_lo.maxDepth_um - targetDepth_um;

% Evaluate at eta_hi (expand if needed)
laser.eta_abs = eta_hi;
res_hi = run_single_case_STABLE_FULL(mat, laser, dom, base.dx, base.dz, base.cfl_safety, ...
    base.forceNtDivisibleBy10, snap_fracs, optsCal, surface_offset, zoom);
d_hi = res_hi.maxDepth_um - targetDepth_um;

% Expand eta_hi until the target is bracketed or we hit a safety limit
expandCount = 0;
while sign(d_lo) == sign(d_hi) && eta_hi < 50
    eta_hi = eta_hi * 2;
    laser.eta_abs = eta_hi;
    res_hi = run_single_case_STABLE_FULL(mat, laser, dom, base.dx, base.dz, base.cfl_safety, ...
        base.forceNtDivisibleBy10, snap_fracs, optsCal, surface_offset, zoom);
    d_hi = res_hi.maxDepth_um - targetDepth_um;
    expandCount = expandCount + 1;
    fprintf('[expand] eta_hi -> %.6f | depth=%.4f um\n', eta_hi, res_hi.maxDepth_um);
    if expandCount > 12, break; end
end

if sign(d_lo) == sign(d_hi)
    fprintf('\nERROR: Could not bracket target depth. Try different eta_lo/eta_hi.\n');
    fprintf('d_lo=%.4f um, d_hi=%.4f um\n', d_lo, d_hi);
    return;
end

% Binary search
eta_best = NaN;
res_best = [];
for it = 1:maxIter
    eta_mid = 0.5*(eta_lo + eta_hi);
    laser.eta_abs = eta_mid;

    res_mid = run_single_case_STABLE_FULL(mat, laser, dom, base.dx, base.dz, base.cfl_safety, ...
        base.forceNtDivisibleBy10, snap_fracs, optsCal, surface_offset, zoom);

    d_mid = res_mid.maxDepth_um - targetDepth_um;

    fprintf('[%02d] eta=%.8f | depth=%.4f um | err=%.4f um\n', it, eta_mid, res_mid.maxDepth_um, d_mid);

    eta_best = eta_mid;
    res_best = res_mid;

    if abs(d_mid) <= tol_um
        break;
    end

    if sign(d_mid) == sign(d_lo)
        eta_lo = eta_mid;
        d_lo = d_mid;
    else
        eta_hi = eta_mid;
        d_hi = d_mid;
    end
end

%% =================== 3) PRINT RESULT ===========================
fprintf('\nRESULT (plug this into your main code):\n');
fprintf('  NEW alpha_abs = %.3e 1/m\n', alpha_new);
fprintf('  TUNED eta_abs = %.8f\n', eta_best);
fprintf('  depth matched = %.4f um (target %.4f um)\n\n', res_best.maxDepth_um, targetDepth_um);

end
