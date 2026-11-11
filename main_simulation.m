%  Optimal Control of Multidrug-Resistant Tuberculosis Dynamics with Treatment Non-Adherence and Relapse
clear; clc; close all;
% Figures 
set(groot,'defaultAxesFontName','Times New Roman','defaultTextFontName','Times New Roman', ...
          'defaultLegendFontName','Times New Roman','defaultAxesFontSize',10, ...
          'defaultTextFontSize',10,'defaultLegendFontSize',10,'defaultLineLineWidth',1.4,'defaultFigureColor','w');

%% 1. SPECIFIED PARAMETERS AND CALIBRATION TARGETS ------------------------
N0 = 1e5;                   % per 100,000 population (Proposition 1)
fx.mu    = 1/74.15;         % life expectancy at birth 74.15 years, 2024 [BPS 2024]
fx.CFR_I = 118000/1080000;  % TB deaths (HIV-negative) / incident TB, 2024 [Kemenkes 2025, p. 2]
fx.fL    = 0.083;           % DS-TB lost to follow-up 8.3%, outcomes reported 2024 [Kemenkes 2025, p. 31]
fx.CFR_D = 16/(16+59);      % MDR/RR outcomes, cohort 2022: died 16%, success 59% [Kemenkes 2025, pp. 34, 37]
fx.D_D   = 20/12;           % MDR episode 20 months, longer regimens of 18-20 months [WHO 2025]
fx.rho   = 0.0941*12;       % loss of immunity 0.0941/month [Alhusna et al. 2025]
fx.eS    = 0.0106*12;       % relapse 0.0106/month [Alhusna et al. 2025]
fx.kE    = 1.0;             % eps_D = kE*eps_S (assumption, tested with kE = 2, 3)
fx.kappa = 1 - 1/4.33;      % DOT efficacy from PMO odds ratio 4.33 [Amala & Cahyati 2021]
fx.xi2   = 1/0.65 - 1;      % MDR/RR enrollment 65% in 2022 (cohort of CFR_D) -> 100% [Kemenkes 2025, p. 29]

% targets per 100,000: incidence 382 (1,080,000 cases) in 2024 [Kemenkes 2025, p. 2],
% prevalence 759 [Noviyani et al. 2021], MDR/RR among new cases 1.8% in 2024 [WHO global TB database],
% MDR/RR incidence 26,000 cases in 2024 [Kemenkes 2025, p. 2], i.e. 382*26000/1080000 per 100,000,
% latent infection 35.51% of the population in 2019 [Ding et al. 2022, Table 1]
pop = 281603799;            % population 2024 [Kemenkes 2025, p. 4]; national counts = per 100,000 x 2,816
dData = 41;                 % TB deaths (HIV-negative) per 100,000, 2024 [Kemenkes 2025, p. 2] (validation)
rrpData = 0.18;             % MDR/RR among previously treated cases, 2024 [WHO global TB database] (validation)
tg = [382, 759, 0.018, 382*26000/1080000, 0.3551*N0];

%% 2. CALIBRATION (Levenberg-Marquardt on log(theta)) ---------------------
th0 = [2.36 1.04 0.023 0.0083 1.96];
[th, maxres] = calibrate(fx, tg, th0, N0);
P  = derived(th, fx, N0);
X0 = equilibrium(P, N0);
ind = indicators(X0, P);
banner('CALIBRATION  (Section 2.3, Table 2, Table 3)');
fprintf('Estimated parameters (Eq. 6), max |log residual| = %.1e\n', maxres);
fprintf('  %-8s %10s   %s\n','beta',sprintf('%.4f',th(1)),'year^-1');
fprintf('  %-8s %10s   %s\n','phi',sprintf('%.4f',th(2)),'-');
fprintf('  %-8s %10s   %s\n','psi',sprintf('%.4f',th(3)),'-');
fprintf('  %-8s %10s   %s\n','tau',sprintf('%.5f',th(4)),'year^-1');
fprintf('  %-8s %10s   %s\n','D_I',sprintf('%.4f',th(5)),'years');
fprintf('Derived rates (Eq. 5), year^-1\n');
fprintf('  gamma = %.4f   sigma = %.4f   alpha_I = %.4f   delta = %.4f   alpha_D = %.4f\n', P.gam, P.sig, P.aI, P.dl, P.aD);
fprintf('\nFit to the calibration targets (Table 3, per 100,000)\n');
fprintf('  %-30s %10s %10s\n','Indicator','Target','Model');
fprintf('  %-30s %10.1f %10.1f\n','TB incidence', tg(1), ind.inc);
fprintf('  %-30s %10.1f %10.1f\n','TB prevalence', tg(2), ind.prev);
fprintf('  %-30s %9.2f%% %9.2f%%\n','MDR among new cases', 100*tg(3), 100*ind.rrn);
fprintf('  %-30s %10.2f %10.2f\n','MDR/RR-TB incidence', tg(4), ind.rrinc);
fprintf('  %-30s %10.0f %10.0f\n','Latent infection', tg(5), ind.ltbi);
fprintf('\nValidation, not used in the fit (Table 3, limitation iii)\n');
fprintf('  %-46s %8s %8s\n','Indicator','Data','Model');
fprintf('  %-46s %8.1f %8.1f\n','TB deaths, disease-induced (per 100,000)', dData, ind.deaths);
fprintf('  %-46s %8.1f %8.1f\n','TB deaths, all deaths during an episode', dData, (P.aI+P.mu)*X0(4)+(P.aD+P.mu)*X0(5));
fprintf('  %-46s %7.1f%% %7.1f%%\n','MDR among previously treated cases', 100*rrpData, 100*ind.rrp);
fprintf('  Previously treated inflow: loss to follow-up %.1f, relapse %.1f per 100,000 per year\n', ...
    P.sig*X0(4), P.eS*X0(6)+P.eD*X0(7));
KRS = P.eS+P.rho+P.mu;
fprintf('\nOther quantities quoted in the text\n');
fprintf('  2-year relapse probability               %6.2f%%\n', 100*P.eS/KRS*(1-exp(-2*KRS)));
fprintf('  Incidence from latent progression        %6.1f%%\n', 100*P.tau*(X0(2)+X0(3))/ind.inc);
fprintf('  MDR/RR incidence from acquired resistance %5.1f%%\n', 100*ind.acqshare);

%% 3. THRESHOLD ANALYSIS --------------------------------------------------
[RSS, RDD] = R0s(P);
banner('THRESHOLD AND EQUILIBRIUM  (Section 3.2)');
fprintf('  R_SS = %.4f    R_DD = %.4f    R_0 = max = %.4f\n', RSS, RDD, max(RSS,RDD));
fprintf('\nEndemic equilibrium (per 100,000), N = %.1f\n', sum(X0));
fprintf('  %10s %10s %10s %10s %10s %10s %10s\n','S','E_S','E_D','I','D','R_S','R_D');
fprintf('  %10.1f %10.1f %10.1f %10.2f %10.2f %10.2f %10.2f\n', X0);
ev = eig(numjac(@(x) rhs(x,0,0,P,N0), X0)); [~,ie] = sort(real(ev)); ev = ev(ie);
fprintf('\nEigenvalues of the Jacobian at the equilibrium (year^-1)\n');
for j = 1:numel(ev), fprintf('  %+.6f %+.6fi\n', real(ev(j)), imag(ev(j))); end
Pb = P; Pb.beta = P.beta/RSS;                      % beta* with R_SS(beta*) = 1
A  = numjac(@(x) rhs(x,0,0,Pb,N0), [N0 0 0 0 0 0 0]');
[VR,DR] = eig(A);  [~,k] = min(abs(diag(DR)));  w = real(VR(:,k)); if w(4)<0, w=-w; end
[VL,DL] = eig(A'); [~,k] = min(abs(diag(DL)));  v = real(VL(:,k)); v = v/(v'*w);
b_bif = v(2)*w(4) + Pb.phi*v(3)*w(5);
a_bif = 2*Pb.beta*w(1)/N0*b_bif*N0;               % proportional variables
fprintf('\nBifurcation coefficients at R_SS = 1 (Remark 1): a = %.4f, b = %.4f\n', a_bif, b_bif);

%% 4. OPTIMAL CONTROL: CASES ---------------------------------------------
% Weights of the original formulation: A2 = 100*A1 and W2 = 100*W1
oc.T = 5; oc.M = 1000; oc.A1 = 0.025; oc.A2 = 2.5; oc.W1 = 0.125; oc.W2 = 12.5;
oc.u1max = 1; oc.u2max = 1; oc.om = 0.5; oc.tol = 1e-6; oc.init = 0; oc.maxit = 1000;
names = {'No control','Case I','Case II','Case III'}; use = [0 0; 1 0; 0 1; 1 1];
R = cell(1,4); Mt = cell(1,4);
for c = 1:4
    o = oc; o.use1 = use(c,1); o.use2 = use(c,2);
    R{c} = solveOC(X0, P, o, N0);  Mt{c} = metrics(R{c}, P, o);
end
banner(sprintf('OUTCOMES OVER %d YEARS, per 100,000  (Table 4)', oc.T));
fprintf('  %-11s %9s %9s %9s %9s %8s %6s %6s %8s %5s\n','Strategy','PY DS-TB','PY MDR','MDR/RR','Acquired','Deaths','C1','C2','J','Iter');
for c = 1:4
    m = Mt{c};
    fprintf('  %-11s %9.1f %9.1f %9.1f %9.2f %8.1f %6.2f %6.2f %8.2f %5d\n', names{c}, ...
        m.PY_I, m.PY_D, m.CumRR, m.CumAcq, m.Deaths, m.C1, m.C2, m.J, m.it);
end
fprintf('\nReductions relative to no additional control (negative = increase)\n');
fprintf('  %-11s %9s %9s %9s %9s %11s\n','Strategy','DS-TB','MDR-TB','MDR/RR','Acquired','Incidence');
for c = 2:4
    fprintf('  %-11s %8.1f%% %8.1f%% %8.1f%% %8.1f%% %+10.2f%%\n', names{c}, ...
        100*(1-Mt{c}.PY_I/Mt{1}.PY_I), 100*(1-Mt{c}.PY_D/Mt{1}.PY_D), 100*(1-Mt{c}.CumRR/Mt{1}.CumRR), ...
        100*(1-Mt{c}.CumAcq/Mt{1}.CumAcq), 100*(Mt{c}.CumInc/Mt{1}.CumInc-1));
end
fprintf('\nCase III averts %.1f deaths and %.2f acquired-resistance events per 100,000\n', ...
    Mt{1}.Deaths-Mt{4}.Deaths, Mt{1}.CumAcq-Mt{4}.CumAcq);
fprintf('  nationally (x 2,816): %.0f deaths and %.0f acquired-resistance events\n', ...
    pop/1e5*(Mt{1}.Deaths-Mt{4}.Deaths), pop/1e5*(Mt{1}.CumAcq-Mt{4}.CumAcq));

banner('VALUES QUOTED IN SECTIONS 3.4-3.6');
for c = 2:4
    X = R{c}.X; t = R{c}.t;
    [Imin,ki] = min(X(:,4)); [Dmin,kd] = min(X(:,5)); [RSmax,ks] = max(X(:,6)); [RDmax,kr] = max(X(:,7));
    fprintf('%s\n', names{c});
    fprintf('  min I   = %7.2f at t = %.2f      min D   = %6.2f at t = %.2f\n', Imin, t(ki), Dmin, t(kd));
    fprintf('  max R_S = %7.2f at t = %.2f      max R_D = %6.2f at t = %.2f\n', RSmax, t(ks), RDmax, t(kr));
    if use(c,1), fprintf('  u1 leaves its upper bound at t = %.2f\n', t(find(R{c}.u1 < 0.999,1))); end
    if use(c,2), fprintf('  u2 leaves its upper bound at t = %.2f, u2(T/2) = %.3f\n', t(find(R{c}.u2 < 0.999,1)), R{c}.u2(round(end/2))); end
end
fprintf('max |u1(Case I) - u1(Case III)| = %.1e\n', max(abs(R{2}.u1-R{4}.u1)));
Lc3 = backwardRK4(R{4}.X, R{4}.u1, R{4}.u2, P, oc.T/oc.M, oc.A1, oc.A2, N0);
b1 = (1-P.psi)*Lc3(1,4)-Lc3(1,6); b2 = P.psi*Lc3(1,5);
fprintf('Bracket of u1* at t = 0 (Case III): DS-cure part %.4f, acquired-resistance part %.4f (%.1f%%)\n', b1, b2, 100*b2/(b1+b2));
cn = {'S','E_S','E_D','I','D','R_S','R_D'};
banner('COMPARTMENTS AT MONTH 60, per 100,000  (Table 5)');
fprintf('  %-5s %10s %20s %20s %20s\n','','No control','Case I','Case II','Case III');
for j = 1:7
    fprintf('  %-5s %10.2f', cn{j}, R{1}.X(end,j));
    for c = 2:4, fprintf('  %10.2f (%+6.2f%%)', R{c}.X(end,j), 100*(R{c}.X(end,j)/R{1}.X(end,j)-1)); end
    fprintf('\n');
end

% ---- Figure 3: trajectories ----
sty = {'-','--','-.','-'}; col = {[0 0 0],[0 0.45 0.74],[0.85 0.33 0.10],[0.47 0.67 0.19]};
lw  = [2.0 1.3 1.3 2.4]; ord = [1 4 3 2];        % thick lines underneath, thin lines on top
f3 = figure(3); clf(f3); set(f3,'WindowStyle','normal','Name','Figure 3: Trajectories','NumberTitle','off','Units','centimeters','Position',[1 1 16 16.5]); tiledlayout(3,3,'TileSpacing','compact');
ttl = {'(a) Susceptible, S','(b) Latent DS, E_S','(c) Latent MDR, E_D', ...
       '(d) Active DS-TB, I','(e) Active MDR-TB, D','(f) Recovered DS, R_S', ...
       '(g) Recovered MDR, R_D','(h) Acquired-resistance rate','(i) Cumulative MDR/RR cases'};
hL = gobjects(1,4);
for p = 1:9
    nexttile; hold on;
    for c = ord
        t = R{c}.t(:); X = R{c}.X; u1 = R{c}.u1(:);
        acq = P.psi*P.sig*(1-P.kappa*u1).*X(:,4);
        if p <= 7, y = X(:,p); elseif p == 8, y = acq; else, y = cumtrapz(t, P.tau*X(:,3)+P.eD*X(:,7)+acq); end
        hL(c) = plot(t, y, sty{c}, 'Color', col{c}, 'LineWidth', lw(c));
    end
    title(ttl{p}); xlabel('Time (years)'); xlim([0 oc.T]); xticks(0:1:oc.T);
    if p == 8, ylabel('per 100,000 per year'); else, ylabel('per 100,000'); end
    grid on; box on;
end
lg = legend(hL, names, 'Orientation','horizontal','FontSize',10); lg.Layout.Tile = 'south';
% exported like Figures 4-6 (cropped, at most 16 cm wide), so that it is not reduced in the paper;
% the size is set again just before export in case the window was resized or did not fit the screen
set(f3,'WindowState','normal','Units','centimeters'); f3.Position(3:4) = [16 16.5]; drawnow;
exportgraphics(f3,'Figure3_Trajectories.pdf','ContentType','vector');
fprintf('  Figure 3 exported from a %.1f x %.1f cm window (should be 16.0 x 16.5)\n', f3.Position(3:4));
if abs(f3.Position(3)-16) > 0.2, warning('Figure 3 window is not 16 cm wide: undock it, do not resize it, and run this section again.'); end

% ---- Figure 4: controls ----
f4 = figure(4); clf(f4); set(f4,'WindowStyle','normal','Name','Figure 4: Optimal controls','NumberTitle','off','Units','centimeters','Position',[2 2 16 7]); tiledlayout(1,2,'TileSpacing','compact');
nexttile; hold on;
h3 = plot(R{4}.t,R{4}.u1,'-','Color',col{4},'LineWidth',2.4);
h1 = plot(R{2}.t,R{2}.u1,'--','Color',col{2},'LineWidth',1.3);
xlim([0 oc.T]); xticks(0:1:oc.T); ylim([0 1.05]); title('(a) DOT supervision u_1^*(t)'); xlabel('Time (years)'); ylabel('Control intensity');
grid on; box on; legend([h1 h3],names([2 4]),'Location','southwest');
nexttile; hold on;
h3 = plot(R{4}.t,R{4}.u2,'-','Color',col{4},'LineWidth',2.4);
h2 = plot(R{3}.t,R{3}.u2,'-.','Color',col{3},'LineWidth',1.3);
xlim([0 oc.T]); xticks(0:1:oc.T); ylim([0 1.05]); title('(b) Second-line treatment u_2^*(t)'); xlabel('Time (years)'); ylabel('Control intensity');
grid on; box on; legend([h2 h3],names([3 4]),'Location','southwest');
drawnow; exportgraphics(f4,'Figure4_Optimal_controls.pdf','ContentType','vector');

%% 5. CONVERGENCE DIAGNOSTICS (Section 3.8, Table 6) ---------------------
banner('NUMERICAL CONVERGENCE, Case III  (Section 3.8, Table 6)');
fprintf('(a) Step-size refinement, u(0) = 0, omega = 0.5\n');
fprintf('  %8s %12s %14s\n','h (yr)','J','rel. change');
Ms = [250 500 1000 2000]; Jg = zeros(size(Ms));
for k = 1:numel(Ms)
    o = oc; o.M = Ms(k); o.use1 = 1; o.use2 = 1; mk = metrics(solveOC(X0,P,o,N0),P,o); Jg(k) = mk.J;
end
for k = 1:numel(Ms)
    fprintf('  %8.4f %12.5f %14.1e\n', oc.T/Ms(k), Jg(k), abs(Jg(k)-Jg(end))/Jg(end));
end
Jinit = [];
fprintf('(b) Initial guess and relaxation, h = 0.01\n');
fprintf('  %8s %7s %12s %6s\n','u(0)','omega','J','Iter');
for init = [0 0.5 1]
    for om = [0.3 0.5 0.8]
        o = oc; o.M = 500; o.use1 = 1; o.use2 = 1; o.init = init; o.om = om;
        r = solveOC(X0,P,o,N0); m = metrics(r,P,o); Jinit(end+1,:) = [init om m.J r.it]; %#ok<SAGROW>
        fprintf('  %5.1f*umax %5.1f %12.5f %6d\n', init, om, m.J, r.it);
    end
end

%% 6. PARAMETER UNCERTAINTY AND PRCC (Section 3.9, Table 7, Figure 5) --
banner('PARAMETER UNCERTAINTY  (Section 3.9, Figure 5, Table 7)  -- several minutes');
Nd = 200;
% columns: inc prev rrn rrinc ltbi | CFR_I fL CFR_D D_D rho eS kappa xi2
b25 = @(c) [0.75*c 1.25*c];          % +-25% where the source reports no interval
rngs = [b25(tg(1)); 589 961; b25(tg(3)); b25(tg(4)); b25(tg(5)); ...
        b25(fx.CFR_I); b25(fx.fL); b25(fx.CFR_D); b25(fx.D_D); b25(fx.rho); b25(fx.eS); ...
        1-1/1.15 1-1/16.32; 1/0.79-1 1/0.61-1];   % kappa: OR 95% CI; xi2: enrollment 61-79% in 2021-2024
if isfile('lhs_samples.csv')
    fid = fopen('lhs_samples.csv','r'); first = fgetl(fid); frewind(fid);
    hl = double(any(isletter(first)));                                 % 1 if the file has a header line
    C = textscan(fid, repmat('%f',1,13), 'Delimiter', {',',' ','\t'}, 'MultipleDelimsAsOne', true, 'HeaderLines', hl);
    fclose(fid); S = [C{:}];                                            % samples used in the paper
    if all(S(:) >= 0 & S(:) <= 1)                                       % unit Latin hypercube design (0-1):
        U = round(S, 5, 'significant');                                 % the design used in the paper (5 digits)
        S = rngs(:,1)' + U.*(rngs(:,2)-rngs(:,1))';                     % scaled to the parameter ranges
        fprintf('  LHS design in [0,1] scaled to the parameter ranges\n');
    end
    assert(isequal(size(S),[Nd 13]) && all(isfinite(S(:))), ...
        'lhs_samples.csv must contain 200 rows and 13 numeric columns (found %d x %d)', size(S,1), size(S,2));
    fprintf('  LHS file read: %d x %d; incidence %.1f-%.1f, latent %.0f-%.0f per 100,000\n', ...
        size(S,1), size(S,2), min(S(:,1)), max(S(:,1)), min(S(:,5)), max(S(:,5)));
else
    warning('lhs_samples.csv not found: a new sample is drawn, numbers will differ slightly.');
    rng(20261027); K = size(rngs,1); U = zeros(Nd,K);
    for j = 1:K, U(:,j) = (randperm(Nd)' - rand(Nd,1))/Nd; end
    S = rngs(:,1)' + U.*(rngs(:,2)-rngs(:,1))';
end
parNames = {'\beta','\phi','\psi','\tau','D_I','CFR_I','f_L','CFR_D','D_D','\rho','\epsilon_S','\kappa','\xi_2'};
keep = false(Nd,1); PAR = nan(Nd,13); OUT = nan(Nd,9); EXT = nan(Nd,6); REJ = nan(Nd,6);
o = oc; o.M = 200;
tLHS = tic;
for n = 1:Nd
    if mod(n,20) == 0                                  % progress report, no effect on results
        fprintf('  sample %3d of %d, %d accepted so far  (%.0f s)\n', n, Nd, sum(keep), toc(tLHS)); drawnow;
    end
    f = fx; f.CFR_I = S(n,6); f.fL = S(n,7); f.CFR_D = S(n,8); f.D_D = S(n,9);
    f.rho = S(n,10); f.eS = S(n,11); f.kappa = S(n,12); f.xi2 = S(n,13);
    [thn, err] = calibrate(f, S(n,1:5), th, N0);
    Pn = derived(thn, f, N0); Xn = equilibrium(Pn, N0); if isempty(Xn), continue; end
    [rs, rd] = R0s(Pn);
    mm = runCases(Xn, Pn, o, use, N0); mn = mm{1};
    rD = zeros(1,3); Jn = zeros(1,3); rA = zeros(1,3);
    for c = 2:4
        m = mm{c};
        rD(c-1) = 100*(1-m.PY_D/mn.PY_D); Jn(c-1) = m.J; rA(c-1) = 100*(1-m.CumAcq/mn.CumAcq);
    end
    if err > 1e-6                                        % no exact solution of Eq. (6):
        REJ(n,:) = [thn(3), rd/rs, err, rD];             % kept only as a diagnostic, at the best approximate fit
        continue;
    end
    keep(n) = true; PAR(n,:) = [thn, S(n,6:13)];
    OUT(n,:) = [rs rd rD rA(3) Jn];      % RSS RDD rD_I rD_II rD_III rAcq_III J_I J_II J_III
    % added value of Case III: extra MDR-TB reduction vs II (points), C2(III)/C2(II), C1(III)/C1(I),
    % DS-TB reduction (III), extra deaths averted vs II, extra MDR-TB reduction vs I (points)
    EXT(n,:) = [100*(mm{3}.PY_D-mm{4}.PY_D)/mn.PY_D, mm{4}.C2/mm{3}.C2, mm{4}.C1/mm{2}.C1, ...
                100*(1-mm{4}.PY_I/mn.PY_I), mm{3}.Deaths-mm{4}.Deaths, 100*(mm{2}.PY_D-mm{4}.PY_D)/mn.PY_D];
end
PAR = PAR(keep,:); OUT = OUT(keep,:); EXT = EXT(keep,:); nk = sum(keep);
if nk < 20                                             % diagnostic instead of an indexing error
    [~, e1] = calibrate(fx, S(1,1:5), th, N0);
    error(['Only %d of %d LHS samples were accepted. Residual of sample 1 with the baseline parameters: %.1e. ' ...
           'Check that lhs_samples.csv is the file sent with this code.'], nk, Nd, e1);
end
q = @(x) [median(x) prctile_(x,2.5) prctile_(x,97.5)];
fprintf('  Admissible samples: %d of %d (%d rejected)\n', nk, Nd, Nd-nk);
fprintf('  %-26s %8s %18s\n','Quantity','Median','95% interval');
fprintf('  %-26s %8.2f %8.2f - %7.2f\n','R_SS', q(OUT(:,1)));
fprintf('  %-26s %8.2f %8.2f - %7.2f\n','R_DD', q(OUT(:,2)));
for c = 1:3
    fprintf('  %-26s %7.1f%% %7.1f%% - %6.1f%%   (min %.1f%%)\n', ['MDR-TB reduction, ' names{c+1}], q(OUT(:,2+c)), min(OUT(:,2+c)));
end
fprintf('  P(R_DD > 1) = %.2f\n', mean(OUT(:,2)>1));
fprintf('  Case III has the lowest J in %.0f%% of samples (guaranteed, check only); J(II) < J(I) in %.0f%%\n', ...
    100*mean(OUT(:,9)<=OUT(:,7) & OUT(:,9)<=OUT(:,8)), 100*mean(OUT(:,8)<OUT(:,7)));
fprintf('\nAdded value of Case III (health benefits and control efforts, Section 3.9)\n');
fprintf('  extra MDR-TB reduction vs Case II  %5.1f points (%.1f - %.1f), min %.2f\n', q(EXT(:,1)), min(EXT(:,1)));
fprintf('  extra MDR-TB reduction vs Case I   %5.1f points (%.1f - %.1f)\n', q(EXT(:,6)));
fprintf('  C2(III)/C2(II) %.3f (%.3f - %.3f), max %.3f;  C1(III)/C1(I) %.4f - %.4f\n', q(EXT(:,2)), max(EXT(:,2)), min(EXT(:,3)), max(EXT(:,3)));
fprintf('  acquired-resistance reduction %.1f%% (%.1f - %.1f%%)\n', q(OUT(:,6)));
fprintf('  DS-TB reduction %.1f%% (%.1f - %.1f%%); extra deaths averted vs II %.1f (%.1f - %.1f)\n', q(EXT(:,4)), q(EXT(:,5)));

fprintf('\nExcluded samples (Section 3.9)\n');
rj = ~keep & ~isnan(REJ(:,1));
fprintf('  %d of %d rejected samples have an approximate equilibrium: psi <= %.2e, R_DD/R_SS >= %.5f\n', sum(rj), Nd-nk, max(REJ(rj,1)), min(REJ(rj,2)));
fprintf('  largest relative residual %.1e - %.3f (median %.3f)\n', min(REJ(rj,3)), max(REJ(rj,3)), median(REJ(rj,3)));
inNames = {'inc','prev','rrn','rrinc','ltbi','CFR_I','fL','CFR_D','D_D','rho','eS','kappa','xi2'};
fprintf('  %-6s %12s %12s %10s   (Mann-Whitney, normal approximation)\n','Input','retained','rejected','p');
for j = 1:13
    fprintf('  %-6s %12.4g %12.4g %10.1e\n', inNames{j}, median(S(keep,j)), median(S(~keep,j)), mwp_(S(keep,j), S(~keep,j)));
end
for c = 1:3
    fprintf('  rejected, %-9s MDR-TB reduction %6.2f%% (%.2f - %.2f%%), range %.3f - %.3f\n', names{c+1}, q(REJ(rj,3+c)), min(REJ(rj,3+c)), max(REJ(rj,3+c)));
end
for c = 1:3
    fprintf('  all samples, %-9s MDR-TB reduction %5.1f%% (%.1f - %.1f%%)\n', names{c+1}, q([OUT(:,2+c); REJ(rj,3+c)]));
end
pr = prcc(PAR, OUT(:,5));
dof = nk - 13 - 1;                                     % partial correlation with 12 controls
tcrit = fzero(@(t) betainc(dof/(dof+t^2), dof/2, 0.5) - 0.01, [1 10]);   % two-sided p = 0.01
rcrit = tcrit/sqrt(dof + tcrit^2);
fprintf('\nPRCC with the Case III MDR-TB reduction; significant if |PRCC| > %.3f (p < 0.01, %d d.o.f.)\n', rcrit, dof);
fprintf('  %-12s %8s   %-22s %s\n','Parameter','PRCC','95% range of samples','');
[~,ix] = sort(abs(pr),'descend');
for j = ix
    if abs(pr(j)) > rcrit, sg = 'significant'; else, sg = 'n.s.'; end
    fprintf('  %-12s %+8.3f   %9.4g - %-10.4g %s\n', strrep(strrep(parNames{j},'\',''),'_',''), pr(j), prctile_(PAR(:,j),2.5), prctile_(PAR(:,j),97.5), sg);
end
% Dependencies among the PRCC inputs (Section 3.9)
Rk = zeros(nk,13); for j = 1:13, Rk(:,j) = tiedrank_(PAR(:,j)); end
RS = corrcoef(Rk);                                     % Spearman rank correlations
[ii,jj] = find(triu(true(13),1)); [~,o2] = sort(abs(RS(sub2ind([13 13],ii,jj))),'descend');
fprintf('\nLargest Spearman correlations among the PRCC inputs\n');
for k = o2(1:6)'
    fprintf('  %-8s %-8s %+.2f\n', strrep(strrep(parNames{ii(k)},'\',''),'_',''), strrep(strrep(parNames{jj(k)},'\',''),'_',''), RS(ii(k),jj(k)));
end
VIF = zeros(1,13); R2 = zeros(1,13);
for j = 1:13
    Z = [ones(nk,1) Rk(:,setdiff(1:13,j))]; e = Rk(:,j) - Z*(Z\Rk(:,j));
    R2(j) = 1 - sum(e.^2)/sum((Rk(:,j)-mean(Rk(:,j))).^2); VIF(j) = 1/(1-R2(j));
end
fprintf('  Rank R^2 of each input on the other twelve inputs\n');
for j = 1:13, fprintf('    %-8s R2 = %.3f  VIF = %5.2f   simple rank corr. with outcome %+.3f\n', strrep(strrep(parNames{j},'\',''),'_',''), R2(j), VIF(j), corr_(Rk(:,j), tiedrank_(OUT(:,5)))); end
pr2 = prcc(S(keep,:), OUT(:,5));                      % PRCC with the independently sampled inputs
Rs2 = corrcoef(cell2mat(arrayfun(@(j) tiedrank_(S(keep,j)), 1:13, 'UniformOutput', false)));
fprintf('PRCC with the 13 independently sampled inputs (max |rank corr.| %.2f)\n', max(abs(Rs2(~eye(13)))));
[~,ix2] = sort(abs(pr2),'descend');
for j = ix2
    if abs(pr2(j)) > rcrit, sg = 'significant'; else, sg = 'n.s.'; end
    fprintf('  %-8s %+8.3f %s\n', inNames{j}, pr2(j), sg);
end

% ---- Figure 5: uncertainty ----
f5 = figure(5); clf(f5); set(f5,'WindowStyle','normal','Name','Figure 5: Uncertainty and PRCC','NumberTitle','off','Units','centimeters','Position',[2 2 16 13]); tiledlayout(2,2,'TileSpacing','compact');
nexttile(1); hold on;                                % (a) top left
ha = histogram(OUT(:,1),20,'FaceAlpha',0.6); hb = histogram(OUT(:,2),20,'FaceAlpha',0.6);
xline(1,'k--','LineWidth',1.5); ylim([0 1.3*max([ha.Values hb.Values])]);   % headroom for the legend
legend({'R_{SS}','R_{DD}','threshold'},'Orientation','horizontal','Location','north','FontSize',9,'Box','off');   % one row above the bars, covers no data
xlabel('Reproduction number'); ylabel('Samples');
title(sprintf('(a) P(R_{DD} > 1) = %.2f', mean(OUT(:,2)>1))); box on;
[~,ordp] = sort(abs(pr),'ascend'); vv = pr(ordp); vv = vv(:); isSig = abs(vv) > rcrit;
C = repmat([0.70 0.70 0.70],13,1);
C(vv>0 & isSig,:) = repmat([0.16 0.47 0.85],sum(vv>0 & isSig),1);
C(vv<0 & isSig,:) = repmat([0.84 0.24 0.24],sum(vv<0 & isSig),1);
nexttile(2,[2 1]); hold on;                          % (b) right column, full height
bh = barh(1:13, vv, 0.65, 'FaceColor','flat', 'EdgeColor','none'); bh.CData = C;
xline(0,'k-'); xline([-rcrit rcrit],'--','Color',[0.4 0.4 0.4],'LineWidth',1);
for j = 1:13
    lab = sprintf('%+.2f', vv(j)); if ~isSig(j), lab = [lab ' n.s.']; end %#ok<AGROW>
    xl = max(abs(vv(j)), rcrit) + 0.04;
    if vv(j) >= 0, text(xl, j, lab,'HorizontalAlignment','left','FontSize',9);
    else, text(-xl, j, lab,'HorizontalAlignment','right','FontSize',9); end
end
set(gca,'YTick',1:13,'YTickLabel',parNames(ordp),'TickLabelInterpreter','tex');
ylim([0.4 13.6]); xlim([-1.45 1.45]); xticks(-1:0.5:1);
xlabel('PRCC with MDR-TB reduction (Case III)'); title('(b) Global sensitivity'); grid on; box on;
nexttile(3); hold on;                                % (c) bottom left
h1 = histogram(OUT(:,3),15,'FaceAlpha',0.5); h2 = histogram(OUT(:,4),15,'FaceAlpha',0.5); h3 = histogram(OUT(:,5),15,'FaceAlpha',0.5);
ylim([0 1.3*max([h1.Values h2.Values h3.Values])]);                 % headroom for the legend
legend(names(2:4),'Orientation','horizontal','Location','north','FontSize',9,'Box','off'); xlabel('Reduction in MDR-TB person-years (%)'); ylabel('Samples');
title('(c) Benefit under uncertainty'); box on;
drawnow; exportgraphics(f5,'Figure5_Uncertainty_and_PRCC.pdf','ContentType','vector');

%% 7. SENSITIVITY TO WEIGHTS, EFFICACIES AND BOUNDS (Section 3.10, Figure 6)
banner('WEIGHTS, EFFICACIES, BOUNDS AND HORIZON  (Section 3.10, Figure 6)');
o = oc; o.M = 250; o.use1 = 0; o.use2 = 0; ncR = metrics(solveOC(X0,P,o,N0),P,o);
W2s = oc.W2*[0.1 1 3 10 30 100]; W1s = oc.W1*[0.1 1 3 10 30 100];
xis = [0.27 0.40 0.54 0.64]; kaps = [0.13 0.5 0.77 0.94];
HA = zeros(numel(xis),numel(W2s)); HB = zeros(numel(kaps),numel(W1s)); HBI = HB;
for i = 1:numel(xis)
    for j = 1:numel(W2s)
        Pi = P; Pi.xi2 = xis(i); o = oc; o.M = 250; o.use1 = 1; o.use2 = 1; o.W2 = W2s(j);
        m = metrics(solveOC(X0,Pi,o,N0),Pi,o); HA(i,j) = 100*(1-m.PY_D/ncR.PY_D);
    end
end
for i = 1:numel(kaps)
    for j = 1:numel(W1s)
        Pi = P; Pi.kappa = kaps(i); o = oc; o.M = 250; o.use1 = 1; o.use2 = 1; o.W1 = W1s(j);
        m = metrics(solveOC(X0,Pi,o,N0),Pi,o);
        HB(i,j) = 100*(1-m.CumAcq/ncR.CumAcq); HBI(i,j) = 100*(1-m.PY_I/ncR.PY_I);
    end
end
fprintf('  MDR-TB reduction over the (W2, xi2) grid:          %5.1f%% - %5.1f%%\n', min(HA(:)), max(HA(:)));
fprintf('  Acquired-resistance reduction over (W1, kappa):   %5.1f%% - %5.1f%%   (max DS-TB reduction %.1f%%)\n', ...
    min(HB(:)), max(HB(:)), max(HBI(:)));
JA2 = []; Jcor = []; Hz = []; Umax = []; KE = [];
fprintf('\nObjective J by case\n  %-24s %9s %9s %9s\n','Weights','Case I','Case II','Case III');
ADD = [];                          % added value of Case III for each weight setting
for A2 = oc.A1*[1 10 100]
    Jc = zeros(1,3); Mk = cell(1,3);
    for c = 2:4
        o = oc; o.M = 250; o.A2 = A2; o.use1 = use(c,1); o.use2 = use(c,2); mk = metrics(solveOC(X0,P,o,N0),P,o); Jc(c-1) = mk.J; Mk{c-1} = mk;
    end
    ADD(end+1,:) = addedValue(Mk, ncR); %#ok<SAGROW>
    JA2(end+1,:) = [A2/oc.A1 Jc]; %#ok<SAGROW>
    fprintf('  %-24s %9.2f %9.2f %9.2f\n', sprintf('A2/A1 = %d',round(A2/oc.A1)), Jc);
end
for w1 = oc.W1*[0.1 10]
    for w2 = oc.W2*[0.1 10]
        Jc = zeros(1,3); Mk = cell(1,3);
        for c = 2:4
            o = oc; o.M = 250; o.W1 = w1; o.W2 = w2; o.use1 = use(c,1); o.use2 = use(c,2);
            mk = metrics(solveOC(X0,P,o,N0),P,o); Jc(c-1) = mk.J; Mk{c-1} = mk;
        end
        ADD(end+1,:) = addedValue(Mk, ncR); %#ok<SAGROW>
        Jcor(end+1,:) = [w1 w2 Jc]; %#ok<SAGROW>
        fprintf('  %-24s %9.2f %9.2f %9.2f\n', sprintf('W1 = %g, W2 = %g',w1,w2), Jc);
    end
end
fprintf('\nAdded value of Case III in the seven weight settings above (rows in the same order)\n');
fprintf('  %10s %10s %10s %10s %12s %12s\n','MDR I','MDR II','MDR III','extra pts','acq. red.','C2 III/II');
for k = 1:size(ADD,1), fprintf('  %9.2f%% %9.2f%% %9.2f%% %10.2f %11.1f%% %12.4f\n', ADD(k,:)); end
fprintf('\nHorizon T (Case III reductions)\n  %6s %9s %9s %10s\n','T (yr)','MDR-TB','DS-TB','Acquired');
for T = [3 5 10]
    o = oc; o.T = T; o.M = 80*T; o.use1 = 0; o.use2 = 0; mn = metrics(solveOC(X0,P,o,N0),P,o);
    o.use1 = 1; o.use2 = 1; m = metrics(solveOC(X0,P,o,N0),P,o);
    Hz(end+1,:) = [T, 100*(1-m.PY_D/mn.PY_D), 100*(1-m.PY_I/mn.PY_I), 100*(1-m.CumAcq/mn.CumAcq)]; %#ok<SAGROW>
    fprintf('  %6d %8.1f%% %8.1f%% %9.1f%%\n', T, Hz(end,2:4));
end
fprintf('\nControl bounds (Case III)\n  %5s %32s %26s\n','bound','acquired-resistance (u1max)','MDR-TB (u2max)');
for um = [0.5 0.7 1.0]
    o = oc; o.M = 250; o.use1 = 1; o.use2 = 1; o.u1max = um; m1 = metrics(solveOC(X0,P,o,N0),P,o);
    o = oc; o.M = 250; o.use1 = 1; o.use2 = 1; o.u2max = um; m2 = metrics(solveOC(X0,P,o,N0),P,o);
    Umax(end+1,:) = [um, 100*(1-m1.CumAcq/ncR.CumAcq), 100*(1-m2.PY_D/ncR.PY_D)]; %#ok<SAGROW>
    fprintf('  %5.1f %31.1f%% %25.1f%%\n', um, Umax(end,2), Umax(end,3));
end
fprintf('\nRelapse after MDR cure (re-calibrated)\n  %-14s %8s %8s %10s %10s\n','eps_D','psi','R_DD','MDR red.','residual');
for kE = [2 3]
    f2 = fx; f2.kE = kE; [th2,e2] = calibrate(f2, tg, th, N0); P2 = derived(th2,f2,N0); X2 = equilibrium(P2,N0);
    o = oc; o.M = 250; o.use1 = 0; o.use2 = 0; n2 = metrics(solveOC(X2,P2,o,N0),P2,o);
    o.use1 = 1; o.use2 = 1; m2 = metrics(solveOC(X2,P2,o,N0),P2,o); [~,rd2] = R0s(P2);
    KE(end+1,:) = [kE, th2(3), rd2, 100*(1-m2.PY_D/n2.PY_D)]; %#ok<SAGROW>
    fprintf('  %-14s %8.4f %8.3f %9.1f%% %10.1e\n', sprintf('%d x eps_S',kE), th2(3), rd2, KE(end,4), e2);
end

%% 8. RE-ESTIMATION WITH THE VALIDATION INDICATORS (Section 3.11, limitation iii)
banner('RE-ESTIMATION WITH THE VALIDATION INDICATORS  (Section 3.11)');
VAL = [];
o = oc; o.M = 250;
fprintf('  %-26s %7s %7s %7s %8s %8s %7s %7s %7s %8s %8s %8s %8s\n','Re-estimation','psi','phi','R_DD','MDR/RR','acq.sh.', ...
        'MDR I','MDR II','MDR III','J(I)','J(II)','C2(III)','DS III');
for rrp = rrpData       % 18% MDR/RR among previously treated cases, 2024 [WHO global TB database]
    tgV = tg; tgV(4) = rrp;
    [thV, eV] = calibrate(fx, tgV, th, N0, true); PV = derived(thV, fx, N0); XV = equilibrium(PV, N0);
    iV = indicators(XV, PV); [~, rdV] = R0s(PV); mV = runCases(XV, PV, o, use, N0);
    rV = arrayfun(@(c) 100*(1-mV{c}.PY_D/mV{1}.PY_D), 2:4);
    VAL(end+1,:) = [rrp, thV(3), thV(2), rdV, iV.rrinc, iV.acqshare, rV, mV{2}.J, mV{3}.J, mV{4}.C2, 100*(1-mV{4}.PY_I/mV{1}.PY_I), eV]; %#ok<SAGROW>
    fprintf('  %-26s %7.3f %7.3f %7.3f %8.2f %8.3f %6.1f%% %6.1f%% %6.1f%% %8.2f %8.2f %8.2f %7.1f%%  (res. %.1e)\n', ...
            sprintf('MDR prev. treated = %g%%',100*rrp), VAL(end,2:end));
end
cfrM = fzero(@(c) deathsFor(c, fx, tg, th, N0) - dData, [fx.CFR_I 0.2]);  % disease-induced deaths = data
fM = fx; fM.CFR_I = cfrM; thM = calibrate(fM, tg, th, N0); PM = derived(thM, fM, N0); XM = equilibrium(PM, N0);
mM = runCases(XM, PM, o, use, N0); rM = arrayfun(@(c) 100*(1-mM{c}.PY_D/mM{1}.PY_D), 2:4);
fprintf('  Mortality %.1f: CFR_I = %.4f, MDR-TB reductions %.1f%%, %.1f%%, %.1f%%\n', dData, cfrM, rM);
mB = runCases(X0, P, o, use, N0); rB = arrayfun(@(c) 100*(1-mB{c}.PY_D/mB{1}.PY_D), 2:4);
fprintf('  Baseline (same grid):          MDR-TB reductions %.1f%%, %.1f%%, %.1f%%\n', rB);

% ---- Figure 6: robustness heat maps ----
f6 = figure(6); clf(f6); set(f6,'WindowStyle','normal','Name','Figure 6: Robustness','NumberTitle','off','Units','centimeters','Position',[2 2 16 12]); tiledlayout(2,1,'TileSpacing','compact');
heat(nexttile, HA, W2s, xis, 'Second-line weight W_2', 'Efficacy \xi_2', '(a) MDR-TB person-year reduction (%)');
heat(nexttile, HB, W1s, kaps, 'DOT weight W_1', 'Efficacy \kappa', '(b) Acquired-resistance reduction (%)');
drawnow; exportgraphics(f6,'Figure6_Robustness.pdf','ContentType','vector');

save('results_for_check.mat','fx','tg','th','P','X0','ind','RSS','RDD','ev','a_bif','b_bif','R','Mt','oc', ...
     'Ms','Jg','Jinit','HA','HB','HBI','JA2','Jcor','Hz','Umax','KE','OUT','PAR','nk','Nd','pr','rcrit', ...
     'EXT','REJ','keep','RS','R2','VIF','pr2','ADD','VAL','cfrM','rM','rB');
banner('DONE');
fprintf('Figures 3-6 saved as vector PDF; all results saved to results_for_check.mat\n');

%% =========================================================================
%                              LOCAL FUNCTIONS
%% =========================================================================
function banner(txt)
    fprintf('\n%s\n %s\n%s\n', repmat('=',1,78), txt, repmat('=',1,78));
    drawnow;                                           % show the output immediately
end

function P = derived(th, fx, N0)
    P.beta = th(1); P.phi = th(2); P.psi = th(3); P.tau = th(4); P.DI = th(5);
    P.mu = fx.mu; KI = 1/th(5);
    P.aI  = fx.CFR_I*KI - fx.mu;                       % Eq. (5)
    P.sig = fx.fL*KI/(1 - fx.fL*(1-th(3)));
    P.gam = KI*(1-fx.CFR_I) - th(3)*P.sig;
    KD = 1/fx.D_D; P.aD = fx.CFR_D*KD - fx.mu; P.dl = KD*(1-fx.CFR_D);
    P.eS = fx.eS; P.eD = fx.kE*fx.eS; P.rho = fx.rho;
    P.kappa = fx.kappa; P.xi2 = fx.xi2; P.pi = fx.mu*N0;
end

function [RSS, RDD] = R0s(P)
    KRS = P.eS+P.rho+P.mu; KRD = P.eD+P.rho+P.mu;
    KI = P.gam+P.psi*P.sig+P.mu+P.aI; KD = P.dl+P.mu+P.aD; f = P.tau/(P.tau+P.mu);
    RSS = P.beta*f*KRS/(KI*KRS - P.gam*P.eS);
    RDD = P.phi*P.beta*f*KRD/(KD*KRD - P.dl*P.eD);
end

function X = equilibrium(P, N0)                        % Proposition 3
    X = [];
    [RSS, RDD] = R0s(P); if RSS <= 1 || RDD >= RSS, return; end
    KRS = P.eS+P.rho+P.mu; KRD = P.eD+P.rho+P.mu; KD = P.dl+P.mu+P.aD;
    s  = 1/RSS;
    cd = P.psi*P.sig/((KD - P.eD*P.dl/KRD)*(1 - RDD*s));
    den = P.beta*s*(1+P.phi*cd) - P.rho*(P.gam/KRS + P.dl*cd/KRD);
    if den <= 0, return; end
    i = P.mu*(1-s)/den; d = cd*i;
    X = N0*[s; P.beta*s*i/(P.tau+P.mu); P.phi*P.beta*s*d/(P.tau+P.mu); i; d; P.gam*i/KRS; P.dl*d/KRD];
end

function ind = indicators(X, P)
    acq = P.psi*P.sig*X(4);
    ind.inc    = P.tau*(X(2)+X(3)) + P.eS*X(6) + P.eD*X(7);
    ind.prev   = X(4)+X(5);
    ind.rrn    = X(3)/(X(2)+X(3));
    ind.rrinc  = P.tau*X(3) + P.eD*X(7) + acq;
    ind.ltbi   = X(2)+X(3);
    ind.deaths = P.aI*X(4) + P.aD*X(5);
    ind.rrp    = (acq + P.eD*X(7))/(P.sig*X(4) + P.eS*X(6) + P.eD*X(7));
    ind.acqshare = acq/ind.rrinc;
end

function r = calres(lt, fx, tg, N0, useRrp)
    % useRrp = true replaces the MDR/RR-TB incidence target by MDR among previously treated cases
    th = exp(lt); P = derived(th, fx, N0);
    if P.aI <= 0 || P.gam <= 0 || P.aD <= 0, r = 10*ones(5,1); return; end
    X = equilibrium(P, N0); if isempty(X), r = 10*ones(5,1); return; end
    ind = indicators(X, P);
    if useRrp, v4 = ind.rrp; else, v4 = ind.rrinc; end
    r = log([ind.inc; ind.prev; ind.rrn; v4; ind.ltbi]./tg(:));
end

function [x, res] = lm(fun, x0)                        % Levenberg-Marquardt
    % Near-singular steps are expected for infeasible trial points and are
    % rejected by the damping rule below, so their warnings are suppressed.
    w1 = warning('off','MATLAB:nearlySingularMatrix'); w2 = warning('off','MATLAB:singularMatrix');
    c1 = onCleanup(@() warning(w1)); c2 = onCleanup(@() warning(w2));
    x = x0(:); r = fun(x); mu = 1e-3; n = numel(x);
    for it = 1:300
        J = zeros(numel(r), n);
        for j = 1:n
            h = 1e-7*max(1,abs(x(j))); e = zeros(n,1); e(j) = h;
            J(:,j) = (fun(x+e) - fun(x-e))/(2*h);
        end
        A = J'*J; g = J'*r;
        while true
            step = -(A + mu*diag(diag(A)+1e-300))\g;
            rn = fun(x+step);
            if all(isfinite(rn)) && rn'*rn < r'*r
                x = x+step; r = rn; mu = max(mu/3,1e-15); break;
            end
            mu = 4*mu; if mu > 1e12, break; end
        end
        if max(abs(r)) < 1e-13 || mu > 1e12, break; end
    end
    res = max(abs(r));
end

function [th, maxres] = calibrate(fx, tg, th0, N0, useRrp)
    if nargin < 5, useRrp = false; end
    starts = [th0; 2.0 0.6 0.05 0.009 2.0; 1.5 0.4 0.1 0.012 1.5; 3.0 0.8 0.03 0.007 2.5; 2.5 0.5 0.02 0.01 3.0];
    maxres = inf; th = th0;
    for k = 1:size(starts,1)
        [x, e] = lm(@(lt) calres(lt,fx,tg,N0,useRrp), log(starts(k,:)));
        if e < maxres, maxres = e; th = exp(x(:)'); end
        if maxres < 1e-10, break; end
    end
end

function dX = rhs(X, u1, u2, P, N0)                   % Eq. (3)
    S=X(1); ES=X(2); ED=X(3); I=X(4); D=X(5); RS=X(6); RD=X(7);
    acq = P.psi*P.sig*(1-P.kappa*u1); cureD = P.dl*(1+P.xi2*u2);
    dX = [P.pi + P.rho*(RS+RD) - P.beta*(I+P.phi*D)/N0*S - P.mu*S;
          P.beta*I/N0*S - (P.tau+P.mu)*ES;
          P.phi*P.beta*D/N0*S - (P.tau+P.mu)*ED;
          P.tau*ES + P.eS*RS - (P.gam + P.kappa*P.sig*u1 + acq + P.mu + P.aI)*I;
          P.tau*ED + P.eD*RD + acq*I - (cureD + P.mu + P.aD)*D;
          (P.gam + P.kappa*P.sig*u1)*I - (P.eS+P.rho+P.mu)*RS;
          cureD*D - (P.eD+P.rho+P.mu)*RD];
end

function dL = adjrhs(L, X, u1, u2, P, A1, A2, N0)      % Eq. (8)
    S=X(1); I=X(4); D=X(5);
    lS=L(1); lES=L(2); lED=L(3); lI=L(4); lD=L(5); lRS=L(6); lRD=L(7);
    acq = P.psi*P.sig*(1-P.kappa*u1); cureD = P.dl*(1+P.xi2*u2);
    QI = P.gam + P.kappa*P.sig*u1 + acq + P.mu + P.aI; QD = cureD + P.mu + P.aD;
    dL = [lS*(P.beta*(I+P.phi*D)/N0 + P.mu) - lES*P.beta*I/N0 - lED*P.phi*P.beta*D/N0;
          lES*(P.tau+P.mu) - lI*P.tau;
          lED*(P.tau+P.mu) - lD*P.tau;
          -A1 + (lS-lES)*P.beta*S/N0 + lI*QI - lD*acq - lRS*(P.gam + P.kappa*P.sig*u1);
          -A2 + (lS-lED)*P.phi*P.beta*S/N0 + lD*QD - lRD*cureD;
          -lS*P.rho - lI*P.eS + lRS*(P.eS+P.rho+P.mu);
          -lS*P.rho - lD*P.eD + lRD*(P.eD+P.rho+P.mu)];
end

function X = forwardRK4(X0, u1, u2, P, h, N0)
    M = numel(u1)-1; X = zeros(M+1,7); X(1,:) = X0(:)';
    for k = 1:M
        m1 = 0.5*(u1(k)+u1(k+1)); m2 = 0.5*(u2(k)+u2(k+1)); x = X(k,:)';
        k1 = rhs(x,u1(k),u2(k),P,N0);      k2 = rhs(x+0.5*h*k1,m1,m2,P,N0);
        k3 = rhs(x+0.5*h*k2,m1,m2,P,N0);   k4 = rhs(x+h*k3,u1(k+1),u2(k+1),P,N0);
        X(k+1,:) = (x + h/6*(k1+2*k2+2*k3+k4))';
    end
end

function L = backwardRK4(X, u1, u2, P, h, A1, A2, N0)
    M = numel(u1)-1; L = zeros(M+1,7);                 % transversality lambda(T) = 0
    for k = M+1:-1:2
        xm = 0.5*(X(k,:)+X(k-1,:))'; m1 = 0.5*(u1(k)+u1(k-1)); m2 = 0.5*(u2(k)+u2(k-1)); l = L(k,:)';
        k1 = adjrhs(l,X(k,:)',u1(k),u2(k),P,A1,A2,N0);  k2 = adjrhs(l-0.5*h*k1,xm,m1,m2,P,A1,A2,N0);
        k3 = adjrhs(l-0.5*h*k2,xm,m1,m2,P,A1,A2,N0);    k4 = adjrhs(l-h*k3,X(k-1,:)',u1(k-1),u2(k-1),P,A1,A2,N0);
        L(k-1,:) = (l - h/6*(k1+2*k2+2*k3+k4))';
    end
end

function r = solveOC(X0, P, o, N0)                     % Forward-Backward Sweep
    h = o.T/o.M; t = linspace(0,o.T,o.M+1);
    u1 = o.init*o.u1max*o.use1*ones(1,o.M+1); u2 = o.init*o.u2max*o.use2*ones(1,o.M+1);
    rel = @(a,b) sum(abs(a(:)-b(:)))/max(sum(abs(a(:))),1e-300);
    Xo = []; Lo = []; it = 0;
    if o.use1 || o.use2
        for it = 1:o.maxit
            X = forwardRK4(X0,u1,u2,P,h,N0);
            L = backwardRK4(X,u1,u2,P,h,o.A1,o.A2,N0);
            v1 = P.kappa*P.sig*X(:,4)'.*((1-P.psi)*L(:,4)' + P.psi*L(:,5)' - L(:,6)')/o.W1;   % Eq. (10)
            v2 = P.xi2*P.dl*X(:,5)'.*(L(:,5)' - L(:,7)')/o.W2;                               % Eq. (11)
            v1 = min(o.u1max, max(0, v1))*o.use1;  v2 = min(o.u2max, max(0, v2))*o.use2;
            n1 = o.om*v1 + (1-o.om)*u1;  n2 = o.om*v2 + (1-o.om)*u2;                        % relaxation
            if isempty(Xo), e = 1; else
                e = max([rel(n1,u1)*o.use1, rel(n2,u2)*o.use2, rel(X,Xo), rel(L,Lo)]);
            end
            u1 = n1; u2 = n2; Xo = X; Lo = L;
            if e < o.tol, break; end
        end
    end
    r.t = t; r.u1 = u1; r.u2 = u2; r.it = it; r.X = forwardRK4(X0,u1,u2,P,h,N0);
end

function m = metrics(r, P, o)
    t = r.t(:); X = r.X; u1 = r.u1(:); u2 = r.u2(:); I = X(:,4); D = X(:,5);
    acq = P.psi*P.sig*(1-P.kappa*u1).*I;
    m.PY_I = trapz(t,I); m.PY_D = trapz(t,D);
    m.CumRR  = trapz(t, P.tau*X(:,3) + P.eD*X(:,7) + acq);
    m.CumAcq = trapz(t, acq);
    m.CumInc = trapz(t, P.tau*(X(:,2)+X(:,3)) + P.eS*X(:,6) + P.eD*X(:,7));
    m.Deaths = trapz(t, P.aI*I + P.aD*D);
    m.C1 = trapz(t, 0.5*o.W1*u1.^2); m.C2 = trapz(t, 0.5*o.W2*u2.^2);
    m.J  = o.A1*m.PY_I + o.A2*m.PY_D + m.C1 + m.C2; m.it = r.it;
end

function J = numjac(f, x)
    n = numel(x); J = zeros(n);
    for j = 1:n
        hh = 1e-4*max(1,abs(x(j))); e = zeros(n,1); e(j) = hh;
        J(:,j) = (f(x+e) - f(x-e))/(2*hh);
    end
end

function heat(ax, H, xv, yv, xl, yl, ttl)
    imagesc(ax, H); colormap(ax, parula); colorbar(ax); set(ax,'YDir','normal', ...
        'XTick',1:numel(xv),'XTickLabel',compose('%g',xv),'YTick',1:numel(yv),'YTickLabel',compose('%g',yv));
    xlabel(ax, xl); ylabel(ax, yl); title(ax, ttl);
    cl = caxis(ax); mid = mean(cl);
    for i = 1:size(H,1)
        for j = 1:size(H,2)
            if H(i,j) > mid, tc = 'k'; else, tc = 'w'; end
            text(ax, j, i, sprintf('%.0f',H(i,j)),'HorizontalAlignment','center','Color',tc,'FontWeight','bold','FontSize',10);
        end
    end
end

function p = prctile_(x, q)                            % percentile, linear interpolation
    x = sort(x(:)); n = numel(x); pos = 1 + (n-1)*q/100;
    lo = floor(pos); hi = ceil(pos); p = x(lo) + (pos-lo)*(x(hi)-x(lo));
end

function r = tiedrank_(x)
    [~,idx] = sort(x); r = zeros(size(x)); r(idx) = 1:numel(x);
end

function pr = prcc(Xp, y)
    [n,k] = size(Xp); Rk = zeros(n,k);
    for j = 1:k, Rk(:,j) = tiedrank_(Xp(:,j)); end
    ry = tiedrank_(y); pr = zeros(1,k);
    for j = 1:k
        Z = [ones(n,1) Rk(:,setdiff(1:k,j))];
        ex = Rk(:,j) - Z*(Z\Rk(:,j)); ey = ry - Z*(Z\ry);
        c = corrcoef(ex, ey); pr(j) = c(1,2);
    end
end

function m = runCases(X0, P, o, use, N0)               % no control and Cases I-III
    m = cell(1,4);
    for c = 1:4
        o.use1 = use(c,1); o.use2 = use(c,2); m{c} = metrics(solveOC(X0,P,o,N0),P,o);
    end
end

function a = addedValue(Mk, nc)                        % Mk = {Case I, Case II, Case III}
    r = cellfun(@(m) 100*(1-m.PY_D/nc.PY_D), Mk);
    a = [r, r(3)-max(r(1:2)), 100*(1-Mk{3}.CumAcq/nc.CumAcq), Mk{3}.C2/Mk{2}.C2];
end

function d = deathsFor(cfr, fx, tg, th0, N0)           % disease-induced deaths after re-estimation
    f = fx; f.CFR_I = cfr; th = calibrate(f, tg, th0, N0); P = derived(th, f, N0);
    ind = indicators(equilibrium(P, N0), P); d = ind.deaths;
end

function p = mwp_(a, b)                                % two-sided Mann-Whitney test, normal approximation
    n1 = numel(a); n2 = numel(b); r = tiedrank_([a(:); b(:)]);
    U = sum(r(1:n1)) - n1*(n1+1)/2;
    z = (U - n1*n2/2)/sqrt(n1*n2*(n1+n2+1)/12); p = erfc(abs(z)/sqrt(2));
end

function c = corr_(a, b)
    C = corrcoef(a, b); c = C(1,2);
end
