%% fourbar_symbolic.m  -  Symbolic kinematics & dynamic force analysis (Symbolic Math Toolbox)
% Same mechanism and frame as fourbar_cad_v2.m (Fusion "Crank-Rocker_final"):
%   O2 = (0,0) disc axle, O4 = (-g,0) rocker pivot, all angles counter-clockwise from +x.
%   theta2 = crank O2->A, theta3 = coupler A->B, theta4 = rocker O4->B.
% Derives position, velocity, acceleration, transmission angle, limit positions, time ratio,
% working space, AND planar Newton-Euler dynamic forces/torque symbolically, then evaluates them
% for the Fusion CAD dimensions/inertias (PLA 1240 kg/m^3) and compares with Simscape Multibody.
clear;
here = fileparts(mfilename('fullpath')); if ~isempty(here), cd(here); end   % run from any folder

%% 1) Symbols and vector loop  r2 + r3 = r1 + r4   (r1 = O2->O4 = (-g, 0))
syms a b c g positive
syms th2 th3 th4 w2 real
fx = a*cos(th2) + b*cos(th3) - (-g + c*cos(th4));
fy = a*sin(th2) + b*sin(th3) -  c*sin(th4);
f  = [fx; fy];
disp('Loop-closure equations (= 0):'); disp(f)

%% 2) Eliminate theta3:  |r3|^2 = b^2  ->  P cos(theta4) + Q sin(theta4) + R = 0
E = (-g + c*cos(th4) - a*cos(th2))^2 + (c*sin(th4) - a*sin(th2))^2 - b^2;
P = -2*c*(g + a*cos(th2));
Q = -2*a*c*sin(th2);
R =  a^2 - b^2 + c^2 + g^2 + 2*a*g*cos(th2);
syms t real
chk = @(ex) double(subs(ex, [a b c g th2 th4 t], [0.05 0.15 0.14 0.1 rand*2*pi rand*2*pi rand-0.5]));
assert(abs(chk(E - (P*cos(th4) + Q*sin(th4) + R))) < 1e-12, 'Derivation check of P, Q, R failed');
disp('Freudenstein form  P*cos(th4) + Q*sin(th4) + R = 0  with'); disp([P; Q; R])

% Tangent half-angle t = tan(theta4/2) -> (R-P) t^2 + 2Q t + (P+R) = 0
quadEq = (R - P)*t^2 + 2*Q*t + (P + R);
assert(abs(chk(subs(P*cos(th4) + Q*sin(th4) + R, th4, 2*atan(t))*(1 + t^2) - quadEq)) < 1e-12, ...
       'Half-angle substitution check failed');
Dq = P^2 + Q^2 - R^2;                                % discriminant
th4_open  = 2*atan((-Q + sqrt(Dq))/(R - P));         % circuit 1 (B above line A-O4, as in Fusion)
th4_cross = 2*atan((-Q - sqrt(Dq))/(R - P));         % circuit 2 (mirror assembly)
th3_expr  = atan2(c*sin(th4) - a*sin(th2), -g + c*cos(th4) - a*cos(th2));

%% 3) Velocity: differentiate the loop in time,  J * [w3; w4] = -df/dth2 * w2
Jq = jacobian(f, [th3 th4]);
Jt = jacobian(f, th2);
wv = simplify(-(Jq \ Jt) * w2);
w3 = wv(1);  w4 = wv(2);
disp('Coupler angular velocity w3 ='); pretty(w3)
disp('Rocker angular velocity  w4 ='); pretty(w4)

%% 4) Acceleration (constant crank speed, alpha2 = 0):  d/dt = w2 d/dth2 + w3 d/dth3 + w4 d/dth4
alv = jacobian(wv, [th2 th3 th4]) * [w2; w3; w4];
al3 = alv(1);  al4 = alv(2);

%% 5) Transmission angle mu (between coupler and rocker at B), Grashof, limit positions
d2   = a^2 + g^2 + 2*a*g*cos(th2);                   % |A O4|^2
mu   = acos((b^2 + c^2 - d2)/(2*b*c));
muMin = acos((b^2 + c^2 - (g - a)^2)/(2*b*c));       % theta2 = 180 deg
muMax = acos((b^2 + c^2 - (g + a)^2)/(2*b*c));       % theta2 = 0 deg
% Rocker limit positions: crank and coupler collinear, |O2 B| = b + a (extended) or b - a (folded)
th4_ext  = acos((g^2 + c^2 - (b + a)^2)/(2*g*c));
th4_fold = acos((g^2 + c^2 - (b - a)^2)/(2*g*c));

%% 6) Symbolic Dynamic Force & Torque Analysis (Newton-Euler planar equations, no friction)
syms w3s w4s al3s al4s positive_dummy
clear positive_dummy
syms w3s w4s al3s al4s real
syms m2d r2d mAh m3 r3c IG3 m4 r4c IG4 real

e2 = [cos(th2); sin(th2)];   n2 = [-sin(th2); cos(th2)];
e3 = [cos(th3); sin(th3)];   n3 = [-sin(th3); cos(th3)];
e4 = [cos(th4); sin(th4)];   n4 = [-sin(th4); cos(th4)];

% CoM accelerations (alpha2 = 0, constant crank speed)
aA   = -a*w2^2*e2;
aG2d = -r2d*w2^2*e2;
aG3  = aA + r3c*(al3s*n3 - w3s^2*e3);
aG4  = r4c*(al4s*n4 - w4s^2*e4);

% Position vectors from each link CoM to its joints
r2_O2 = -r2d*e2;       r2_A = (a - r2d)*e2;
r3_A  = -r3c*e3;       r3_B = (b - r3c)*e3;
r4_O4 = -r4c*e4;       r4_B = (c - r4c)*e4;

% Unknowns: xF = [F12x; F12y; F32x; F32y; F43x; F43y; F14x; F14y; T2]
%   F12 = force of ground on crank at O2; F32 = force of coupler on crank at A
%   F43 = force of rocker on coupler at B; F14 = force of ground on rocker at O4; T2 = input crank torque
cross2 = @(r, F) r(1)*F(2) - r(2)*F(1);
syms F12x F12y F32x F32y F43x F43y F14x F14y T2 real
F12 = [F12x; F12y]; F32 = [F32x; F32y]; F43 = [F43x; F43y]; F14 = [F14x; F14y];
xF  = [F12x; F12y; F32x; F32y; F43x; F43y; F14x; F14y; T2];

eqF = [
    F12 + F32 - (m2d*aG2d + mAh*aA);
    cross2(r2_O2, F12) + cross2(r2_A, F32) + T2 - mAh*cross2(r2_A, aA);
   -F32 + F43 - m3*aG3;
    cross2(r3_A, -F32) + cross2(r3_B, F43) - IG3*al3s;
    F14 - F43 - m4*aG4;
    cross2(r4_O4, F14) + cross2(r4_B, -F43) - IG4*al4s
];
[AF_sym, bF_sym] = equationsToMatrix(eqF == 0, xF);

%% 7) Numbers from the Fusion model (SI units, PLA rho = 1240 kg/m^3)
num   = struct('a', 0.050, 'b', 0.150, 'c', 0.140, 'g', 0.100, 'w2', 2*pi);
symv  = [a b c g w2];
numv  = [num.a num.b num.c num.g num.w2];

% Exact inertial properties from Fusion STEP CAD (PLA 1240 kg/m^3)
m2_disc = 0.1761919;  r2_disc = -0.00017952582;
m_pinA  = 1240 * pi * 0.004^2 * 0.032;
m_hdl   = 0.0088383954;
m_pinB  = 1240 * pi * 0.004^2 * 0.022;
I_pinB  = 0.5 * m_pinB * 0.004^2;
m3_link = 0.042048803; I3_link = 0.00010462162;
m3_tot  = m3_link + m_pinB;
r3_com  = (m3_link*(num.b/2) + m_pinB*num.b)/m3_tot;
IG3_tot = I3_link + m3_link*(r3_com - num.b/2)^2 + I_pinB + m_pinB*(num.b - r3_com)^2;
m4_tot  = 0.03956871; r4_com = num.c/2; IG4_tot = 8.7591844e-5;

dynSym = [m2d r2d mAh m3 r3c IG3 m4 r4c IG4];
dynNum = [m2_disc r2_disc (m_pinA + m_hdl) m3_tot r3_com IG3_tot m4_tot r4_com IG4_tot];

L = sort([num.a num.b num.c num.g]);
fprintf('\nLinks: crank a = %g, coupler b = %g, rocker c = %g, ground g = %g mm\n', [num.a num.b num.c num.g]*1e3);
fprintf('Grashof: s + l = %g <= p + q = %g -> %s; shortest link is the crank -> crank-rocker\n', ...
        (L(1)+L(4))*1e3, (L(2)+L(3))*1e3, string(L(1)+L(4) <= L(2)+L(3)));

fTh4 = matlabFunction(subs(th4_open,  symv, numv), 'Vars', {th2});
fTh4x= matlabFunction(subs(th4_cross, symv, numv), 'Vars', {th2});
fTh3 = matlabFunction(subs(th3_expr,  symv, numv), 'Vars', {th2, th4});
fW   = matlabFunction(subs([w3; w4],  symv, numv), 'Vars', {th2, th3, th4});
fAl  = matlabFunction(subs([al3; al4],symv, numv), 'Vars', {th2, th3, th4});
fMu  = matlabFunction(subs(mu,        symv, numv), 'Vars', {th2});
fAF  = matlabFunction(subs(AF_sym, [symv dynSym], [numv dynNum]), 'Vars', {th2, th3, th4});
fbF  = matlabFunction(subs(bF_sym, [symv dynSym], [numv dynNum]), 'Vars', {th2, th3, th4, w3s, w4s, al3s, al4s});

% Verify circuit against Fusion assembly pose (cr2_assembly.step)
O2f = [124.5 80];
xd  = [0.939187956345515 0.34340352743607];  ad = atan2(xd(2), xd(1));
Af  = [124.5 80] - O2f + 50*[-sin(ad) cos(ad)];
xr  = [0.960399846635695 0.278625437787258]; ar = atan2(xr(2), xr(1));
Bf  = [21.6483788313323 71.8840738846138] - O2f + [5 147]*[cos(ar) sin(ar); -sin(ar) cos(ar)];
th0   = atan2(Af(2), Af(1));
th4F  = atan2(Bf(2), Bf(1) + num.g*1e3);
fprintf('Fusion pose: theta2 = %.4f deg, theta4 = %.4f deg | circuit 1 gives %.4f, circuit 2 gives %.4f deg\n', ...
        rad2deg([th0 th4F fTh4(th0) fTh4x(th0)]));
assert(abs(fTh4(th0) - th4F) < 1e-6, 'Circuit 1 does not match the Fusion assembly');

% Extreme positions, stroke, time ratio, transmission angle limits
lim = double(subs([th4_fold th4_ext muMin muMax], [a b c g], numv(1:4)));
B_fold = [-num.g 0] + num.c*[cos(lim(1)) sin(lim(1))];
B_ext  = [-num.g 0] + num.c*[cos(lim(2)) sin(lim(2))];
th2_fold = atan2(B_fold(2), B_fold(1)) + pi;          % A opposite to B (folded)
th2_ext  = atan2(B_ext(2),  B_ext(1));                % A towards B (extended)
fwd  = mod(th2_ext - th2_fold, 2*pi);                 % crank turn while rocker swings out
TR   = fwd/(2*pi - fwd);
stroke_deg = rad2deg(lim(2)-lim(1));
stroke_arc = num.c * (lim(2) - lim(1)) * 1e3;         % mm
th2_mu40   = acos((num.b^2 + num.c^2 - 2*num.b*num.c*cos(deg2rad(40)) - num.a^2 - num.g^2) / (2*num.a*num.g));
fprintf('Rocker limits: %.3f deg (folded, theta2 = %.2f deg) and %.3f deg (extended, theta2 = %.2f deg)\n', ...
        rad2deg(lim(1)), mod(rad2deg(th2_fold),360), rad2deg(lim(2)), mod(rad2deg(th2_ext),360));
fprintf('Rocker stroke = %.3f deg (arc length at pin B = %.2f mm), time ratio = %.4f (%.2f / %.2f deg)\n', ...
        stroke_deg, stroke_arc, TR, rad2deg(fwd), 360 - rad2deg(fwd));
fprintf('Transmission angle: min %.3f deg (theta2=180), max %.3f deg (theta2=0); mu < 40 deg for theta2 in [%.2f, %.2f] deg\n', ...
        rad2deg(lim(3:4)), rad2deg(th2_mu40), 360 - rad2deg(th2_mu40));

%% 8) Evaluate kinematics + forces over one crank revolution
q2 = deg2rad((0:1:360)');
q4 = fTh4(q2);  q3 = fTh3(q2, q4);
W  = fW(q2.', q3.', q4.').';     AL = fAl(q2.', q3.', q4.').';
MU = fMu(q2);
XF = zeros(numel(q2), 9);
for k = 1:numel(q2)
    XF(k,:) = (fAF(q2(k), q3(k), q4(k)) \ fbF(q2(k), q3(k), q4(k), W(k,1), W(k,2), AL(k,1), AL(k,2))).';
end
F_O2_mag = hypot(XF(:,1), XF(:,2));
F_A_mag  = hypot(XF(:,3), XF(:,4));
F_B_mag  = hypot(XF(:,5), XF(:,6));
F_O4_mag = hypot(XF(:,7), XF(:,8));
T2_req   = XF(:,9);

Tall = table(rad2deg(q2), rad2deg(q3), rad2deg(q4), W(:,1), W(:,2), AL(:,1), AL(:,2), rad2deg(MU), ...
             F_O2_mag, F_A_mag, F_B_mag, F_O4_mag, T2_req, ...
    'VariableNames', {'theta2_deg','theta3_deg','theta4_deg','w3_rad_s','w4_rad_s', ...
                      'alpha3_rad_s2','alpha4_rad_s2','mu_deg', ...
                      'F_O2_N','F_A_N','F_B_N','F_O4_N','T2_Nm'});
writetable(Tall, 'fourbar_symbolic_table.csv');
disp(Tall(1:30:end, :));
fprintf('Peak dynamic forces (1 rev/s): |F_O2|=%.4f N, |F_A|=%.4f N, |F_B|=%.4f N, |F_O4|=%.4f N, |T2|=%.5f N*m\n', ...
        max(F_O2_mag), max(F_A_mag), max(F_B_mag), max(F_O4_mag), max(abs(T2_req)));

% Working space bounding box of pin centers
Ax = num.a*1e3*cos(q2); Ay = num.a*1e3*sin(q2);
Bx = -num.g*1e3 + num.c*1e3*cos(q4); By = num.c*1e3*sin(q4);
fprintf('Joint working space (mm): x in [%.2f, %.2f], y in [%.2f, %.2f]\n', ...
        min([Ax; Bx; -num.g*1e3]), max([Ax; Bx]), min([Ay; By]), max([Ay; By]));

figure('Name', 'Symbolic kinematics');
subplot(2,2,1); plot(rad2deg(q2), rad2deg([q3 q4])); grid on; xlim([0 360]);
xlabel('\theta_2 (deg)'); ylabel('angle (deg)'); legend('\theta_3 coupler', '\theta_4 rocker', 'Location', 'best');
subplot(2,2,2); plot(rad2deg(q2), W); grid on; xlim([0 360]);
xlabel('\theta_2 (deg)'); ylabel('\omega (rad/s)'); legend('\omega_3', '\omega_4', 'Location', 'best');
subplot(2,2,3); plot(rad2deg(q2), AL); grid on; xlim([0 360]);
xlabel('\theta_2 (deg)'); ylabel('\alpha (rad/s^2)'); legend('\alpha_3', '\alpha_4', 'Location', 'best');
subplot(2,2,4); plot(rad2deg(q2), rad2deg(MU)); hold on; yline(40, '--'); grid on; xlim([0 360]);
xlabel('\theta_2 (deg)'); ylabel('\mu (deg)'); legend('transmission angle', '40 deg guideline', 'Location', 'best');

figure('Name', 'Symbolic dynamic forces and input torque');
subplot(2,1,1);
plot(rad2deg(q2), [F_O2_mag F_A_mag F_B_mag F_O4_mag], 'LineWidth', 1.3); grid on; xlim([0 360]);
xlabel('\theta_2 (deg)'); ylabel('Joint reaction force (N)');
legend('|F_{O2}|', '|F_A|', '|F_B|', '|F_{O4}|', 'Location', 'best');
subplot(2,1,2);
plot(rad2deg(q2), T2_req*1e3, 'LineWidth', 1.5); grid on; xlim([0 360]);
xlabel('\theta_2 (deg)'); ylabel('Input crank torque T_2 (mN\cdot m)');

%% 9) Equations for the report
fid = fopen('fourbar_symbolic_equations.txt', 'w');
fprintf(fid, 'P = %s\nQ = %s\nR = %s\n\n', char(P), char(Q), char(R));
fprintf(fid, 'theta4 (circuit 1) = %s\n\n', char(th4_open));
fprintf(fid, 'theta3 = %s\n\n', char(th3_expr));
fprintf(fid, 'w3 = %s\n\nw4 = %s\n\n', char(w3), char(w4));
fprintf(fid, 'mu = %s\n\n', char(mu));
fprintf(fid, 'Newton-Euler 9x9 system A_F * x_F = b_F, x_F = [F12x F12y F32x F32y F43x F43y F14x F14y T2]^T\n');
fprintf(fid, 'LaTeX:\n w_3 = %s\n w_4 = %s\n \\mu = %s\n', latex(w3), latex(w4), latex(mu));
fclose(fid);

%% 10) Compare with the Simscape Multibody run (fourbar_cad_v2.m)
simMat = fullfile('..', 'multibody', 'fourbar_v2_sim.mat');
if ~isfile(simMat), simMat = 'fourbar_v2_sim.mat'; end
if isfile(simMat)
    S = load(simMat);
    s4  = fTh4(S.th);  s3 = fTh3(S.th, s4);
    sw  = fW(S.th.', s3.', s4.').';
    sal = fAl(S.th.', s3.', s4.').';
    XF_sim = zeros(numel(S.t), 9);
    for k = 1:numel(S.t)
        XF_sim(k,:) = (fAF(S.th(k), s3(k), s4(k)) \ fbF(S.th(k), s3(k), s4(k), sw(k,1), sw(k,2), sal(k,1), sal(k,2))).';
    end
    e4  = rad2deg(angle(exp(1i*(S.phi - s4))));
    ew4 = S.w4 - sw(:,2);
    fprintf('\nSymbolic vs Simscape over %d steps:\n', numel(S.t));
    fprintf('  max |d theta4| = %.2e deg, max |d w4| = %.2e rad/s\n', max(abs(e4)), max(abs(ew4)));
    if isfield(S, 'T1sig') && all(isfinite(S.T1sig))
        fprintf('  max |d T2|     = %.2e N*m\n', max(abs(S.T1sig - XF_sim(:,9))));
        fprintf('  max |d |F_O2|| = %.2e N, max |d |F_A|| = %.2e N, max |d |F_B|| = %.2e N, max |d |F_O4|| = %.2e N\n', ...
                max(abs(vecnorm(S.FO2(:,1:2),2,2) - hypot(XF_sim(:,1), XF_sim(:,2)))), ...
                max(abs(vecnorm(S.FA(:,1:2),2,2)  - hypot(XF_sim(:,3), XF_sim(:,4)))), ...
                max(abs(vecnorm(S.FB(:,1:2),2,2)  - hypot(XF_sim(:,5), XF_sim(:,6)))), ...
                max(abs(vecnorm(S.FO4(:,1:2),2,2) - hypot(XF_sim(:,7), XF_sim(:,8)))));
    end
    figure('Name', 'Symbolic vs Simscape');
    subplot(3,1,1); plot(S.t, rad2deg(S.phi), '.', S.t, rad2deg(s4), '-'); grid on;
    xlabel('t (s)'); ylabel('\theta_4 (deg)'); legend('Simscape', 'Symbolic', 'Location', 'best');
    subplot(3,1,2); plot(S.t, S.w4, '.', S.t, sw(:,2), '-'); grid on;
    xlabel('t (s)'); ylabel('\omega_4 (rad/s)'); legend('Simscape', 'Symbolic', 'Location', 'best');
    subplot(3,1,3); plot(S.t, S.T1sig*1e3, '.', S.t, XF_sim(:,9)*1e3, '-'); grid on;
    xlabel('t (s)'); ylabel('T_2 (mN\cdot m)'); legend('Simscape', 'Symbolic', 'Location', 'best');
else
    fprintf('\n(fourbar_v2_sim.mat not found: run fourbar_cad_v2.m first to compare with Simscape)\n');
end
