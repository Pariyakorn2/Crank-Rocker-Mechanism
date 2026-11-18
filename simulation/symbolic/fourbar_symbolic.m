%% fourbar_symbolic.m  -  Symbolic kinematics (Symbolic Math Toolbox)
% Same mechanism and frame as fourbar_cad_v2.m (Fusion "Crank-Rocker_final"):
%   O2 = (0,0) disc axle, O4 = (-g,0) rocker pivot, all angles counter-clockwise from +x.
%   theta2 = crank O2->A, theta3 = coupler A->B, theta4 = rocker O4->B.
% Derives position, velocity, acceleration, transmission angle, limit positions, time ratio and
% working space symbolically, evaluates them for the Fusion CAD dimensions (angles in deg, deg/s,
% deg/s^2; theta2 axis 0..360 deg) and compares with Simscape Multibody.
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

%% 6) Numbers from the Fusion model (SI units)
num   = struct('a', 0.050, 'b', 0.150, 'c', 0.140, 'g', 0.100, 'w2', 2*pi);
symv  = [a b c g w2];
numv  = [num.a num.b num.c num.g num.w2];

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

%% 7) Evaluate kinematics over one crank revolution
q2 = deg2rad((0:1:360)');
q4 = fTh4(q2);  q3 = fTh3(q2, q4);
W  = fW(q2.', q3.', q4.').';     AL = fAl(q2.', q3.', q4.').';
MU = fMu(q2);
% angular quantities in degrees: deg, deg/s, deg/s^2
Tall = table(rad2deg(q2), rad2deg(q3), rad2deg(q4), rad2deg(W(:,1)), rad2deg(W(:,2)), ...
             rad2deg(AL(:,1)), rad2deg(AL(:,2)), rad2deg(MU), ...
    'VariableNames', {'theta2_deg','theta3_deg','theta4_deg','w3_deg_s','w4_deg_s', ...
                      'alpha3_deg_s2','alpha4_deg_s2','mu_deg'});
writetable(Tall, 'fourbar_symbolic_table.csv');
disp(Tall(1:30:end, :));

% Working space bounding box of pin centers
Ax = num.a*1e3*cos(q2); Ay = num.a*1e3*sin(q2);
Bx = -num.g*1e3 + num.c*1e3*cos(q4); By = num.c*1e3*sin(q4);
fprintf('Joint working space (mm): x in [%.2f, %.2f], y in [%.2f, %.2f]\n', ...
        min([Ax; Bx; -num.g*1e3]), max([Ax; Bx]), min([Ay; By]), max([Ay; By]));

% All angular quantities in degrees (deg, deg/s, deg/s^2); theta2 axis 0..360 deg, ticks every 90 deg
figure('Name', 'Symbolic kinematics');
ax1 = subplot(2,2,1); plot(rad2deg(q2), rad2deg([q3 q4]), 'LineWidth', 1.5);
ylabel('angle (deg)'); legend('\theta_3 coupler', '\theta_4 rocker', 'Location', 'best');
ax2 = subplot(2,2,2); plot(rad2deg(q2), rad2deg(W), 'LineWidth', 1.5);
ylabel('\omega (deg/s)'); legend('\omega_3', '\omega_4', 'Location', 'best');
ax3 = subplot(2,2,3); plot(rad2deg(q2), rad2deg(AL), 'LineWidth', 1.5);
ylabel('\alpha (deg/s^2)'); legend('\alpha_3', '\alpha_4', 'Location', 'best');
ax4 = subplot(2,2,4); plot(rad2deg(q2), rad2deg(MU), 'LineWidth', 1.5); hold on; yline(40, '--');
ylabel('\mu (deg)'); legend('transmission angle', '40 deg guideline', 'Location', 'best');
for axk = [ax1 ax2 ax3 ax4]
    grid(axk, 'on'); xlim(axk, [0 360]); xticks(axk, 0:90:360); xlabel(axk, '\theta_2 (deg)');
end

%% 8) Equations for the report
fid = fopen('fourbar_symbolic_equations.txt', 'w');
fprintf(fid, 'P = %s\nQ = %s\nR = %s\n\n', char(P), char(Q), char(R));
fprintf(fid, 'theta4 (circuit 1) = %s\n\n', char(th4_open));
fprintf(fid, 'theta3 = %s\n\n', char(th3_expr));
fprintf(fid, 'w3 = %s\n\nw4 = %s\n\n', char(w3), char(w4));
fprintf(fid, 'mu = %s\n\n', char(mu));
fprintf(fid, 'LaTeX:\n w_3 = %s\n w_4 = %s\n \\mu = %s\n', latex(w3), latex(w4), latex(mu));
fclose(fid);

%% 9) Compare with the Simscape Multibody run (fourbar_cad_v2.m)
simMat = fullfile('..', 'multibody', 'fourbar_v2_sim.mat');
if ~isfile(simMat), simMat = 'fourbar_v2_sim.mat'; end
if isfile(simMat)
    S = load(simMat);
    s4  = fTh4(S.th);  s3 = fTh3(S.th, s4);
    sw  = fW(S.th.', s3.', s4.').';
    sal = fAl(S.th.', s3.', s4.').';
    e4  = rad2deg(angle(exp(1i*(S.phi - s4))));
    ew4 = S.w4 - sw(:,2);
    fprintf('\nSymbolic vs Simscape over %d steps:\n', numel(S.t));
    fprintf('  max |d theta4| = %.2e deg, max |d w4| = %.2e deg/s\n', max(abs(e4)), max(abs(rad2deg(ew4))));
    figure('Name', 'Symbolic vs Simscape');
    subplot(2,1,1); plot(S.t, rad2deg(S.phi), '.', S.t, rad2deg(s4), '-'); grid on;
    xlabel('t (s)'); ylabel('\theta_4 (deg)'); legend('Simscape', 'Symbolic', 'Location', 'best');
    subplot(2,1,2); plot(S.t, rad2deg(S.w4), '.', S.t, rad2deg(sw(:,2)), '-'); grid on;
    xlabel('t (s)'); ylabel('\omega_4 (deg/s)'); legend('Simscape', 'Symbolic', 'Location', 'best');
else
    fprintf('\n(fourbar_v2_sim.mat not found: run fourbar_cad_v2.m first to compare with Simscape)\n');
end
