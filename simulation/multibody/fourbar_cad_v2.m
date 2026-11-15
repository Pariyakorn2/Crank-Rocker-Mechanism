%% fourbar_cad_v2.m  -  Crank-rocker (Fusion 360 "Crank-Rocker_final") in Simscape Multibody
% Real CAD (STEP exported from Fusion) placed exactly as in the Fusion assembly, closed-loop
% kinematics driven at constant crank speed, verified against the closed-form position solution.
% Needs in the current folder: cr2_link2.step cr2_link3.step cr2_link4.step cr2_handle.step cr2_ground.step
% Then run fourbar_symbolic.m for the symbolic (Symbolic Math Toolbox) analysis and comparison.
%
% Link lengths measured from the hole centres in the Fusion parts (mm):
%   crank  a = 50   (Link 4 disc: axle hole (0,0) -> crank-pin hole (0,50))
%   coupler b = 150 (Link 3 holes at y = 7 and 157)
%   rocker c = 140  (Link 2 holes at y = 7 and 147)
%   ground g = 100  (Ground: disc axle (124.5,80) -> rocker pivot (24.5,80))
% World frame = Fusion assembly frame shifted so the disc axle O2 is at the origin
% (x right, y up, z out of the mechanism plane). O4 is therefore at (-100, 0).
clear; clc; close all;
mdl = 'fourbar_v2_model';
dirScript = fileparts(mfilename('fullpath'));    % folder of this script (works from any Current Folder)
if isempty(dirScript), dirScript = pwd; end
dirCAD = fullfile(dirScript, '..', '..', 'hardware', 'mechanical');
if isfolder(dirCAD) && isfile(fullfile(dirCAD, 'cr2_link4.step'))
    dirSTEP = char(java.io.File(dirCAD).getCanonicalPath());
else
    dirSTEP = dirScript;
end
stepFiles = {'cr2_link2.step','cr2_link3.step','cr2_link4.step','cr2_handle.step','cr2_ground.step'};
if ~all(cellfun(@(f) isfile(fullfile(dirSTEP, f)), stepFiles))
    error(['The STEP files were not found in (%s).\n' ...
           'If you opened it from inside a .zip: extract all files first, ' ...
           'then run RUN_ME.m from the repository root.'], dirSTEP);
end
cd(dirScript);                                  % results are written next to the script

%% ---- parameters (SI) ----------------------------------------------------------------------
a = 0.050; b = 0.150; c = 0.140; g = 0.100;     % crank, coupler, rocker, ground
gx = -g;                                        % rocker pivot O4 at x = -g
omega = 2*pi;                                   % crank speed [rad/s] (1 rev/s, counter-clockwise)
tStop = 3;                                      % [s] = 3 revolutions
rho = 1240;                                     % PLA filament [kg/m^3] (solid, 100% infill)
mm  = 1e-3;
O2f = [124.5 80];                               % disc axle in the Fusion frame [mm]

% Part placements read from cr2_assembly.step (Fusion frame, mm): origin, x-axis direction
F.disc    = struct('p', [124.5            80               -7.5], 'x', [0.939187956345515  0.34340352743607 ]);
F.coupler = struct('p', [110.099015571377 118.814977888157  2.5], 'x', [0.583310537674698  0.812249233079147]);
F.rocker  = struct('p', [21.6483788313323 71.8840738846138 12.5], 'x', [0.960399846635695  0.278625437787258]);
F.handle  = struct('p', [107.329823628197 126.959397817276 12.5], 'x', [0.812249233079144 -0.583310537674695]);
F.ground  = struct('p', [0 0 0],                                 'x', [1 0]);

% crank angle of the saved Fusion pose: crank-pin hole (0,50) of the disc, seen from O2
ad   = atan2(F.disc.x(2), F.disc.x(1));
Af   = F.disc.p(1:2) - O2f + 50*[-sin(ad) cos(ad)];
th0  = atan2(Af(2), Af(1));                                  % = 110.08 deg
[phi0, cp0, A0, B0] = solveRocker(th0, a, b, c, gx);         % same circuit as Fusion

col.rocker  = [0.96 0.55 0.15];                 % orange     (Link 2)
col.coupler = [0.45 0.75 0.30];                 % green      (Link 3)
col.disc    = [0.62 0.62 0.62];                 % grey       (Link 4, large disc)
col.ground  = [0.55 0.72 0.90];                 % light blue (Ground)
col.handle  = [0.08 0.08 0.08];                 % black      (Handle, small disc)
col.pin     = [0.45 0.45 0.50];                 % pins (visual, PLA mass)

%% ---- model --------------------------------------------------------------------------------
if bdIsLoaded(mdl), close_system(mdl, 0); end
if isfile([mdl '.slx']), delete([mdl '.slx']); end   % rebuild from scratch
new_system(mdl); open_system(mdl);
notes('reset');

lib.world  = 'sm_lib/Frames and Transforms/World Frame';
lib.mech   = 'sm_lib/Utilities/Mechanism Configuration';
lib.solver = 'nesl_utility/Solver Configuration';
lib.rev    = 'sm_lib/Joints/Revolute Joint';
lib.rt     = 'sm_lib/Frames and Transforms/Rigid Transform';
lib.file   = 'sm_lib/Body Elements/File Solid';
lib.cyl    = 'sm_lib/Body Elements/Cylindrical Solid';
lib.s2ps   = 'nesl_utility/Simulink-PS Converter';
lib.ramp   = 'simulink/Sources/Ramp';
lib.const  = 'simulink/Sources/Constant';

W  = addb(mdl, lib.world,  'World', 30, 200);
MC = addb(mdl, lib.mech,   'MechConfig', 30, 300);
SC = addb(mdl, lib.solver, 'Solver', 30, 380);
link2(mdl, port(W,'R',1), port(MC,'any',1));
link2(mdl, port(W,'R',1), port(SC,'any',1));

% joints (all about +Z). Joint frames lie in z = 0; the parts carry their own z offsets.
J1 = addb(mdl, lib.rev, 'J1_crank_axle',   200, 100);   % O2: ground - disc
J2 = addb(mdl, lib.rev, 'J2_pin_A',        500, 100);   % A : disc - coupler
J3 = addb(mdl, lib.rev, 'J3_pin_B',        800, 100);   % B : coupler - rocker
J4 = addb(mdl, lib.rev, 'J4_rocker_pivot', 800, 300);   % O4: rocker - ground

% link frames: crank x-axis O2->A, coupler x-axis A->B, rocker x-axis O4->B
RTa = addrt(mdl, lib, 'RT_crank_a',   350, 100, [ a  0 0], 0);   % crank pin A
RTb = addrt(mdl, lib, 'RT_coupler_b', 650, 100, [ b  0 0], 0);   % coupler pin B
RTc = addrt(mdl, lib, 'RT_rocker_c',  950, 200, [-c  0 0], 0);   % rocker: B -> O4
RTg = addrt(mdl, lib, 'RT_ground_g',  650, 300, [gx  0 0], 0);   % ground: O2 -> O4

link2(mdl, port(W,'R',1),  port(J1,'L',1));
link2(mdl, port(J1,'R',1), port(RTa,'L',1));
link2(mdl, port(RTa,'R',1),port(J2,'L',1));
link2(mdl, port(J2,'R',1), port(RTb,'L',1));
link2(mdl, port(RTb,'R',1),port(J3,'L',1));
link2(mdl, port(J3,'R',1), port(RTc,'L',1));
link2(mdl, port(RTc,'R',1),port(J4,'R',1));
link2(mdl, port(W,'R',1),  port(RTg,'L',1));
link2(mdl, port(RTg,'R',1),port(J4,'L',1));      % closes the loop

% ---- CAD bodies ------------------------------------------------------------------------------
% Each part's offset from the link frame it rides on is computed from its Fusion placement and
% the link frame pose in the same (Fusion) configuration, so the CAD sits exactly as in Fusion.
% Fusion heights (z, mm): axle boss -17.5..-7.5 | disc -7.5..2.5 | coupler 2.5..12.5 |
%                         rocker 12.5..22.5 | handle on pin A 12.5..37.5 | pivot boss 22.5..32.5
Dk = addcad(mdl, lib, 'Link4_disc',    'cr2_link4.step',  dirSTEP, rho, col.disc);
Hd = addcad(mdl, lib, 'Handle',        'cr2_handle.step', dirSTEP, rho, col.handle);
Cp = addcad(mdl, lib, 'Link3_coupler', 'cr2_link3.step',  dirSTEP, rho, col.coupler);
Rk = addcad(mdl, lib, 'Link2_rocker',  'cr2_link2.step',  dirSTEP, rho, col.rocker);
Gr = addcad(mdl, lib, 'Ground',        'cr2_ground.step', dirSTEP, rho, col.ground);

[v, r] = relPose(F.disc,    [0 0],      th0, O2f);  T1 = addrt(mdl, lib, 'T_disc',    350,  20, v*mm, r);
[v, r] = relPose(F.handle,  A0*1e3,     th0, O2f);  T2 = addrt(mdl, lib, 'T_handle',  350, -40, v*mm, r);
[v, r] = relPose(F.coupler, A0*1e3,     cp0, O2f);  T3 = addrt(mdl, lib, 'T_coupler', 650,  20, v*mm, r);
[v, r] = relPose(F.rocker,  B0*1e3,     phi0, O2f); T4 = addrt(mdl, lib, 'T_rocker',  950,  20, v*mm, r);
[v, r] = relPose(F.ground,  [0 0],      0,    O2f); T5 = addrt(mdl, lib, 'T_ground',  350, 450, v*mm, r);

link2(mdl, port(J1,'R',1),  port(T1,'L',1));  link2(mdl, port(T1,'R',1), port(Dk,'any',1));
% The handle is a grip: it turns freely on pin A (joint J5) and the hand keeps it from rotating,
% so J5 is driven with q5 = -(theta2 - theta2_0): the handle keeps its Fusion orientation in space.
J5 = addb(mdl, lib.rev, 'J5_handle_grip', 350, -100);
link2(mdl, port(RTa,'R',1), port(J5,'L',1));
link2(mdl, port(J5,'R',1),  port(T2,'L',1));  link2(mdl, port(T2,'R',1), port(Hd,'any',1));
link2(mdl, port(J2,'R',1),  port(T3,'L',1));  link2(mdl, port(T3,'R',1), port(Cp,'any',1));
link2(mdl, port(J3,'R',1),  port(T4,'L',1));  link2(mdl, port(T4,'R',1), port(Rk,'any',1));
link2(mdl, port(W,'R',1),   port(T5,'L',1));  link2(mdl, port(T5,'R',1), port(Gr,'any',1));

% steel pins (visual only, r = 4 mm in the 4.05 mm holes): centre z and length in mm
pinO2 = addpin(mdl, lib, 'Pin_O2', [0 0 -7.5]*mm, 22*mm, 4*mm, col.pin, 450, -160);
pinA  = addpin(mdl, lib, 'Pin_A',  [0 0  7.5]*mm, 32*mm, 4*mm, col.pin, 450, -100);
pinB  = addpin(mdl, lib, 'Pin_B',  [0 0 12.5]*mm, 22*mm, 4*mm, col.pin, 750, -100);
pinO  = addpin(mdl, lib, 'Pin_O4', [0 0 22.5]*mm, 22*mm, 4*mm, col.pin, 550, 450);
link2(mdl, port(W,'R',1),   port(pinO2.T,'L',1)); link2(mdl, port(pinO2.T,'R',1), port(pinO2.S,'any',1));
link2(mdl, port(RTa,'R',1), port(pinA.T,'L',1));  link2(mdl, port(pinA.T,'R',1),  port(pinA.S,'any',1));
link2(mdl, port(RTb,'R',1), port(pinB.T,'L',1));  link2(mdl, port(pinB.T,'R',1),  port(pinB.S,'any',1));
link2(mdl, port(RTg,'R',1), port(pinO.T,'L',1));  link2(mdl, port(pinO.T,'R',1),  port(pinO.S,'any',1));

% ---- drive: crank angle = th0 + omega*t, with its derivatives supplied ------------------------
ramp = addb(mdl, lib.ramp, 'Crank_angle', 30, 40);
set_param(ramp, 'slope', num2str(omega, 12), 'start', '0', 'X0', num2str(th0, 12));
conv = addb(mdl, lib.s2ps, 'Angle_to_PS', 120, 40);
setp(conv, 'Unit', 'rad');
setp(conv, 'FilteringAndDerivatives', 'provide');
setp(conv, 'UdotUserProvided', '2');
cw = addb(mdl, lib.const, 'Crank_speed', 30, 100);  set_param(cw, 'Value', num2str(omega, 12));
ca = addb(mdl, lib.const, 'Crank_accel', 30, 150);  set_param(ca, 'Value', '0');
ph = get_param(conv, 'PortHandles');
add_line(mdl, get_param(ramp,'PortHandles').Outport(1), ph.Inport(1));
add_line(mdl, get_param(cw,'PortHandles').Outport(1),   ph.Inport(2));
add_line(mdl, get_param(ca,'PortHandles').Outport(1),   ph.Inport(3));
setEnum(J1, 'MotionActuationMode', 'Input');
setEnum(J1, 'TorqueActuationMode', 'Computed');
link2(mdl, port(conv,'R',1), port(J1,'L',2));

ramp5 = addb(mdl, lib.ramp, 'Grip_angle', 30, -120);
set_param(ramp5, 'slope', num2str(-omega, 12), 'start', '0', 'X0', '0');
conv5 = addb(mdl, lib.s2ps, 'Grip_to_PS', 120, -120);
setp(conv5, 'Unit', 'rad');
setp(conv5, 'FilteringAndDerivatives', 'provide');
setp(conv5, 'UdotUserProvided', '2');
cw5 = addb(mdl, lib.const, 'Grip_speed', 30, -60);  set_param(cw5, 'Value', num2str(-omega, 12));
ca5 = addb(mdl, lib.const, 'Grip_accel', 30, -20);  set_param(ca5, 'Value', '0');
ph5 = get_param(conv5, 'PortHandles');
add_line(mdl, get_param(ramp5,'PortHandles').Outport(1), ph5.Inport(1));
add_line(mdl, get_param(cw5,'PortHandles').Outport(1),   ph5.Inport(2));
add_line(mdl, get_param(ca5,'PortHandles').Outport(1),   ph5.Inport(3));
setEnum(J5, 'MotionActuationMode', 'Input');
setEnum(J5, 'TorqueActuationMode', 'Computed');
link2(mdl, port(conv5,'R',1), port(J5,'L',2));

% ---- joint force & torque sensing -----------------------------------------------------------
lib.ps2s = 'nesl_utility/PS-Simulink Converter';
lib.out  = 'built-in/Outport';
set_param(J1, 'SenseTorqueForce', 'on', 'SenseConstraintForce', 'on');
set_param(J2, 'SenseConstraintForce', 'on');
set_param(J3, 'SenseConstraintForce', 'on');
set_param(J4, 'SenseConstraintForce', 'on');

addSenseOut(mdl, lib, J1, 2, 'T1_torque',  'N*m', 200, -20);
addSenseOut(mdl, lib, J1, 3, 'F_O2',       'N',   200,  30);
addSenseOut(mdl, lib, J2, 2, 'F_A',        'N',   500,  30);
addSenseOut(mdl, lib, J3, 2, 'F_B',        'N',   800,  30);
addSenseOut(mdl, lib, J4, 2, 'F_O4',       'N',   800, 380);

% ---- assembly targets (same closure circuit as the Fusion pose) ------------------------------
setTarget(J1, th0); setTarget(J2, cp0 - th0); setTarget(J4, phi0);

set_param(mdl, 'Solver', 'ode23t', 'StopTime', num2str(tStop), 'MaxStep', '1e-3', ...
               'RelTol', '1e-6', 'AbsTol', '1e-8');
set_param(mdl, 'SimscapeLogType', 'all', 'SimscapeLogName', 'simlog');
try, Simulink.BlockDiagram.arrangeSystem(mdl); catch, end   % tidy block layout
save_system(mdl);
notes('print');

%% ---- simulate and verify against the closed-form solution -----------------------------------
simOut = sim(mdl);
sl  = simOut.simlog;
q4  = sl.J4_rocker_pivot.Rz.q.series;
t   = q4.time;  phi = q4.values('rad');
w4  = sl.J4_rocker_pivot.Rz.w.series.values('rad/s');
q2  = sl.J2_pin_A.Rz.q.series.values('rad');
th  = th0 + omega*t;
psi = th + q2;                                              % coupler angle (from +x)
phiRef = arrayfun(@(x) solveRocker(x, a, b, c, gx), th);
err = rad2deg(phi - phiRef);
T1sig = getOutData(simOut, 'T1_torque', t);
FO2   = getOutData(simOut, 'F_O2', t);
FA    = getOutData(simOut, 'F_A',  t);
FB    = getOutData(simOut, 'F_B',  t);
FO4   = getOutData(simOut, 'F_O4', t);
q5  = sl.J5_handle_grip.Rz.q.series.values('rad');
gripDev = rad2deg(max(abs(th + q5 - th0)));                % handle rotation in space (should be 0)

fprintf('\nCrank-rocker  a=%g b=%g c=%g g=%g mm, crank speed %.4f rad/s\n', [a b c g]*1e3, omega);
fprintf('Rocker angle : Simscape %.3f .. %.3f deg | closed form %.3f .. %.3f deg\n', ...
        rad2deg([min(phi) max(phi) min(phiRef) max(phiRef)]));
fprintf('Max |Simscape - closed form| = %.2e deg over %d logged steps\n', max(abs(err)), numel(t));
fprintf('Handle rotation in space during the run: %.2e deg (held still by the hand)\n', gripDev);
if all(isfinite(T1sig))
    fprintf('Crank drive torque (PLA, 1 rev/s, no friction): %.3e .. %.3e N*m\n', min(T1sig), max(T1sig));
    fprintf('Peak joint reaction forces: |F_O2|=%.4f N, |F_A|=%.4f N, |F_B|=%.4f N, |F_O4|=%.4f N\n', ...
            max(vecnorm(FO2(:,1:2),2,2)), max(vecnorm(FA(:,1:2),2,2)), ...
            max(vecnorm(FB(:,1:2),2,2)),  max(vecnorm(FO4(:,1:2),2,2)));
end

rel = (0:30:360)';  tt = rel/360/(omega/(2*pi));
T = table(mod(rad2deg(th0)+rel, 360), ...
          rad2deg(interp1(t, phi, tt, 'spline')), ...
          rad2deg(arrayfun(@(x) solveRocker(th0+deg2rad(x), a, b, c, gx), rel)), ...
          interp1(t, w4, tt, 'spline'), ...
          interp1(t, T1sig, tt, 'spline'), ...
          'VariableNames', {'crank_deg','rocker_sim_deg','rocker_theory_deg','rocker_speed_sim_rad_s','crank_torque_sim_Nm'});
disp(T); writetable(T, 'fourbar_v2_results.csv');
save('fourbar_v2_sim.mat', 't', 'th', 'phi', 'psi', 'w4', 'T1sig', 'FO2', 'FA', 'FB', 'FO4', 'a', 'b', 'c', 'g', 'omega', 'th0');

figure('Name', 'Simscape vs closed form');
subplot(2,1,1); plot(rad2deg(th), rad2deg(phi), '.', rad2deg(th), rad2deg(phiRef), '-');
grid on; xlabel('Crank angle \theta_2 (deg, continuous)'); ylabel('Rocker angle \theta_4 (deg)');
legend('Simscape', 'Closed form', 'Location', 'best');
subplot(2,1,2); plot(t, err); grid on; xlabel('t (s)'); ylabel('Difference (deg)');

%% ============================ local functions =============================================
function [phi, cp, A, B] = solveRocker(theta, a, b, c, gx)
% closed-form position solution, circuit with B above the line A-O4 (as in Fusion)
    A = a*[cos(theta) sin(theta)];
    d = [gx 0] - A;  D = norm(d);  e = d/D;
    k = (b^2 - c^2 + D^2)/(2*D);  h = sqrt(b^2 - k^2);
    M = A + k*e;
    B1 = M + h*[-e(2) e(1)];  B2 = M - h*[-e(2) e(1)];
    if B1(2) >= B2(2), B = B1; else, B = B2; end
    phi = atan2(B(2), B(1) - gx);
    cp  = atan2(B(2) - A(2), B(1) - A(1));
end

function [v, rotDeg] = relPose(P, frameOrigin_mm, frameAngle, O2f)
% offset [mm] and z-rotation [deg] of a Fusion-placed part relative to a link frame that sits at
% frameOrigin_mm (world, mm) with x-axis at frameAngle (rad), world = Fusion - O2f
    al = atan2(P.x(2), P.x(1));
    d  = P.p(1:2) - O2f - frameOrigin_mm;
    R  = [cos(frameAngle) sin(frameAngle); -sin(frameAngle) cos(frameAngle)];  % world -> frame
    v  = [(R*d(:))' P.p(3)];
    rotDeg = mod(rad2deg(al - frameAngle) + 180, 360) - 180;
end

function h = addb(mdl, lib, name, x, y)
    h = add_block(lib, [mdl '/' name], 'Position', [x y x+70 y+40]);
end

function h = addrt(mdl, lib, name, x, y, v, rotZdeg)
    h = addb(mdl, lib.rt, name, x, y);
    setp(h, 'TranslationMethod', 'Cartesian');
    setp(h, 'TranslationCartesianOffset', sprintf('[%.12g %.12g %.12g]', v));
    setp(h, 'TranslationCartesianOffsetUnits', 'm');
    if abs(rotZdeg) > 1e-9
        setp(h, 'RotationMethod', 'StandardAxis');
        setp(h, 'RotationStandardAxis', '+Z');
        setp(h, 'RotationAngle', num2str(rotZdeg, 12));
        setp(h, 'RotationAngleUnits', 'deg');
    end
end

function h = addcad(mdl, lib, name, file, dirSTEP, rho, rgb)
    h = addb(mdl, lib.file, name, 600, 600);
    f = fullfile(dirSTEP, file);
    if ~isfile(f), error('Missing STEP file: %s', f); end
    setp(h, 'ExtGeomFileName', f);
    setp(h, 'BasedOnType', 'Density');
    setp(h, 'InertiaType', 'CalculateFromGeometry');
    setp(h, 'Density', num2str(rho));
    setColor(h, rgb);
end

function s = addpin(mdl, lib, name, off, len, rad, rgb, x, y)
    s.T = addrt(mdl, lib, ['T_' name], x, y, off, 0);
    s.S = addb(mdl, lib.cyl, name, x+150, y);
    setp(s.S, 'CylinderRadius', num2str(rad));
    setp(s.S, 'CylinderLength', num2str(len));
    setp(s.S, 'Density', '1240');
    setColor(s.S, rgb);
end

function setColor(h, rgb)
    dp = get_param(h, 'DialogParameters');  n = fieldnames(dp);
    k = find(strcmp(n, 'Color'), 1);
    if isempty(k), k = find(contains(n, 'Color') & ~contains(n, {'Specular','Ambient','Emissive','Shininess'}), 1); end
    try
        vp = n(strcmp(n, 'VisualProperties')); if ~isempty(vp), set_param(h, vp{1}, 'Simple'); end
        set_param(h, n{k}, sprintf('[%.3f %.3f %.3f]', rgb));
    catch
        notes('add', sprintf('%s: set colour [%.2f %.2f %.2f] manually', get_param(h,'Name'), rgb));
    end
end

function P = port(h, side, idx)
    ph = get_param(h, 'PortHandles');
    switch side
        case 'L', P = ph.LConn(idx);
        case 'R', P = ph.RConn(idx);
        otherwise
            if ~isempty(ph.LConn), P = ph.LConn(idx); else, P = ph.RConn(idx); end
    end
end

function link2(mdl, a, b), add_line(mdl, a, b, 'autorouting', 'on'); end

function ok = setp(h, name, val)
    try, set_param(h, name, val); ok = true;
    catch ME
        ok = false;
        notes('add', sprintf('%s.%s = %s  (%s)', get_param(h,'Name'), name, val, ME.message(1:min(end,80))));
    end
end

function setEnum(h, namePat, valPat)
    dp = get_param(h, 'DialogParameters');  names = fieldnames(dp);
    hit = names(contains(names, namePat));
    for i = 1:numel(hit)
        en = dp.(hit{i}).Enum;  j = find(contains(lower(en), lower(valPat)), 1);
        if ~isempty(j), set_param(h, hit{i}, en{j}); return, end
    end
    notes('add', sprintf('%s: set %s to %s', get_param(h,'Name'), namePat, valPat));
end

function setTarget(h, q)
    dp = get_param(h, 'DialogParameters');  n = fieldnames(dp);
    s = n(contains(n, 'PositionTargetSpecify'));
    v = n(startsWith(n, 'PositionTargetValue') & ~contains(n, 'Unit') & ~contains(n, 'Priority'));
    u = n(contains(n, 'PositionTargetValueUnit'));
    try
        set_param(h, s{1}, 'on');
        if ~isempty(u), set_param(h, u{1}, 'rad'); end
        set_param(h, v{1}, num2str(q, 12));
    catch
        notes('add', sprintf('%s: set position target %.4f rad', get_param(h,'Name'), q));
    end
end

function notes(cmd, msg)
    persistent list
    if isempty(list), list = {}; end
    switch cmd
        case 'reset', list = {};
        case 'add',   list{end+1} = msg; %#ok<AGROW>
        case 'print'
            if isempty(list), fprintf('Model built, no warnings.\n'); return, end
            fprintf('\nWARNINGS (%d):\n', numel(list));
            for i = 1:numel(list), fprintf('  %d. %s\n', i, list{i}); end
    end
end

function addSenseOut(mdl, lib, jBlk, rIdx, outName, unitStr, x, y)
    cv = addb(mdl, lib.ps2s, ['Conv_' outName], x + 90, y);
    setp(cv, 'Unit', unitStr);
    op = addb(mdl, lib.out, outName, x + 190, y);
    add_line(mdl, port(jBlk, 'R', rIdx), port(cv, 'L', 1), 'autorouting', 'on');
    add_line(mdl, get_param(cv, 'PortHandles').Outport(1), get_param(op, 'PortHandles').Inport(1), 'autorouting', 'on');
end

function v = getOutData(simOut, outName, tRef)
    ds = simOut.yout;
    sig = [];
    for k = 1:ds.numElements
        el = ds{k};
        bp = el.BlockPath.getBlock(1);
        if endsWith(bp, ['/' outName])
            sig = el; break;
        end
    end
    if isempty(sig), error('Output signal %s not found in simOut.yout', outName); end
    ts  = sig.Values;
    d   = squeeze(ts.Data);
    if size(d, 1) ~= numel(ts.Time) && size(d, 2) == numel(ts.Time), d = d.'; end
    v   = interp1(ts.Time, double(d), tRef, 'linear', 'extrap');
end
