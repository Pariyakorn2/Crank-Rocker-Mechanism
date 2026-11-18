function fourbar_realtime(playbackRevPerSec)
%FOURBAR_REALTIME  Real-time symbolic kinematics of the crank-rocker with an animated kinematic diagram.
%   fourbar_realtime        plays at 0.25 rev/s (slow enough to follow)
%   fourbar_realtime(1)     plays in true real time (1 rev/s, the analysis speed)
%
% Same mechanism and frame as fourbar_cad_v2.m / fourbar_symbolic.m (Fusion "Crank-Rocker_final"):
%   O2 = (0,0) disc axle, O4 = (-g,0) rocker pivot, angles counter-clockwise from +x.
%   a = 50 (crank), b = 150 (coupler), c = 140 (rocker), g = 100 (ground) mm, crank 1 rev/s.
% Points tracked: A (crank pin), B (coupler-rocker pin), P (coupler midpoint).
% Position is solved in closed form with the Symbolic Math Toolbox; velocity and acceleration are
% the symbolic time derivatives  v = dr/dtheta2 * w2,  a = d2r/dtheta2^2 * w2^2  (w2 constant).
% Left: x, y and magnitude of position / velocity / acceleration against crank angle, with a cursor
% that moves with the animation. Right: kinematic diagram driven by the same values.

if nargin < 1, playbackRevPerSec = 0.25; end

%% ---- symbolic model ------------------------------------------------------------------------
fprintf('Deriving the kinematics symbolically...\n');
th2 = sym('th2', 'real');                       % (syms cannot be used in a function with nested functions)
a = 50; b = 150; c = 140; g = 100;             % mm
w2 = 2*pi;                                      % rad/s (1 rev/s, counter-clockwise)
O2 = [0; 0];  O4 = [-g; 0];

% Freudenstein form P cos th4 + Q sin th4 + R = 0, open circuit (same as Fusion assembly)
Pq = -2*c*(g + a*cos(th2));
Qq = -2*a*c*sin(th2);
Rq =  a^2 - b^2 + c^2 + g^2 + 2*a*g*cos(th2);
th4 = 2*atan((-Qq + sqrt(Pq^2 + Qq^2 - Rq^2))/(Rq - Pq));
rA  = O2 + a*[cos(th2); sin(th2)];
rB  = O4 + c*[cos(th4); sin(th4)];
th3 = atan2(rB(2) - rA(2), rB(1) - rA(1));
rP  = (rA + rB)/2;                              % coupler midpoint

R   = [rA; rB; rP];                             % 6x1 positions        [mm]
V   = diff(R, th2) * w2;                        % 6x1 velocities       [mm/s]
Acc = diff(R, th2, 2) * w2^2;                   % 6x1 accelerations    [mm/s^2]
ang = [th3; th4];
W   = diff(ang, th2) * w2 * 180/pi;             % w3, w4               [deg/s]
AL  = diff(ang, th2, 2) * w2^2 * 180/pi;        % alpha3, alpha4       [deg/s^2]
d2  = (rA(1) - O4(1))^2 + (rA(2) - O4(2))^2;
mu  = acos((b^2 + c^2 - d2)/(2*b*c));           % transmission angle

fprintf('Generating numeric functions...\n');
F = matlabFunction([R; V; Acc; ang; W; AL; mu], 'Vars', th2);   % 25 rows
% rows: 1-6 positions (Ax Ay Bx By Px Py), 7-12 velocities, 13-18 accelerations,
%       19-20 th3 th4, 21-22 w3 w4, 23-24 al3 al4, 25 mu

%% ---- curves over one revolution + self-check -------------------------------------------------
q  = linspace(0, 2*pi, 721);
S  = F(q);
dq = q(2) - q(1);
vNum = gradient(S(1:6,:), dq) * w2;             % finite-difference check of the symbolic velocity
k  = 5:numel(q)-4;
fprintf('Check: symbolic velocity vs finite difference, max error %.2e mm/s\n', ...
        max(abs(vNum(:,k) - S(7:12,k)), [], 'all'));
fprintf('Rocker %.2f..%.2f deg, |vB| max %.1f mm/s, |aB| max %.0f mm/s^2, mu %.2f..%.2f deg\n', ...
        rad2deg(min(S(20,:))), rad2deg(max(S(20,:))), max(hypot(S(9,:), S(10,:))), ...
        max(hypot(S(15,:), S(16,:))), rad2deg(min(S(25,:))), rad2deg(max(S(25,:))));

%% ---- figure ------------------------------------------------------------------------------------
fig = figure('Name', 'Crank-rocker: real-time symbolic kinematics', 'NumberTitle', 'off', ...
             'Color', 'w', 'Units', 'normalized', 'Position', [0.02 0.05 0.96 0.85]);
tl  = tiledlayout(fig, 3, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
tl.OuterPosition = [0 0.05 1 0.95];             % leave a strip at the bottom for the controls
deg = rad2deg(q);
cX = [0 0.447 0.741];  cY = [0.85 0.325 0.098];  cM = [0.93 0.69 0.13];
names  = {'A (crank pin)', 'B (rocker pin)', 'P (coupler midpoint)'};
qty    = {'position', 'velocity', 'acceleration'};
units  = {'mm', 'mm/s', 'mm/s^2'};
rowOff = [0 6 12];                              % offset of pos / vel / acc rows in S
H = struct('cur', {}, 'mk', {});
for p = 1:3
    for k3 = 1:3
        ax = nexttile(tl, (p-1)*5 + k3);  hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
        ix = rowOff(k3) + 2*p - 1;  iy = ix + 1;
        x = S(ix,:);  y = S(iy,:);
        plot(ax, deg, x, 'Color', cX, 'LineWidth', 1.2);
        plot(ax, deg, y, 'Color', cY, 'LineWidth', 1.2);
        if k3 > 1
            plot(ax, deg, hypot(x, y), 'Color', cM, 'LineWidth', 1.2);
            lg = {'x', 'y', 'magnitude'};
        else
            lg = {'x', 'y'};
        end
        H(p,k3).cur = xline(ax, 0, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
        H(p,k3).mk  = plot(ax, zeros(1, numel(lg)), zeros(1, numel(lg)), 'ko', ...
                           'MarkerFaceColor', 'k', 'MarkerSize', 5, 'HandleVisibility', 'off');
        xlim(ax, [0 360]);  xticks(ax, 0:90:360);
        title(ax, sprintf('%s %s', names{p}(1), qty{k3}));
        ylabel(ax, sprintf('%s (%s)', qty{k3}, units{k3}));
        if p == 3, xlabel(ax, '\theta_2 (deg)'); end
        if p == 1, legend(ax, lg, 'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off'); end
    end
end

% ---- kinematic diagram ----
axD = nexttile(tl, 4, [3 2]);  hold(axD, 'on'); grid(axD, 'on'); box(axD, 'on');
title(axD, 'Kinematic diagram (driven by the symbolic solution)');
xlabel(axD, 'x (mm)'); ylabel(axD, 'y (mm)');
plot(axD, S(1,:), S(2,:), '--', 'Color', [0.6 0.6 0.6]);           % path of A
plot(axD, S(3,:), S(4,:), '--', 'Color', [0.96 0.55 0.15]);        % path of B (arc)
plot(axD, S(5,:), S(6,:), '--', 'Color', [0.45 0.75 0.30]);        % coupler curve of P
plot(axD, [O2(1) O4(1)], [O2(2) O4(2)], '-', 'Color', [0.55 0.72 0.90], 'LineWidth', 3);
plot(axD, [O2(1) O4(1)], [O2(2) O4(2)] - 6, '^', 'MarkerSize', 14, ...
     'MarkerFaceColor', [0.55 0.72 0.90], 'MarkerEdgeColor', 'k');
text(axD, O2(1) + 6, -14, 'O_2');  text(axD, O4(1) + 6, -14, 'O_4');
hCrank = plot(axD, [0 0], [0 0], '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 5);
hCoup  = plot(axD, [0 0], [0 0], '-', 'Color', [0.45 0.75 0.30], 'LineWidth', 5);
hRock  = plot(axD, [0 0], [0 0], '-', 'Color', [0.96 0.55 0.15], 'LineWidth', 5);
hJ     = plot(axD, zeros(1,5), zeros(1,5), 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
hLab   = [text(axD, 0, 0, ' A', 'FontWeight', 'bold'), text(axD, 0, 0, ' B', 'FontWeight', 'bold'), ...
          text(axD, 0, 0, ' P', 'FontWeight', 'bold')];
sv = 0.08;  sa = 0.01;                                             % arrow scales: mm per (mm/s), mm per (mm/s^2)
hV = quiver(axD, zeros(3,1), zeros(3,1), zeros(3,1), zeros(3,1), 0, 'Color', [0 0.3 0.9], 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
hA = quiver(axD, zeros(3,1), zeros(3,1), zeros(3,1), zeros(3,1), 0, 'Color', [0.85 0 0],  'LineWidth', 1.5, 'MaxHeadSize', 0.5);
allx = [reshape(S([1 3 5],:), 1, []), O4(1)];  ally = [reshape(S([2 4 6],:), 1, []), 0];
set(axD, 'XLim', [min(allx) - 40, max(allx) + 40], 'YLim', [min(ally) - 150, max(ally) + 30], ...
         'DataAspectRatio', [1 1 1], 'PlotBoxAspectRatioMode', 'auto');
hTxt = text(axD, 0.02, 0.02, '', 'Units', 'normalized', 'VerticalAlignment', 'bottom', ...
            'FontName', 'Consolas', 'FontSize', 9, 'BackgroundColor', [1 1 1 0.8]);
legend(axD, {'path of A', 'path of B', 'coupler curve of P'}, 'Location', 'northwest');

% ---- controls ----
st.playing = true;  st.theta = 0;  st.speed = playbackRevPerSec;  st.t0 = tic;
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Pause', 'Units', 'normalized', ...
          'Position', [0.62 0.005 0.06 0.035], 'Callback', @onPlay, 'Tag', 'play');
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Replay', 'Units', 'normalized', ...
          'Position', [0.69 0.005 0.06 0.035], 'Callback', @onReplay);
uicontrol(fig, 'Style', 'text', 'String', 'speed (rev/s)', 'Units', 'normalized', 'BackgroundColor', 'w', ...
          'Position', [0.76 0.005 0.06 0.03]);
hSpd = uicontrol(fig, 'Style', 'slider', 'Min', 0.05, 'Max', 1, 'Value', min(max(st.speed, 0.05), 1), ...
          'Units', 'normalized', 'Position', [0.82 0.01 0.07 0.025], 'Callback', @onSpeed);
hCkV = uicontrol(fig, 'Style', 'checkbox', 'String', 'velocity', 'Value', 1, 'Units', 'normalized', ...
          'BackgroundColor', 'w', 'Position', [0.90 0.012 0.045 0.025]);
hCkA = uicontrol(fig, 'Style', 'checkbox', 'String', 'accel.', 'Value', 0, 'Units', 'normalized', ...
          'BackgroundColor', 'w', 'Position', [0.95 0.012 0.045 0.025]);
title(tl, sprintf(['Crank-rocker  a=%g  b=%g  c=%g  g=%g mm,  \\omega_2 = %.0f deg/s (1 rev/s)' ...
                   '   |   symbolic solution, animation slowed to the chosen playback speed'], a, b, c, g, w2*180/pi));

update(0);
tmr = timer('ExecutionMode', 'fixedSpacing', 'Period', 0.03, 'TimerFcn', @tick, ...
            'BusyMode', 'drop', 'Name', 'fourbar_realtime');
fig.CloseRequestFcn = @onClose;
start(tmr);

%% ---- nested callbacks ----------------------------------------------------------------------------
    function tick(~, ~)
        if ~isvalid(fig), return, end
        dt = toc(st.t0);  st.t0 = tic;
        if st.playing
            st.theta = mod(st.theta + 2*pi*st.speed*dt, 2*pi);
            try
                update(st.theta);
            catch err
                stop(tmr);  fprintf(2, 'Animation stopped: %s\n', err.message);
            end
        end
    end

    function update(th)
        s  = F(th);
        dg = rad2deg(th);
        for pp = 1:3
            for kk = 1:3
                ix2 = rowOff(kk) + 2*pp - 1;
                yy = [s(ix2) s(ix2+1)];
                if kk > 1, yy(3) = hypot(s(ix2), s(ix2+1)); end
                H(pp,kk).cur.Value = dg;
                set(H(pp,kk).mk, 'XData', dg*ones(size(yy)), 'YData', yy);
            end
        end
        A = s(1:2);  B = s(3:4);  Pm = s(5:6);
        set(hCrank, 'XData', [O2(1) A(1)], 'YData', [O2(2) A(2)]);
        set(hCoup,  'XData', [A(1) B(1)],  'YData', [A(2) B(2)]);
        set(hRock,  'XData', [O4(1) B(1)], 'YData', [O4(2) B(2)]);
        set(hJ, 'XData', [O2(1) O4(1) A(1) B(1) Pm(1)], 'YData', [O2(2) O4(2) A(2) B(2) Pm(2)]);
        pts = [A B Pm];
        for ii = 1:3, hLab(ii).Position = [pts(:,ii)' 0]; end
        vv = reshape(s(7:12), 2, 3);  aa = reshape(s(13:18), 2, 3);
        set(hV, 'XData', pts(1,:)', 'YData', pts(2,:)', 'UData', sv*vv(1,:)', 'VData', sv*vv(2,:)', ...
                'Visible', onoff(hCkV.Value));
        set(hA, 'XData', pts(1,:)', 'YData', pts(2,:)', 'UData', sa*aa(1,:)', 'VData', sa*aa(2,:)', ...
                'Visible', onoff(hCkA.Value));
        hTxt.String = sprintf(['\\theta_2 = %6.1f deg   \\theta_3 = %6.1f   \\theta_4 = %6.1f deg\n' ...
                               '\\omega_3 = %6.1f   \\omega_4 = %6.1f deg/s\n' ...
                               '\\alpha_3 = %6.0f   \\alpha_4 = %6.0f deg/s^2\n' ...
                               '|v_B| = %6.1f mm/s   |a_B| = %7.1f mm/s^2\n' ...
                               '\\mu = %5.1f deg (transmission angle)'], ...
                              dg, rad2deg(mod(s(19), 2*pi)), rad2deg(s(20)), s(21), s(22), s(23), s(24), ...
                              hypot(s(9), s(10)), hypot(s(15), s(16)), rad2deg(s(25)));
        drawnow limitrate
    end

    function onPlay(src, ~)
        st.playing = ~st.playing;  st.t0 = tic;
        if st.playing, src.String = 'Pause'; else, src.String = 'Play'; end
    end

    function onReplay(~, ~)
        st.theta = 0;  st.playing = true;  st.t0 = tic;
        set(findobj(fig, 'Tag', 'play'), 'String', 'Pause');
        update(0);
    end

    function onSpeed(src, ~)
        st.speed = src.Value;
    end

    function onClose(~, ~)
        try, stop(tmr); delete(tmr); catch, end
        delete(fig);
    end
end

function s = onoff(v)
    if v, s = 'on'; else, s = 'off'; end
end
