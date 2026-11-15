%% RUN_ME.m  -  Crank-rocker four-bar: run everything with one click
% 1) Open this file in MATLAB and press Run (or type RUN_ME in the Command Window).
%
% What happens (about 2 minutes):
%   - fourbar_cad_v2   builds the Simulink / Simscape Multibody model from the Fusion CAD,
%                      simulates 3 crank turns and opens Multibody Explorer (3D animation:
%                      press the Play button at the bottom of that window).
%   - fourbar_symbolic Symbolic Math Toolbox analysis: equations, dynamic forces/torque,
%                      table and plots, compared with the Simscape result.
%   - fourbar_realtime real-time dashboard: position / velocity / acceleration of A, B, P
%                      with a moving cursor, plus the animated kinematic diagram.
% Needs: Simulink, Simscape, Simscape Multibody, Symbolic Math Toolbox.

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
cadScript = fullfile(here, 'simulation', 'multibody', 'fourbar_cad_v2.m');
symScript = fullfile(here, 'simulation', 'symbolic',  'fourbar_symbolic.m');
stepCheck = fullfile(here, 'hardware',   'mechanical','cr2_link4.step');

if ~isfile(stepCheck) || ~isfile(cadScript) || ~isfile(symScript)
    error(['Required repository files are missing relative to RUN_ME.m (%s).\n' ...
           'Make sure the full repository (hardware/, simulation/) is extracted.'], here);
end

need = {'Simulink', 'Simscape', 'Simscape Multibody', 'Symbolic Math Toolbox'};
v = ver;  have = {v.Name};
missing = need(~ismember(need, have));
if ~isempty(missing)
    error('This computer is missing: %s', strjoin(missing, ', '));
end

run(cadScript);
run(symScript);
addpath(fullfile(here, 'simulation', 'symbolic'));
fourbar_realtime;

fprintf(['\nDone.\n' ...
         ' - 3D animation: go to the "Multibody Explorer" window and press Play (bottom bar).\n' ...
         '   Tip: Y up -> Iso -> Fit gives the same view as Fusion.\n' ...
         ' - Simulink diagram: window "fourbar_v2_model".\n' ...
         ' - Real-time dashboard: window "Crank-rocker: real-time symbolic kinematics".\n' ...
         ' - Results: simulation/multibody/fourbar_v2_results.csv, simulation/symbolic/fourbar_symbolic_table.csv\n']);
