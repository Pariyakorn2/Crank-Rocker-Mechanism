%% record_simscape_avi.m
% Records the Simscape Multibody 3D animation (Mechanics Explorer) to an .avi file.
% Usage:
%   1) In MATLAB, set the current folder to simulation/multibody (where fourbar_v2_model.slx is),
%      or run RUN_ME.m once first so the model exists.
%   2) Run this file. Output: crank_rocker.avi (in the current folder).
% Needs: Simscape Multibody (function smwritevideo, R2021b or newer).
% If smwritevideo is not available: in Mechanics Explorer use the Record button
% (camera icon in the toolbar) and save the video as AVI.

mdl = 'fourbar_v2_model';
out = fullfile(pwd, 'crank_rocker.avi');

if ~bdIsLoaded(mdl), load_system(mdl); end

% Simulate once so Mechanics Explorer holds the animation data
sim(mdl);

% Write the 3D animation (real time, 30 fps, Motion-JPEG AVI)
smwritevideo(mdl, out, ...
    'PlaybackSpeedRatio', 1, ...
    'FrameRate',     30, ...
    'VideoFormat',    'motion jpeg avi');

fprintf('Saved: %s\n', out);
