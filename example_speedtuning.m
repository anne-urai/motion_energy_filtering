% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% If you use the Software for your own research, cite the paper.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
% DEALINGS IN THE SOFTWARE.
%
% Anne Urai and Klaus Wimmer, 2016
% anne.urai@gmail.com / kwimmer@clinic.ub.es

clear all; clc; close all;

sigmas = 0.25:0.05:1; % both sigma_g and sigma_c
spds = [5.5:1:30.5]; % degrees/s
motionenergy = nan(numel(spds), numel(sigmas));

for sigma = sigmas,
    for spd = spds,
        
        % ======================================================= %
        % SIMULATE DATA
        % ======================================================= %
        
        display = struct('frameRate', 60, 'width', 42, 'height', 32, 'dist', 65, ...
            'res', struct('width', 1280, 'height', 1024, 'hz', 60, 'pixelSize', 32), ...
            'center', [640 512]);
        display.ppd     = deg2pix(display, 1);                   % pix
        dots            = setupDots(display);
        dots.direction  = 0;
        dots.coherence  = 1;
        dots.speed      = deg2pix(display, spd); % VARY SPEED
        dots.nDots      = 5; % same as Klaus
        dots.nvar       = 3;
        dots.lifetime   = 1045;
        setup           = struct('nframes', 30);
        
        % OPTION 2. USE ANNE'S FUNCTION
        coordA          = dots_limitedlifetime(setup, display, dots);
        stimMovieA      = coord2stim(display, coordA, -dots.direction);
        
        % ======================================================= %
        % MAKE APPROPIATE FILTERS, APPLY
        % ======================================================= %
        
        cfg.frameRate   = display.frameRate;
        cfg.k           = 60;
        cfg.sigma_c     = sigma;
        cfg.sigma_g     = sigma;
        cfg.ppd         = display.ppd;
        
        % ======================================================= %
        % COMPUTE AND SAVE A BUNCH OF THINGS
        % ======================================================= %
        
        % 1. stimulus spectogram, should turn with increasing speed
        % spectogram(squeeze(sum(stimMovieA,2))', cfg);
        % title(sprintf('%2.2f deg/s', spd / dots.nvar * cfg.frameRate / cfg.ppd));
        
        % uning curve: output of ME filter for each speed
        [f1, f2] = makeSpatialFilters(cfg);
        [g1, g2] = makeTemporalFilters(cfg);
        
        energy = applyFilters(stimMovieA, f1, f2, g1, g2);
        motionenergy(find(spd == spds), find(sigma==sigmas)) = ...
            squeeze(sum(sum(sum(energy(:, :, 20:end)))));
    end
end

% ======================================================= %
% PLOT TUNING CURVES
% ======================================================= %

clf; subplot(3,3,[1 2]);
set(gca, 'NextPlot','replacechildren', 'ColorOrder', viridis(length(sigmas)));
hold on;
x = spds / dots.nvar;

l = plot(x, motionenergy, '-');
[val, idx] = max(motionenergy);
scatter(x(idx), val, 40, viridis(length(sigmas)), 'filled');
xlabel('Speed (deg/s)');
l = legend(l, strread(num2str(sigmas),'%s'), 'Location', 'EastOutside'); legend boxoff;
text(max(x)*1.1, nanmax(motionenergy(:)), 'Sigma');
ylabel('Motion energy (a.u.)');

subplot(3,3,3);
% for each sigma, what's the preferred speed?
[val, idx] = max(motionenergy, [], 1);
scatter(sigmas, x(idx), 40, 'filled');
p = polyfit(x(idx), sigmas, 1);
bestSigma = polyval(p, 3.8);
xlabel('Generated dot speed');
ylabel('Optimal sigma');
lsline;
