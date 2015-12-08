%% 3. visualize results

clear all; close all; clc;
cd ~/Data/MotionEnergy/

subjects = 15;
theta = 1:360;

for sj = unique(subjects),
    
    files = dir(sprintf('motionenergy_P%d_s*.mat', sj));
    for f = 1, %:length(files),
        
        load(files(f).name);
        allenergy = motionenergy;
        trialscalculated = find(~isnan(nanmean(nanmean(allenergy.one, 3), 2)));
        
        for trial = 1; %unique(trialscalculated)',
            close all;
            figure; set(gcf, 'DefaultAxesFontSize', 8);
            
            for int = 1:2,
                switch int
                    case 1
                        motionenergy = squeeze(allenergy.one(trial, :, :))';
                        realdir = realdirection.one(trial);
                    case 2
                        motionenergy = squeeze(allenergy.two(trial, :, :))';
                        realdir = realdirection.two(trial);
                end
                
                subplot(3, 3, 1+(int-1)*3);
                imagesc((1:size(motionenergy,1))*(1/60), theta, motionenergy');
                colormap bone;
                set(gca, 'TickDir', 'out');
                xlabel('seconds'); ylabel('direction (degree)');
                % cb = colorbar; cb.Label.String = 'motion energy (a.u.)';
                
                subplot(3, 3, 2+(int-1)*3);
                plotmePolar = {};
                for f = 1:size(motionenergy,1),
                    plotmePolar = [plotmePolar {theta*(pi/180)}];
                    plotmePolar = [plotmePolar {motionenergy(f,:)}];
                    plotCol(f).TTickSign = '+';
                    
                end
                
                hold on;
                try
                    set(0,'DefaultAxesColorOrder', parula(size(motionenergy,1))); % plot each line more yellow
                catch
                    set(0,'DefaultAxesColorOrder', hot(size(motionenergy,1))); % plot each line more yellow
                end
                m = mmpolar(plotmePolar{:}, plotCol, 'fontsize', 4, 'border', 'off', 'RGridLineWidth', 0.1, 'TGridLineWidth', 0.1); % plot the motion energy over time

                hold on;
                % add the real direction
                mmpolar(realdir*(pi/180)*ones(1, 2), ...
                    [min(motionenergy(:)) max(motionenergy(:))], 'r-');
                
                subplot(3, 3, 3+(int-1)*3);
                
                hold on;
                % compute population vector for each frame
                circmean = circ_mean(repmat(theta*(pi/180), size(motionenergy, 1), 1),...
                    motionenergy, 2) / (pi/180);
                circmean(circmean<0) = circmean(circmean<0) + 360;

                [~, ind] = max(motionenergy, [], 2);
                hold on;
                % plot those 2 values
                plot((1:size(motionenergy,1))*(1/60), theta(ind), '.b'); % maximum motion energy
                plot((1:size(motionenergy,1))*(1/60), circmean, '.g'); % population vector - mean in vector space

                 ylabel('Direction (deg)');
                r = refline(0, realdir); set(r, 'Color', 'r');
                axis tight; box off;
                if int == 1,
                    title('Energy over time');
                elseif int == 2,
                    l = legend('Max', 'Circular mean', 'Location', 'SouthOutside');
                    legend boxoff;
                    s = subplot(339); spos = get(gca, 'Position');
                    set(l, 'Position', spos); axis off;
                    xlabel('Seconds');
                end
                
            end
            saveas(gcf, sprintf('~/Data/MotionEnergy/motionenergy_P%02d-S%d_trial%d.eps', sj, setup.session, trial), 'epsc');
        end % trial
        
    end % file
end % sj


