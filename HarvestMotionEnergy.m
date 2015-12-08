% harvest motion energy and write to a nice file that we want
% Anne Urai, 11 June 2015

clear; close all; clc;
dbstop if error;
subjects = 0:15;
% subjects = 0;
%subjects = 1;

for sj = subjects,
    clearvars -except sj subjects;
    
    % load the datafile that we already have
    load(sprintf('~/Data/Commitment/Pupil/P%02d_data.mat', sj));
    
    % correct for weirdness
    if sj == 5,
        dat(208:276, 11) = 4;
        dat(277:345, 11) = 5;
    elseif sj == 8,
        dat(484:552, 11) = 3;
    end
    
    % start a new file with motionenergy
    allmotion.one   = nan(length(dat), 33);
    allmotion.two   = nan(length(dat), 33);
    
    theta       = 1:360; % for circular mean
    
    % nr of sessions can differ, 11 for most sj
    sessions = unique(dat(:, 12));
    for session = sessions',
        for block = 1:5,
            
            % correct for weirdness
            if sj == 5 && session == 1 && block > 3,
                iblock = block - 3; % block 1 and 2 were run again, should be 4 and 5
            else
                iblock = block;
            end
            
            % get the file
            load(sprintf('~/Data/MotionEnergy/motionenergy_P%02d_s%d_b%d.mat', sj, session, block));
            fprintf('~/Data/MotionEnergy/motionenergy_P%02d_s%d_b%d.mat \n', sj, session, block);
            
            for trial = 1:setup.ntrials,
                
                try
                    % match the dat file with the motionenergy file
                    datTrl = find(dat(:, 12) == session & dat(:, 11) == block & dat(:, 10) == trial);
                    if isempty(datTrl), disp('skipping trial'); continue; end % if we removed this trial from the dat file (perhaps no pupil data?) skip here
                    assert(length(datTrl) == 1, 'something went horribly wrong'); % this should match only 1 trial!
                    
                    % check the match
                    assert(isequal(setup.direction.stimone(iblock, trial) - setup.refdirec(iblock, trial), ...
                        dat(datTrl, 1)), 'no match for direction int1');
                    assert(isequaln(results.binary.response(iblock, trial),dat(datTrl, 3)), 'no match for binary response');
                    assert(isequaln(results.binary.correct(iblock, trial),dat(datTrl, 4)), 'no match for binary correct');
                catch
                    disp('skipping trial, no match found');
                    continue;
                end
                
                % check if this was a trial with an evaluation
                if setup.choice(iblock, trial) == -1 || dat(datTrl, 2) == -1,
                    if ~(setup.choice(iblock, trial) == -1 && dat(datTrl, 2) == -1), disp('nochoice info not consistent'); end
                    continue;
                end
                
                % do further checks on the second interval
                assert(isequal(setup.direction.stimtwo(iblock, trial) - setup.refdirec(iblock, trial), ...
                    dat(datTrl, 6)), 'no match for direction int2');
                if ~isnan(results.estimation.RT(iblock, trial)) && isnan(dat(datTrl, 9)),
                    assert(ismembertol(results.estimation.RT(iblock, trial), dat(datTrl, 9), 0.1), 'no match for RT');
                end
                
                for int = 1:2,
                    
                    % loop over the two fields
                    switch int
                        case 1
                            intfld = 'one';
                        case 2
                            intfld = 'two';
                    end
                    
                    % extract what we need
                    thisenergy  = squeeze(motionenergy.(intfld)(trial + (iblock-1)*69, :, :))';
                    
                    % compute the circular mean, represents the MT population
                    % vector response
                    circmean = rad2deg(circ_mean(repmat(deg2rad(theta), size(thisenergy, 1), 1),...
                        thisenergy, 2));
                    circmean = circmean - setup.refdirec(iblock, trial);
                    % correct all directions with the reference
                    
                    % correct for unit circle discontinuities
                    if any(circmean < -100),
                        circmean(circmean<-100)    = circmean(circmean<-100) + 360;
                    end
                    
                    % add to allmotion matrix
                    allmotion.(intfld)(datTrl, 1:length(circmean)) = circmean;
                end
            end % trial
        end % block
    end % session

    % save all the data we need
    save(sprintf('~/Data/MotionEnergy/Extracted/P%02d_motiondata.mat', sj), 'dat', 'allmotion');
    
    scatterPlot = false;
    if 1,
        
        % make an overview plot of the fluctuations in evidence
        direcs = [-20 -10 0 10 20];
        cols = linspecer(length(direcs), 'sequential');
        clf;
        
        if scatterPlot,
            for int = 1:2,
                
                switch int
                    case 1
                        intcol = 1; intfld = 'one';
                    case 2
                        intcol = 6; intfld = 'two';
                end
                
                s1 = subplot(3,3,int);
                hold on;
                sp = scatter(repmat(1:size(allmotion.(intfld), 2), 1, size(allmotion.(intfld), 1)), ...
                    reshape(allmotion.(intfld), [1, numel(allmotion.(intfld))]),...
                    1, reshape(repmat(dat(:, intcol), 1, size(allmotion.(intfld), 2)), [1, numel(allmotion.(intfld))]), 'filled');
                alpha(sp, 0.5);
                axis tight;
                ylims = get(gca, 'ylim'); ylim([-max(abs(ylims)) max(abs(ylims))]);
                xlim([0 size(allmotion.(intfld), 2)+1]);
                set(gca, 'box', 'off', 'tickdir', 'out');
                title(['Interval ' intfld]);
                if int == 1, ylabel('Filtered direction'); end
                xlabel('Time');
            end
            
        else
            
            for int = 1:2,
                
                switch int
                    case 1
                        intcol = 1; intfld = 'one';
                    case 2
                        intcol = 6; intfld = 'two';
                end
                
                s1 = subplot(3,3,int);
                hold on;
                
                for d = 1:length(direcs),
                    trls    = find(dat(:,intcol) == direcs(d));
                    clear ptmp;
                    ptmp = plot(allmotion.(intfld)(trls, :)', '-', 'color', cols(d, :), 'LineWidth', 0.5);
                    thishand(d) = ptmp(1);
                    
                    % make all the lines transparant
                    for i = 1:length(ptmp),   ptmp(i).Color(4) = 0.05;    end
                end
                axis tight; ylims = get(gca, 'ylim'); ylim([-max(abs(ylims)) max(abs(ylims))]);
                xlim([0 size(allmotion.(intfld), 2)+1]);
                set(gca, 'box', 'off', 'tickdir', 'out');
                title(['Interval ' intfld]);
                if int == 1, ylabel('Filtered direction'); end
                xlabel('Time');
            end
            
            lh = legend(thishand', {'20', '10', '0', '-10', '-20'});
            s3 = subplot(333); spos = get(s3, 'Position');
            set(lh, 'Position', spos, 'box', 'off');
            axis off;
        end
        
        try
            suplabel(sprintf('P%02d - threshold %.1f %% coherence', sj, setup.coherencelevel*100), 'x');
        catch
            suplabel(sprintf('P%02d - threshold %.1f %% coherence', sj, setup.cohlevel*100), 'x');
        end
        saveas(gcf, sprintf('~/Data/MotionEnergy/Extracted/P%02d_energyPlot.pdf', sj));
    end
    
end
