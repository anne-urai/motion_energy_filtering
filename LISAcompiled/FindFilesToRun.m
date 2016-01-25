params2run = [];

for sj = 0:15,
	for session = 1:11,
	for block = 1:5,
            % ONLY RUN IF THE OUTPUT FILE ISNT THERE YET!
            if ~exist(sprintf('~/Data/MotionEnergy/motionenergy_P%02d_s%d_b%d.mat', sj, session, block), 'file'),
                
                % check if the input is right
                files = dir(sprintf('~/Data/Commitment/Behav/P%d_s%d_20*.mat', sj, session));
                if length(files) == 1,
                    params2run = [params2run; sj session block];
                elseif length(files) == 0,
                    fprintf('~/Data/Commitment/Behav/P%d_s%d_20*.mat \n', sj, session)
                else
                    error('multiple files found!');
                end
                
            end
        end
    end
end
% write to disk
dlmwrite('lastparams',params2run,' ')
