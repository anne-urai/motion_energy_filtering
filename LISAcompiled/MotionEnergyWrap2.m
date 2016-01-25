function [] =  MotionEnergyWrap2(sj, session, block)
% for one specific subject and session, run the motion energy filtering
% after emails with John (June), make sure this compiled script doesnt use multithreading

begin = tic;

% input arguments in the command line will be treated as text
if ischar(sj),          sj      = str2double(sj); end
if ischar(session),     session = str2double(session); end
if ischar(block),       block 	= str2double(block); end

% ONLY RUN IF THE OUTPUT FILE ISNT THERE YET!
if exist(sprintf('~/Data/MotionEnergy/motionenergy_P%02d_s%d_b%d.mat', sj, session, block), 'file'),
    fprintf('file already exists, skipping: motionenergy_P%02d_s%d_b%d.mat', sj, session, block);
    return % end
end


files = dir(sprintf('~/Data/Commitment/Behav/P%d_s%d_20*.mat', sj, session));
assert(length(files)==1);

% get tmpdir scratch space (LISA)
tmppath = getenv('TMPDIR');
loadedFile = 0;
while ~loadedFile,
    try
        randi(60); % wait a bit
        % load from tmpdir if already exists
        load(sprintf('%s/%s', tmppath, files.name));
        loadedFile = 1;
    catch
        % copy file to scratch
        copyfile(sprintf('~/Data/Commitment/Behav/%s', files.name), ...
            sprintf('%s/%s', tmppath, files.name));
        fprintf('copying file to %s/%s \n', tmppath, files.name);
        disp('skipping file copy, error, file already exists')
    end
end

% cheat a bit and skip scratch space for now
load(sprintf('~/Data/Commitment/Behav/%s', files.name));

% to avoid window UI bugs
display.dist        = window.dist;
display.res         = window.res;
display.width       = window.width;
display.frameRate   = window.frameRate;
display.center      = window.center;

clear sound audio

% run this once to find optimal fft algorithm
fftw('planner', 'exhaustive');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute some general things for this file
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg.ppd             = deg2pix(display, 1);  % pixel per degree
cfg.srange          = [-0.7 0.7];           % predetermine size the filter  space will have
cfg.frameRate       = display.frameRate;

% to avoid bugs with frameRates that are slightly below 60 Hz
if cfg.frameRate < 60,
    cfg.frameRate = 60.1;
end

cfg.trange          = [0 0.2];              % temporal range
cfg.filtsize        = [ceil(diff(cfg.srange)*cfg.ppd) ...
    ceil(diff(cfg.srange)*cfg.ppd) ...
    ceil(diff(cfg.trange)*cfg.frameRate)];

% pad with some pixels for valid convolution
cfg.stimpad         = 50;
cfg.stimsize        = [2*dots.radius+1+cfg.stimpad 2*dots.radius+1+cfg.stimpad setup.nframes]; %  predetermine the size the cloud of dots will have
cfg.n_convolution   = cfg.stimsize + cfg.filtsize - 1;  % 'full'  size of the convolution
cfg.nextpow2        = 2.^nextpow2(cfg.n_convolution);
cfg.validsize       = cfg.stimsize - cfg.filtsize + 1;  % 'valid' size of convolution

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate filters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta               = 1:360; % around the unit circle

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. CREATE SPATIAL AND TEMPORAL FILTERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% all parameters from Kiani et al., http://www.jneurosci.org/content/28/12/3017.full

% time and space axes
x = cfg.srange(1) : 1/cfg.ppd : cfg.srange(2);  % range in space, in degree, spaced by the size in degree of one pixel
y = cfg.srange(1) : 1/cfg.ppd : cfg.srange(2);  % range in space, in degree, spaced by the size in degree of one pixel
t = cfg.trange(1) : 1/cfg.frameRate: cfg.trange(2);   % range in time, in s, spaced by framerate

% important: use an uneven nr of points, easier for the later convolution
assert(mod(numel(x), 2) == 1, 'x should have an uneven nr of samples');
assert(mod(numel(y), 2) == 1, 'y should have an uneven nr of samples');
assert(mod(numel(t), 2) == 1, 't should have an uneven nr of samples');

% check if the filter has the size I thought it would get
assert(all([length(x) length(y) length(t)] == cfg.filtsize), 'filter does not have prespecified size');

% constants of the spatial filters
sc      = 0.35;
sg      = 0.05;

% create spatial mesh
[xmesh, ymesh] = meshgrid(x,y);

% constants of the temporal filters
k       = 60;
nslow   = 3;
nfast   = 5;

% TEMPORAL FUNCTIONS
% put them in the third dimension for the outer product
g1(1,1,:) = (k*t).^nslow .* exp(-k*t) .* (1/factorial(nslow) - ((k*t).^2) / factorial(nslow + 2));
g2(1,1,:) = (k*t).^nfast .* exp(-k*t) .* (1/factorial(nfast) - ((k*t).^2) / factorial(nfast + 2));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make a separate filter for each direction
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for thistheta = theta,
    
    % rotate the meshgrid
    yprime  = - xmesh .* sind(thistheta) + ymesh .* cosd(thistheta);
    xprime  =   xmesh .* cosd(thistheta) + ymesh .* sind(thistheta);
    
    % SPATIAL FUNCTIONS two fourth order cauchy functions
    % transpose to match my unit circle directionality
    alpha      = atand(xprime ./ sc);
    f1rot      = (cosd(alpha).^4 .* cosd(4*alpha) .* exp(-yprime.^2/(2*sg^2)))';
    f2rot      = (cosd(alpha).^4 .* sind(4*alpha) .* exp(-yprime.^2/(2*sg^2)))';
    
    % these two linear filters are in space-time quadrature
    filt1   = bsxfun(@times, f1rot, g1) + bsxfun(@times, f2rot, g2);
    filt2   = bsxfun(@times, f2rot, g1) - bsxfun(@times, f1rot, g2);
    
    % save the filters to output
    % already run the filter fft!
    filters(find(thistheta==theta)).one = single(filt1);
    filters(find(thistheta==theta)).two = single(filt2);
    
    if 0,
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot info about the filters
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        figure; colormap bone;
        subplot(441); imagesc(x, y, f1rot); xlabel('x'); ylabel('y'); title('f1');
        subplot(445); imagesc(x, y, f2rot); xlabel('x'); ylabel('y'); title('f2');
        
        subplot(4,4,[2 6]); plot(t, squeeze(g1), t, squeeze(g2)); title('g1 g2'); xlabel('time');
        
        % freq spectrum of filters
        subplot(4,4,3); imagesc(x, t, squeeze(sum((filter.one),2)) .* squeeze(sum((filter.two),2)));
        %xlabel('time'); ylabel('x');
        title('freq along xdir');
        subplot(4,4,7); imagesc(y, t, squeeze(sum((filter.one),1)) .* squeeze(sum((filter.two),1)));
        %xlabel('time'); ylabel('y');
        title('freq along ydir');
        
        % collapsed over x dir
        subplot(4,4,9); imagesc(y, t, squeeze(sum(filt1,1)));
        title('filt1'); xlabel('y'); ylabel('t');
        subplot(4,4,13); imagesc(y, t, squeeze(sum(filt2,1)));
        title('filt2'); xlabel('y'); ylabel('t');
        
        % collapse over y dir
        subplot(4,4,10); imagesc(x, t, squeeze(sum(filt1, 2)));
        title('filt1'); xlabel('x'); ylabel('t');
        subplot(4,4, 14); imagesc(x, t, squeeze(sum(filt2,2)));
        title('filt2'); xlabel('x'); ylabel('t');
        
        % collapse over t dir
        subplot(4,4,11); imagesc(x, y, squeeze(sum(filt1, 3)));
        title('filt1'); xlabel('x'); ylabel('y');
        subplot(4,4, 15); imagesc(x, y, squeeze(sum(filt2,3)));
        title('filt2'); xlabel('x'); ylabel('y');
        
        suplabel([ 'Rotation ' num2str(thistheta) ' degrees'], 't');
        saveas(gcf, sprintf('~/Data/MotionEnergy/filters/filterimg_%d.eps', thistheta), 'epsc');
        close all;
    end
end

fprintf('preparing filters took %.2f seconds \n', toc(begin));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load stimuli and run convolution operation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preallocate the motion energy output
motionenergy.one = single(nan(setup.ntrials, length(theta), cfg.validsize(3)));
motionenergy.two = single(nan(setup.ntrials, length(theta), cfg.validsize(3)));

% go over both intervals
intervals = char('one', 'two');

for trial = 1:setup.ntrials,
    
    % find the block and trial we'll work on
    thistrial  = sub2ind(size(setup.direction.meanestimate'), trial, block);
    
    if setup.choice(block, trial) > -1, % no point doing this for the 1/3rd of choice only trials
        
        for int = 1:2,
            
            tic;
            
            % 1. CHANGE THE stim COORDINATES INTO A REPRESENTATION IN X Y Z SPACE
            stimulus = squeeze(stim.(intervals(int,:))(thistrial,:, :, :));
            
            % take this from the setup struct
            xres    = display.res.width;
            yres    = display.res.height;
            nfr     = setup.nframes;
            
            % preallocate
            stimrep = zeros(xres, yres, nfr);
            
            for f = 1:nfr,
                
                posx = squeeze(stimulus(f, 1, :));
                posy = squeeze(stimulus(f, 2, :));
                
                % normalize by the range of values, so that all position estimates are
                % round to the nearest pixel to create matrix
                posx = round(posx + display.center(1));
                posy = round(posy + display.center(2));
                
                % put those coordinates in the stim representation matrix
                for i = 1:length(posx), % for each dot
                    % Y X T
                    stimrep(posx(i), posy(i), f) = 1; % put a pixel in the matrix
                end
            end
            
            % continue with only those parts of the stimulus in
            % the x and y dir that contain the cloud of dots
            x2use = [display.center(1)-dots.radius-cfg.stimpad/2 : display.center(1)+dots.radius+cfg.stimpad/2];
            y2use = [display.center(2)-dots.radius-cfg.stimpad/2 : display.center(2)+dots.radius+cfg.stimpad/2];
            
            stimrep = stimrep(x2use, y2use, :);
            assert(all(size(stimrep) == cfg.stimsize), 'stimulus size not correctly computed');
            
            % get the fft of the stimulus
            stim_fft = fftn(stimrep, cfg.n_convolution);
            
            clear stimulus stimrep
            % determine the size the output will have
            for s = 1:3,
                size2use(s,:) = [cfg.stimsize(s)-cfg.validsize(s)+1  ...
                    cfg.n_convolution(s) - (cfg.stimsize(s)-cfg.validsize(s)) ];
            end
            
            fprintf('P%d_s%d_b%d_trial%d_int%d took %.3f seconds to prepare \n', sj, session, block, trial, int, toc);
            
            tic;
            
            %% 2. FILTER AT EACH THETA
            for thistheta = theta,
                
                % run multiplication in the frequency domain
                % this is where the bulk of the computation happens
                
                resp1       = ifftn(stim_fft .* fftn(filters(find(thistheta==theta)).one, cfg.n_convolution), cfg.n_convolution);
                resp2       = ifftn(stim_fft .* fftn(filters(find(thistheta==theta)).two, cfg.n_convolution), cfg.n_convolution);
                
                % use only valid part of the result, slightly smaller than the size of the input
                resp1       = resp1(size2use(1,1):size2use(1,2), size2use(2,1):size2use(2,2), size2use(3,1):size2use(3,2));
                resp2       = resp2(size2use(1,1):size2use(1,2), size2use(2,1):size2use(2,2), size2use(3,1):size2use(3,2));
                
                % sum and square the results of the two filters in quadrature, see Adelson & Bergen
                energy      = (resp1.^2 + resp2.^2);
                
                % collapse over the x and y directions, just give the time output
                % take the square root, to normalize the responses
                motionenergy.(intervals(int,:))(thistrial, (thistheta == theta), :) = single(sqrt(squeeze(sum(sum(energy)))))';
                
            end % theta
            
            fprintf('P%d_s%d_b%d_trial%d_int%d took %.3f seconds to filter \n', sj, session, block, trial, int, toc);
            
            % also save the real direction according to setup for easier plotting
            switch int
                case 1
                    realdir = setup.direction.stimone(block, trial);
                case 2
                    realdir = setup.direction.stimtwo(block, trial);
            end
            
            realdirection.(intervals(int,:))(thistrial) = realdir;
            
        end % int
    else
        fprintf('P%d_s%d_b%d_trial%d choice-1 \n', sj, session, block, trial);
        for int = 1:2, realdirection.(intervals(int,:))(thistrial) = NaN; end
    end
    
end % thistrial

% save to scratch space
tmppath = getenv('TMPDIR');
filename = sprintf('motionenergy_P%02d_s%d_b%d.mat', sj, session, block);
save([tmppath '/' filename], '-mat', ...
    'realdirection', 'motionenergy', 'cfg', 'setup', 'display', 'results', 'filters', 'dots');

% copy to home dir
copyfile([tmppath '/' filename], ['~/Data/MotionEnergy/' filename]);

fprintf('SAVED /home/aeurai/Data/MotionEnergy/motionenergy_P%02d_s%d_b%d.mat', sj, session, block)
fprintf('P%d_s%d_b%d took %.3f seconds (= %.3f hours) to fully process \n', sj, session, block, toc(begin), toc(begin)/3600);

end % function end
