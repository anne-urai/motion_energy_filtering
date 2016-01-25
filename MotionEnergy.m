function [motionenergy, filters] =  MotionEnergy(stimulus, display, theta, plotme, filters)
% implements motion energy filtering as described in Adelson & Bergen
% (1985), Kiani et al. (2008), Bollimunta et al. (2012)
% 
% see also these excellent tutorials on the 2D equivalent of these filters
% http://www.georgemather.com/Model.html (George Mather)
% http://mplab.ucsd.edu/~marni/CSHL_Tutorials/MotionEnergy.tar 
%   (Simoncelli, Glimcher, Chichilnisky)
% 
% INPUT
% stimulus, a 3d matrix containing the x-y-t stimulus movie
%  display, structure containing the following fields:
% 	display.frameRate (in Hz)
%   display.ppd (pixels/degree) OR all of the following:
%       display.width (in cm of physical screen)
%       display.res.width (in pixels)
%       display.dist (between eyes and screen, in cm)
% theta, the angles at which motion energy should be computed
%   (rightwards = 0, counterclockwise).
% plotme, generate graphical output (default = false);
% 
% OUTPUT
% motionenergy, length(theta) x size(stimulus, 3) motionenergy filter output
% filters (optional)
% 
% Anne Urai, 2016
% anne.urai@gmail.com / anneurai.net

if ~exist('plotme', 'var'); plotme = 1; end
fftw('planner', 'exhaustive'); % run this once to find optimal fft algorithm

% ======================================================= %
% compute some general things for this file
% ======================================================= %

if ~isfield(display, 'ppd'),
    pixSize = display.width/display.res.width;   %cm/pix
    sz      = 2*display.dist*tan(pi*ang/(2*180));  %cm
    cfg.ppd = round(sz/pixSize);   %pix
else
    cfg.ppd = display.ppd;
end

% spatial range of the filter
cfg.srange          = [-0.7 0.7];    

% temporal range of the filter
cfg.frameRate       = display.frameRate;
% to avoid bugs with frameRates that are slightly below 60 Hz
if cfg.frameRate < 60,
    cfg.frameRate = 60.1;
end
cfg.trange          = [0 0.2];    

% filter size
cfg.filtsize        = [ceil(diff(cfg.srange)*cfg.ppd) ...
    ceil(diff(cfg.srange)*cfg.ppd) ...
    ceil(diff(cfg.trange)*cfg.frameRate)];

% ======================================================= %
% 1. CREATE SPATIAL AND TEMPORAL FILTERS
% ======================================================= %

% all parameters from Kiani et al., http://www.jneurosci.org/content/28/12/3017.full

% time and space axes
x = linspace(cfg.srange(1), cfg.srange(2), cfg.filtsize(1));
y = linspace(cfg.srange(1), cfg.srange(2), cfg.filtsize(2));
t = linspace(cfg.trange(1), cfg.trange(2), cfg.filtsize(3));

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

% ======================================================= %
% make a separate filter for each direction
% ======================================================= %

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
    filters(['theta' num2str(find(thistheta==theta))]).one = single(filt1);
    filters(['theta' num2str(find(thistheta==theta))]).two = single(filt2);
    
    if plotme,
        
        % ======================================================= %
        % plot info about the filters
        % ======================================================= %
        
        figure; colormap bone;
        subplot(441); imagesc(x, y, f1rot); xlabel('x'); ylabel('y'); title('f1');
        subplot(445); imagesc(x, y, f2rot); xlabel('x'); ylabel('y'); title('f2');
        
        subplot(4,4,[2 6]); plot(t, squeeze(g1), t, squeeze(g2)); title('g1 g2'); xlabel('time');
        
        % freq spectrum of filters
        subplot(4,4,3); imagesc(x, t, squeeze(sum((filter.one),2)) .* squeeze(sum((filter.two),2)));
        xlabel('time'); ylabel('x');
        title('freq along xdir');
        subplot(4,4,7); imagesc(y, t, squeeze(sum((filter.one),1)) .* squeeze(sum((filter.two),1)));
        xlabel('time'); ylabel('y');
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
    end
end
fprintf('preparing filters took %.2f seconds \n', toc(begin));

% ======================================================= %
% load stimuli and run convolution operation
% ======================================================= %

% preallocate the motion energy output
motionenergy = single(nan(length(theta), cfg.validsize(3)));

% FILTER AT EACH THETA
for thistheta = theta,
    
    % run multiplication in the frequency domain rather than convolution in
    % the time domain - this is where the bulk of the computation happens
    resp1       = convnfft_light(stimulus, filters.(['theta' num2str(find(thistheta==theta))]).one);
    resp2       = convnfft_light(stimulus, filters.(['theta' num2str(find(thistheta==theta))]).two);

    % sum and square the results of the two filters in quadrature, see Adelson & Bergen
    energy      = (resp1.^2 + resp2.^2);
    
    % collapse over the x and y directions, just give the time output
    % take the square root, to normalize the responses
    motionenergy(find(thistheta == theta), :) = single(sqrt(squeeze(sum(sum(energy)))))';
    
end % theta

end % function end


function A = convnfft_light(A,B)
% light version of http://nl.mathworks.com/matlabcentral/fileexchange/24504-fft-based-convolution
% Anne Urai, 2016
% anne.urai@gmail.com

% AEU: eliminate loop
m = size(A);
n = size(B);

% IFUN function will be used later to truncate the result
% M and N are respectively the length of A and B in some dimension
ifun = @(m,n) (1:m) + ceil((n-1)/2);

% subset of datapoints for same size
subs = arrayfun(ifun, size(A), size(B), 'uniformoutput', 0);

% compute the FFT length
l = m+n-1;

% fftn is faster than looping
A = fftn(A, l);
B = fftn(B, l);

% multiply in frequency domain
A = A.*B;
clear B

% back to time domain
A = ifftn(A);

% Truncate the results
A = A(subs{:});

end % convnfft

