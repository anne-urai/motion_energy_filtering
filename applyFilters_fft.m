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

function energyright_fft = applyFilters_fft(stim, f1, f2, g1, g2)
% apply the filters using 3d convolution
% !!! IMPORTANT !!!
% THIS RETURNS ENERGY ONLY IN THE RIGHTWARD DIRECTION, NOT OPPONENT ENERGY
% AS COMPUTED IN applyFilters.m !

% these two linear filters are in space-time quadrature
g1a(1,1,:)      = g1;
g2a(1,1,:)      = g2;
filt1           = bsxfun(@times, f1, g1a) + bsxfun(@times, f2, g2a); % A1 + B2
filt2           = bsxfun(@times, f1, g2a) - bsxfun(@times, f2, g1a); % A2 - B1

% take fft of the filter itself
filters.one     = fftn(single(filt1), n_convolution);
filters.two     = fftn(single(filt2), n_convolution);

% take fft of the stimulus
stim_fft        = fftn(stim, n_convolution);

% run multiplication in the frequency domain
resp1           = ifftn(stim_fft .* filters.one, n_convolution);
resp2           = ifftn(stim_fft .* filters.two, n_convolution);

if 0,
    
    stimsize        = size(stim);
    n_convolution   = size(stim) + size(filt1) - 1;
    validsize       = size(stim) - size(filt1) + 1;
    filtsize        = size(filt1);
    
    % determine the size the output will have
    size2use = nan(3,2);
    for s = 1:3,
        size2use(s,:) = [(n_convolution(s) - stimsize(s))/2+1  ...
            n_convolution(s) - (n_convolution(s) - stimsize(s))/2];
    end
    
    % use only valid part of the result, slightly smaller than the size of the input
    resp1       = resp1(size2use(1,1):size2use(1,2), size2use(2,1):size2use(2,2), size2use(3,1):size2use(3,2));
    resp2       = resp2(size2use(1,1):size2use(1,2), size2use(2,1):size2use(2,2), size2use(3,1):size2use(3,2));
end

% sum and square the results of the two filters in quadrature, see Adelson & Bergen
energyright_fft      = (resp1.^2 + resp2.^2);



if 0,
    % ===================================================== %
    % 3D CONVOLUTION
    % ===================================================== %
    
    % these two linear filters are in space-time quadrature
    g1a(1,1,:)      = g1;
    g2a(1,1,:)      = g2;
    filt1           = bsxfun(@times, f1, g1a) + bsxfun(@times, f2, g2a); % A1 + B2
    filt2           = bsxfun(@times, f1, g2a) - bsxfun(@times, f2, g1a); % A2 - B1
    
    stimsize        = size(stim);
    n_convolution   = size(stim) + size(filt1) - 1;
    validsize       = size(stim) - size(filt1) + 1;
    filtsize        = size(filt1);
    
    % take fft of the filter itself
    filters.one     = fftn(single(filt1), n_convolution);
    filters.two     = fftn(single(filt2), n_convolution);
    
    % take fft of the stimulus
    stim_fft        = fftn(stim, n_convolution);
    
    % run multiplication in the frequency domain
    resp1           = ifftn(stim_fft .* filters.one, n_convolution);
    resp2           = ifftn(stim_fft .* filters.two, n_convolution);
    
    if 0,
        % determine the size the output will have
        size2use = nan(3,2);
        for s = 1:3,
            size2use(s,:) = [(n_convolution(s) - stimsize(s))/2+1  ...
                n_convolution(s) - (n_convolution(s) - stimsize(s))/2];
        end
        
        % use only valid part of the result, slightly smaller than the size of the input
        resp1       = resp1(size2use(1,1):size2use(1,2), size2use(2,1):size2use(2,2), size2use(3,1):size2use(3,2));
        resp2       = resp2(size2use(1,1):size2use(1,2), size2use(2,1):size2use(2,2), size2use(3,1):size2use(3,2));
    end
    
    % sum and square the results of the two filters in quadrature, see Adelson & Bergen
    energyright_fft      = (resp1.^2 + resp2.^2);
    
    clf;
    % this should be the same as energy_right!
    subplot(223); plot(energyright_fft(:), energy_right(:), '.'); xlabel('fft'); ylabel('filter');
    subplot(222); plot(squeeze(sum(sum(energyright_fft)))); hold on;  title('FFT'); ylim([1540 1640]);
    subplot(221); plot(squeeze(sum(sum(energy_right)))); xlabel('Time'); title('filter right');  ylim([1540 1640]);
    
    % FFT-based method returns the same time-course,
    figure; colormap(parula);
    subplot(221); imagesc(energyright_fft(:, :, 20)); axis image;
    set(gca,'xtick', linspace(50, 850,  17), 'ytick', linspace(50, 850, 17));
    set(gca,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', 'xcolor', 'k', 'ycolor', 'k');
    subplot(222); imagesc(energy_right(:, :, 20)); axis image;
    set(gca,'xtick', linspace(50, 850,  17), 'ytick', linspace(50, 850,  17));
    set(gca,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', 'xcolor', 'k', 'ycolor', 'k');
    
    
    figure; colormap(parula);
    subplot(221); plot(energyright_fft(:, 400, 20));
    subplot(222); plot(energy_right(:, 400, 20));
    
    for s = -100:1:100,
        
        clf; plot(energy_right(:, 400, 20));
        hold on; plot(energyright_fft(:, 400+s, 20));
        title(sprintf('lag = %d', s));
        waitforbuttonpress;
    end
    
end

end