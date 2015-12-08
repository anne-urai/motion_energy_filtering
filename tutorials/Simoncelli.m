%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This tutorial was created for use in the short course
%%% "Computational Neuroscience: Vision", held at Cold Spring Harbor
%%% Laboratories in Summer 1996.  It is provided "as is", without
%%% express or implied warranty, for educational purposes only, and
%%% is not to be used, rewritten, or adapted as the basis of a
%%% commercial software or hardware product.  Please do not
%%% distribute this file without prior permission of the instructors.
%%% Many hours were spent developing these tutorials and we ask your
%%% compliance to assure proper attribution and maintainance.  If you
%%% would like to use this tutorial in a course, please contact one
%%% of the following instructors:
%%% 
%%% Eero Simoncelli <eero.simoncelli@nyu.edu>
%%% Paul Glimcher <glimcher@cns.nyu.edu>
%%% EJ Chichilnisky <ej@salk.edu>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Adelson-Bergen Motion Model Tutorial
%
%This tutorial presents a straightforward implementation of the 
%Spatio-temporal energy model described by Adelson & Bergen
%in 'Spatiotemporal energy models for the perception of motion'
%(JOSA-A,2:284-299).  
%
%The paper, and this tutorial, shows how linear separable filters
%can be added and subtracted to form space-time oriented linear 
%filters.  These oriented filters, which are selective to
%direction of motion, can be combined nonlinearly to produce 
%filters that are not only motion selective, but are also insensitive
%to phase.  Such response properties are typical of V1 complex cells.
%
%It might be useful to have a copy of the article with you as you work
%through this tutorial, since some of the figures generated here are
%much like those in the paper.

%Written by Geoff Boynton rudely during lectures at Cold Spring Harbor
%July, 1996.

% Run by cutting and pasting from the tutorial file into the Matlab 
% Command window.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%First, we'll create some linear filters that are separable in space
%and time: 

nx=20;  %number of spatial samples in the filter
nt=20;  %number of temporal samples

%Define the space axis
max_x =1;     %Half-width of Gabor (pixels)
x_filt=linspace(-max_x,max_x,nx);

%Define the time axis
max_t=1.0; %Duration of impulse response (sec)
t_filt=linspace(0,max_t,nt)';

%Spatial functions:
%Create a pair of Gabor functions; one even and one odd.
sx=0.5;   %standard deviation of Gaussian
sfx=1.1;  %spatial freqency of carrier

gaussx=exp(-x_filt.^2/sx.^2);				%Gaussian envelope
even_x=cos(2*pi*sfx*x_filt).*gaussx;			%Even Gabor
odd_x=sin(2*pi*sfx*x_filt).*gaussx;			%Odd Gabor

%Plot the spatial profiles (See fig 6a)
figure(1)
subplot(1,2,1);
plot(x_filt,even_x,x_filt,odd_x)
xlabel('Space'); 

%Temporal functions:
%Create a pair of temporal functions; one fast and one slow.
k=25;  %Temporal scale factor 
slow_t=temp_imp_resp(5,k,t_filt);
fast_t=temp_imp_resp(3,k,t_filt);

%Plot the temporal profiles (See fig 6b)
subplot(1,2,2)
plot(t_filt,slow_t,t_filt,fast_t);
xlabel('time');

%Create separable linear filters by taking the outer
%products of the four combinations of space and time 
%impulse responses.  See figure 
even_slow= slow_t * even_x;
even_fast= fast_t * even_x ;
odd_slow = slow_t *odd_x ;
odd_fast = fast_t * odd_x ;

%Show the four impulse response functions as images
%See fig 6c (top row)
figure(2)
subplot(2,4,1);show_image(odd_fast,[50,50]);
subplot(2,4,2);show_image(odd_slow,[50,50]);
subplot(2,4,3);show_image(even_fast,[50,50]);
subplot(2,4,4);show_image(even_slow,[50,50]);

%These filters can be added or subtracted in pairs to produce
%new linear filters that are space-time oriented.  For example,
%adding the odd_fast filter to the even_slow filter results in
%a leftward-selective linear filter.  
 
leftward_1=odd_fast+even_slow; 
leftward_2=-odd_slow+even_fast;
rightward_1=-odd_fast+even_slow;
rightward_2=odd_slow+even_fast;
 
%Here are images of the four space-time oriented filters shown
%in A&B's figure 4c.

subplot(2,4,5);show_image(leftward_1,[50,50]);
subplot(2,4,6);show_image(leftward_2,[50,50]);
subplot(2,4,7);show_image(rightward_1,[50,50]);
subplot(2,4,8);show_image(rightward_2,[50,50]);


%Let's measure the response of these filters to various stimuli.
%First, we'll generate the moving bar like the stimulus shown in 
%figure 15a.

stim_dur=6; 	   %Duration of stimulus (seconds)
stim_width=4;      %Width of stimulus (pixels)

%x-axis for stimulus
x_stim=-stim_width:x_filt(2)-x_filt(1):stim_width;
%t-axis for stimulus
t_stim=(0:t_filt(2):stim_dur)';

%This part generates the 2-d (space-time) stimulus, 'stim'.  
t_freq=0.25;  %Temporal frequency of oscillation (Hz)
amplitude= 3;	%Amplitude of oscillation (pixels)
x_profile=x_stim<0;
t_profile=-cos(t_freq*2*pi*t_stim)*amplitude;
stim=(ones(length(t_profile),1) * x_stim) < (t_profile*ones(1,length(x_profile)));
stim = double(stim);  %%% Marni added to update.

%Show the stimulus as a space-time image
figure(3)
subplot(1,2,1)
show_image(stim);

%To compute a filter's response to such a stimulus, we simply perform
%a 2-d convolution of the stimulus with the filter.

%We'll use the first rightward sensitive filter, 'rightward_1'

resp=conv2(stim,rightward_1,'same');

%Plot the result
subplot(1,2,2)
show_image(resp);

%The way to interpret the response image is to consider each column
%of the image as the the time-course of the response of a filter
%centered at that column.  The filter responds when the edge passes
%past the column.  For example, if we consider the response of a filter
%centered on the center column of the response image:

%column=size(stim,2)/2; %center column
column=round(size(stim,2)/2); %center column   (Marni's update)

line([column,column],[0,length(t_stim)],'Color','r');

%We can look at the time course of the response by pulling out the
%center column and plotting:

figure(4)
%plot(t_stim,resp(:,column));
plot(t_stim(1:115),resp(:,column));  %% Marni's revision. 
xlabel('Time (sec)');
ylabel('Response');

%Notice that the intensity of response is greater when the bar is
%moving rightward then for leftward movement.  But the response is
%not zero for leftward movment.  This is because the oriented filters
%are not perfectly oriented.  You can see in the image of the
%oriented filters that the light and dark regions curve downward.

%This can be corrected by using a true oriented Gabor function as
%an oriented filter.  This is not biologically plausible, since such
%a function is not causal in time.

%It's also interesting to study the response of these filters to drifting
%sinewave gratings.  First, we'll generate the stimulus:

%Start with sinewave grating
stim_dur=12;    %seconds
stim_width=4;   %multiples of filter
velocity=60;  %pixels/sec

x_stim=-stim_width:x_filt(2)-x_filt(1):stim_width;
x_profile=sin(2*pi*x_stim/2);
t_stim=(0:t_filt(2):stim_dur)';
t_profile=(t_stim(1:length(t_stim)/4))*velocity;

%Create sawtooth motion of grating
t_profile=[t_profile;flipud(t_profile);t_profile;flipud(t_profile)];

stim=make_stim(x_profile,t_profile);

%We've created a short function that convolves a stimulus with a filter
%and creates images of the stimulus and the reponse:  

figure(3);
resp=plot_stim_and_resp(stim,rightward_1);

%And now we'll plot the time course of the response for a 
%centered filter again.

line([column,column],[0,length(t_stim)],'Color','r');
figure(4)
plot(resp(:,column));

%Ideal complex cells are motion selective but do not respond to the 
%spatial phase of a grating.  Adelson and Bergen show how to build
%spatial phase independent responses by squaring and adding the
%responses of pairs of even and odd direction selective filters.

%Let 'resp_right1' be the response of the 'rightward_1' filter:
resp_right_1=conv2(stim,rightward_1,'same');
%'resp_right2' is the response of the 'rightward_2' filter:
resp_right_2=conv2(stim,rightward_2,'same');

%We now square and sum the two responses:
energy_right=resp_right_1.^2+resp_right_2.^2;

figure(3)
subplot(1,2,2)
show_image(energy_right);

%And the time course of the centered filter is:
figure(4)
plot(energy_right(:,column));

%The response of this 'rightward energy mechanism' still modulates 
%with the phase of the stimulus.  This is due to imperfections in
%the underlying linear filters. 

%In the paper, Adelson and Bergen go on to combine the rightward and 
%leftward energy mechanisms to produce an 'opponent energy' mechanism.
%They also discuss how energy mechanisms filter out velocity 
%discontinuities inherent in movies, televisions, and computer monitors.

%There are a number of ways to go from here.  For example, what is 
%the response of an energy mechanism to moving dot stimuli?  Or, what
%is the spatial and temporal frequency sensitivity of such an energy
%model.  Is this consistent with single-cell physiology of complex
%cells?  Finally, we know that this model cannot predict the
%response to stimuli with multiple components.  For example, the response
%of a complex cell selective to rightward motion is suppressed by leftward
%motion, even though leftward alone motion produces little or no response.
%David Heeger proposes that this is the result of inhibitory effects
%across mechanisms.  The energy selective mechanisms in the Adelson-Bergen
%model form the inputs of Heeger's 'Divisive Inhibition' model.  Can this
%tutorial be modified to predict the response to multiple-component stimuli?



