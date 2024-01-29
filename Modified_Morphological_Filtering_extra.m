%% ECG signal conditioning by morphological filtering
% A modified morphological filtering (MMF) algorithm is proposed for
% baseline correction and noise suppression of ECG signals.

%% Input
files = dir('Data\*.mat');  %load all dataset files
numData = numel(files);     %number of files
size = 3600;                %length of ECG signals
fs = 360;                   %frequency [Hz]

% Butterworth filter
[b,a] = butter(3, [0.5 40]/(fs/2), "bandpass");


%% Noise parameter
snr = [10, 20, 30, 40, 50, 60, 70];                 % Signal-to-noise ratio in dB

%% Pacemaker spike parameter
spike_amplitude = [-45, -25, 5, 25, 45, 65, 85];      % Amplitude of the pacemaker spike

%% Baseline Drift parameters
% SINUSOIDAL SIGNAL
x = linspace(0, 2*pi, size);
A = [0.5, 32, 64, 96, 128, 160, 192];
bf = (14/60)*2*pi;                                  % respiratory rate at rest (rad/seconds)
c = [2, 30, 50, 60, 80, 90, 100];                   % constant

% SLANTED LINE
minDrift = 0;
maxDrift = [-150, -50, 10, 50, 100, 150, 200];

%% Baseline correction parameters
% Bo and Bc are selected as two horizontal line segments of zero amplitude, but with different lengths.
Lo = 0.2*fs;
Lc = 1.5*Lo;
Bo = strel("line",Lo,0);
Bc = strel("line",Lc,0);

%% Noise suppression parameters
B1 = [0 1 5 1 0]; B2 = [1 1 1 1 1];  %two structuring elements


%% evaluation matrix
Ev = zeros(numData*2, 3);    %prealloc
j = 1;                       %row index

%% matrix of averages
Average = zeros(7, 3);      %prealloc

for si = 1:7
    %si=3;
    ai=6;
    mi=3;
    ci=1;

    for f = 1:numData
        load(fullfile("Data\",files(f).name)); % load all data from current file

        for i = 1:2
            ECG = val(1,:);
            ECG = filtfilt(b,a,ECG);     % Zero-phase digital filtering

%% (1) Plot the Input signal
            fig1 = figure(1);
            plot_signal(ECG, size, 1);
            title('Original ECG signal');

%% Add Noise to the ECG signal
            Noise_Ecg = awgn(ECG, snr(si), 'measured');      %white Gaussian noise

    % Plot Noisy signal
            plot_signal(Noise_Ecg, size, 2);
            title('ECG signal with White Gaussian Noise');

%% Add pacemaker spikes to the ECG signal
         %Noise_Ecg = ECG;
         
         %for index = 1:36
          %  spike_location = round(index * 100);
          %  Noise_Ecg(spike_location) = Noise_Ecg(spike_location) + spike_amplitude(si);
         %end

    % Plot Noisy signal
         %plot_signal(Noise_Ecg, size, 2);
         %title('ECG signal with pacemaker spikes');

%% Add Baseline Drift to ECG signal
            sinusoidal_signal = A(ai)*cos(x*bf) + c(ci);
            slanted_line = linspace(minDrift, maxDrift(mi), size);

        % Baseline drift is simulated by adding a slanted line to a sinusoidal signal
            BaselineDrift = sinusoidal_signal + slanted_line;

            Dirty_Signal = Noise_Ecg + BaselineDrift;
            
    % Plot Dirty signal
            plot_signal(Dirty_Signal, size, 3);
            title('Noise + Baseline Drift');

            %saveas(fig1, 'Dirty_Signal.fig');

%% (2) Plot the Dirty signal
            fig2 = figure(2);
            plot_signal(Dirty_Signal, size, 1);
            title('Dirty ECG signal');

%% Baseline correction
    % Apply MMF algorithm
    
    % The signal is first opened by a structuring element Bo for removing peaks in the signal.
            peaks_suppression = imopen(Dirty_Signal, Bo);

    % Then the resultant waveforms with pits are removed by a closing operation
    % using the other structuring element Bc.
            pits_suppression = imclose(peaks_suppression, Bc);

    % The final result is then an estimate of the baseline drift.
    % The correction of the baseline is then done by subtracting the result from the original signal.
            Corrected_Signal = Dirty_Signal - pits_suppression;


    % Plot the corrected signal (Baseline)
            plot_signal(Corrected_Signal, size, 2);
            title('Corrected signal (Baseline)');

%% Noise suppression
    % Apply MMF algorithm

    % Opening operation
            erosion = erosion_function(Corrected_Signal, B2);
            opening = dilation_function(erosion, B1);

    % Closing operation
            dilation = dilation_function(Corrected_Signal, B1);
            closing = erosion_function(dilation, B2);

            MMF = (opening + closing)/2; % average

    % Plot the reconstructed signal
            plot_signal(MMF, size, 3);
            title('Reconstructed signal (Denoised)');

            %saveas(fig2, 'Reconstructed_Signal.fig');

%% Performance evaluation of signal conditioning
    % BCR Baseline correction ratio
    % is defined for measuring the degree of baseline being corrected
            BCR = sum(abs(Corrected_Signal)) / sum(abs(Noise_Ecg))
            Ev(j,1) = BCR;
    % NSR Nosie suppression ratio
    % is defined for measuring the degree of noise being suppressed
            NSR = sum(abs(MMF)) / sum(abs(ECG))
            Ev(j,2) = NSR;
    % SDR signal-distortion ratio
    % is defined for measuring the degree of signal being distorted after the conditioning
            SDR = sum(abs(ECG) - abs(MMF)) / sum(abs(MMF))
            Ev(j,3) = SDR;

            j = j+1;
        end

    end


%% Experimental results
    fig3 = figure(3);
    boxplot(Ev, 'Labels', {'BCR', 'NSR', 'SDR'});

% Calculate the average
    avg = mean(Ev,1)
    Average(si,1) = avg(1);
    Average(si,2) = avg(2);
    Average(si,3) = avg(3);

end

fig4 = figure(4);
plot(Average, '-*');
xticklabels(snr);
xlabel('snr');
yline(1,'--');
yline(0,'--')
legend('BCR', 'NSR', 'SDR');

s = std(Average)
