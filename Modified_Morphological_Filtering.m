%% ECG signal conditioning by morphological filtering
% A modified morphological filtering (MMF) algorithm is proposed for
% baseline correction and noise suppression of ECG signals.

%% Input
files = dir('Data\*.mat');  %load all dataset files
numData = numel(files);     %number of files
frequency = 3600;


%% For Baseline Drift
% SINUSOIDAL SIGNAL
x = linspace(0, 2*pi, frequency);
A = 0.8;    % A = 0.8 | 0.2;
N = 60;     % N = 60 | 40;
sinusoidal_signal = A*cos(x./N);

% SLANTED LINE
minDrift = 0;
maxDrift = 50;  % maxDrift = 50 | 30;
slanted_line = linspace(minDrift, maxDrift, frequency);

% Baseline drift is simulated by adding a slanted line to a sinusoidal signal
BaselineDrift = sinusoidal_signal + slanted_line;

%% For baseline correction
% Bo and Bc are selected as two horizontal line segments of zero amplitude, but with different lengths.
Lo = 0.2*frequency;
Lc = 1.5*Lo;
Bo = strel("line",Lo,0);
Bc = strel("line",Lc,0);

%% For noise suppression
B1 = [0 1 1 1 0]; B2 = [1 1 1 1 1];  %two structuring elements


%% evaluation matrix
Ev = zeros(numData*2, 3);    %prealloc
j = 1;                       %row index


%% Algorithm

for f = 1:numData
    load(fullfile("Data\",files(f).name)); % load all data from current file

    for i = 1:2
        ECG = val(1,:);

%% (1) Plot the Input signal
        fig1 = figure(1);
        plot_signal(ECG, frequency, 1);
        title('Original ECG signal');

%% Add Noise to the ECG signal
        Noise_Ecg = awgn(ECG, 20, 'measured');      %white Gaussian noise   param 20 | 30

    % Plot Noisy signal
        plot_signal(Noise_Ecg, frequency, 2);
        title('ECG signal with White Gaussian Noise');

%% Add Baseline Drift to ECG signal
        Dirty_Signal = Noise_Ecg + BaselineDrift;

    % Plot Dirty signal
        plot_signal(Dirty_Signal, frequency, 3);
        title('Noise + Baseline Drift');

        %saveas(fig1, 'Dirty_Signal.fig');

%% (2) Plot the Dirty signal
        fig2 = figure(2);
        plot_signal(Dirty_Signal, frequency, 1);
        title('Dirty ECG signal');

%% Baseline correction
    % The signal is first opened by a structuring element Bo for removing peaks in the signal.
        peaks_suppression = imopen(Dirty_Signal, Bo);

    % Then the resultant waveforms with pits are removed by a closing operation
    % using the other structuring element Bc.
        pits_suppression = imclose(peaks_suppression, Bc);

    % The final result is then an estimate of the baseline drift.
    % The correction of the baseline is then done by subtracting the result from the original signal.
        Corrected_Signal = Dirty_Signal - pits_suppression;


    % Plot the corrected signal (Baseline)
        plot_signal(Corrected_Signal, frequency, 2);
        title('Corrected signal (Baseline)');

%% Noise suppression
    % Apply MMF algorithm

    % Opening operation
        erosion = imerode(Corrected_Signal, B2);
        opening = imdilate(erosion, B1);

    % Closing operation
        dilatation = imdilate(Corrected_Signal, B1);
        closing = imerode(dilatation, B2);

        MMF = (opening + closing)/2; % average

    % Plot the reconstructed signal
        plot_signal(MMF, frequency, 3);
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
        SDR = sum(abs(ECG - MMF)) / sum(abs(MMF))
        Ev(j,3) = SDR;

        j = j+1;
    end

end


%% Save the results into a CSV file
Names = {'BCR', 'NSR', 'SDR'};
C = [Names; num2cell(Ev)];

writecell(C,'Evaluation of signal conditioning.csv') 
