%% Plot the signal

% INPUT:
% signal = ECG signal
% frequency = the length of ECG signal
% i = row index in the plot

% OUTPUT:
% figure = row plot

function [figure] = plot_signal(signal, frequency, i)
    figure = subplot(3,1,i);
    plot(signal,'g');
    grid on;
    xlim([0 frequency]);
    xlabel('f (Hz)');
    ylabel('amplitude (mV)');