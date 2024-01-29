%% Plot the signal

% INPUT:
% signal = ECG signal
% length = the length of ECG signal
% i = row index in the plot

% OUTPUT:
% figure = row plot

function [figure] = plot_signal(signal, length, i)
    figure = subplot(3,1,i);
    plot(signal, 'g');
    grid on;
    xlim([0 length]);
    xlabel('Time (sec)');
    ylabel('amplitude (Hz)');