function plot_spectrogram(x, dt)
    %% 
    % x - signal
    % dt - integration time constant
    fs = 1000/dt;       % in Hz
    window = 3*fs;           % Divide the signal into sections
    n_overlap = 0.95*window; %window/20;   % Specify number of samples of overlap between adjoining sections
    dft_points = 5*window;   % Evaluate the spectrum at round(dft_points/2 +1) frequencies.
     
    figure;
    spectrogram(x, hamming(window), n_overlap, dft_points, fs, 'yaxis', 'psd');
 
    %%
    ylim([0, 0.012])
    
    % Get yRuler-object
    ax = gca;
    yRuler = ax.YAxis;

    % Loop through TickValues, multiply and insert into TickLabels cell-array
    for n = 1:numel(yRuler.TickValues)
        yRuler.TickLabels{n} = num2str(yRuler.TickValues(n) * 1000);
    end
    
    ylabel('Frequency (Hz)')
    title('Spectrogram')
    %%
    %yticks(yticks*1000)
    %yticks(1:12)
    %yticks
    %yticklabels
    %yticklabels('auto')
    %ylabel('Frequency (Hz)')
    %yticklabels({'0', '2', '4', '6', '8', '10', '12'})
    %yticklabels('auto')
    
    %yticklabels
    %ylabel('Frequency (Hz)')
    
    %ax = gca;
    %ax.YAxis.Limits = [0, 0.015];
    %ax.YAxis.TickValues = ax.YAxis.TickValues * 1000;
    %ax.YAxis.TickLabelFormat = '%.2f';
    %yl = yticklabels
    %ylabel('Frequency (Hz)')
    %yticklabels(ax.YAxis.TickValues)
    %yticklabels('manual')
    
    %yticks(new_yt) % multiply yticks by 1000 to visualize Hz
    %yticks('manual')
    %ylim([0, 0.015])
    %ylim([0, 15])
    
    %spectrogram(x, window, n_overlap, dft_points, fs, 'yaxis', 'psd');
    %spectrogram(x, [], [], [], fs, 'yaxis', 'psd');
end

 

    