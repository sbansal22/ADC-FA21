function [output] = notch_filter(signal, filter_frequency, sampling_frequency)
    omega_0 = filter_frequency/sampling_frequency * 2 * pi; % Angular Frequency
    q = 0.999; 
    
    output = zeros(length(signal) + 2, 1);
    new_signal = signal;

    output(1) = new_signal(1);
    output(2) = new_signal(2);

    for k = 3:length(new_signal)
        output(k) = 2*q*cos(omega_0)*output(k-1)-q^2*output(k-2) +...
            new_signal(k) - 2*cos(omega_0)*new_signal(k-1) + new_signal(k-2);
    end
    
    figure
    subplot(2,1,1)
    plot_FT(signal, sampling_frequency);
    title("Signal without Filtering")
    xlabel("Frequency [Hz]")
    ylabel("|X\_a|")
       
    subplot(2,1,2)
    plot_FT(output, sampling_frequency);
    title("Signal with Filtering")
    xlabel("Frequency [Hz]")
    ylabel("|X\_a|")

end