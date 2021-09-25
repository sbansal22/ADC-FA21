%% Analog Lab
% @authors: Thomas Jagielski and Sparsh Bansal
% September 2021

%%
%%%%% Exercise 1 %%%%%
%% c)
fs = 300e3; % this is the sample rate
fc = 90.8e6; % this is the center frequency

x = zeros(3e6,1); % empty vector to store data

% create object for RTL-SDR receiver
rx = comm.SDRRTLReceiver('CenterFrequency',fc, 'EnableTunerAGC', false, 'TunerGain', 35,  'SampleRate', fs);

counter = 1; % initialize a counter
while(counter < length(x)) % while the buffer for data is not full
    rxdata = rx();   % read from the RTL-SDR
    x(counter:counter + length(rxdata)-1) = rxdata; % save the samples returned
    counter = counter + length(rxdata); % increment counter
end
% the data are returned as complex numbers
% separate real and imaginary part, and remove any DC offset
y_I = real(x)-mean(real(x));
y_Q = imag(x)-mean(imag(x));

figure
plot_FT(y_I, fs);
title('y_{I} in the Frequency Domain Recorded with RTL-SDR, f_{c} = 90.8e6 Hz and f_{s} = 300e3 s^{-1}')
xlabel('Frequency [Hz]')
ylabel('Magnitude y_{I}')

%% d)
time = (1:max(length(y_I)))/fs;

zoom_bound_upper = max(length(y_I))/2 - 2000;
zoom_bound_lower = max(length(y_I))/2 - 2500;

figure
subplot(2,1,1)
hold on
title('Recieved Signal Y_{I}')
xlabel('Time [s]')
ylabel('Magnitude Y_{I}')
plot(time, y_I)
hold off
subplot(2,1,2)
hold on
title('Recieved Signal Y_{I} Zoomed In')
xlabel('Time [s]')
ylabel('Magnitude Y_{I}')
plot(time(zoom_bound_lower:zoom_bound_upper), y_I(zoom_bound_lower:zoom_bound_upper))
hold off

figure
plot(time(zoom_bound_lower:zoom_bound_upper), y_I(zoom_bound_lower:zoom_bound_upper))
title('Recieved Signal Y_{I} Zoomed In')
xlabel('Time [s]')
ylabel('Magnitude Y_{I}')

%% e)
processed_y_i = diff(y_I);
processed_y_i(processed_y_i<0) = 0;
processed_y_i = (processed_y_i)./(max(abs(processed_y_i)));

time = (1:max(length(processed_y_i)))/fs;

figure
plot(time, processed_y_i)
title('Normalized Differentiated y_{I}')
xlabel('Time [s]')
ylabel('Magnitude of Processed y_{I}')

figure
plot(time(zoom_bound_lower:zoom_bound_upper), processed_y_i(zoom_bound_lower:zoom_bound_upper))
title('Normalized Differentiated y_{I}')
xlabel('Time [s]')
ylabel('Magnitude of Processed y_{I}')

%% f)
%processed_y_i = notch_filter(processed_y_i, 0, fs);

t = [-10000:-1, 1:10000];
%omega = 2 * pi * 0.095e6 * (1/fs);
omega = 2 * pi * 0.1e6 * (1/fs);
lpf = sin(omega*t)./(pi*t);

filtered = conv(processed_y_i, lpf);
%filtered = conv(filtered, lpf);

%figure
%subplot(2,1,1)
%plot(y_I)
%hold on
%subplot(2,1,2)
%plot(filtered)
%hold off

figure
subplot(2,1,1)
plot_FT(y_I, fs);
hold on
title('Frequency Domain of Normalized d/dt[y_{I}] Signal')
xlabel('Frequency [Hz]')
ylabel('Magnitude y_{I}')
hold off
subplot(2,1,2)
hold on
title('Frequency Domain of Normalized d/dt[y_{I}] Signal After LPF')
xlabel('Frequency [Hz]')
ylabel('Magnitude y_{I}')
plot_FT(filtered, fs);
axis([-1.5e5 1.5e5 0 7000])
hold off

%% g)
%filtered = filtered(30000:end);
centered_signal = filtered - mean(filtered);
processed_signal = (centered_signal * 3) ./ (max(abs(centered_signal)));
processed_signal = decimate(processed_signal, 4);

time = (1:max(length(processed_signal)))/fs;

figure
plot(time, processed_signal)
title('Filtered, Normalized, and Decimated y_{I} Signal')
xlabel('Time [s]')
ylabel('Magnitude Y_{I}')

sound(processed_signal, 300000/4)

%%
%%%%% Exercise 2 %%%%%
%% b)
fs = 300e3; % this is the sample rate
%fc = 107.9e6; % this is the center frequency
fc = 97.7e6;

x = zeros(3e6,1); % empty vector to store data

% create object for RTL-SDR receiver
rx = comm.SDRRTLReceiver('CenterFrequency',fc, 'EnableTunerAGC', false, 'TunerGain', 35,  'SampleRate', fs);

counter = 1; % initialize a counter
while(counter < length(x)) % while the buffer for data is not full
    rxdata = rx();   % read from the RTL-SDR
    x(counter:counter + length(rxdata)-1) = rxdata; % save the samples returned
    counter = counter + length(rxdata); % increment counter
end
% the data are returned as complex numbers
% separate real and imaginary part, and remove any DC offset
y_I = real(x)-mean(real(x));
y_Q = imag(x)-mean(imag(x));

diff_y_i = diff(y_I);
diff_y_q = diff(y_Q);
%%
message = (diff_y_q .* y_I(2:end)) - (diff_y_i .* y_Q(2:end));

figure
plot_FT(message, fs);
title('y_{I} in the Frequency Domain Recorded with RTL-SDR, f_{c} = 107.9e6 Hz and f_{s} = 300e3 s^{-1}')
xlabel('Frequency [Hz]')
ylabel('Magnitude of Message')

processed_message = diff(message);
processed_message(processed_message<0) = 0;
processed_message = (processed_message)./(max(abs(processed_message)));

time = (1:max(length(processed_message)))/fs;

figure
plot(time, processed_message)
title('Normalized Differentiated Message')
xlabel('Time [s]')
ylabel('Magnitude of Processed Message')

%figure
%plot(time(zoom_bound_lower:zoom_bound_upper), processed_y_i(zoom_bound_lower:zoom_bound_upper))
%title('Normalized Differentiated Message')
%xlabel('Time [s]')
%ylabel('Magnitude of Processed Message')

t = [-10000:-1, 1:10000];
omega = 2 * pi * 0.095e6 * (1/fs);
%omega = 2 * pi * 0.1e6 * (1/fs);
lpf = sin(omega*t)./(pi*t);

processed_message = notch_filter(processed_message, 0, fs);

filtered = conv(processed_message, lpf);
%filtered = conv(filtered, lpf);


figure
subplot(2,1,1)
plot_FT(message, fs);
hold on
title('Frequency Domain of Message Signal')
xlabel('Frequency [Hz]')
ylabel('Magnitude y_{I}')
hold off
subplot(2,1,2)
hold on
title('Frequency Domain of Normalized Message Signal After LPF')
xlabel('Frequency [Hz]')
ylabel('Magnitude y_{I}')
plot_FT(processed_message, fs);
%axis([-1.5e5 1.5e5 0 7000])
hold off


centered_signal = filtered - mean(filtered);
processed_signal = (centered_signal * 3) ./ (max(abs(centered_signal)));
processed_signal = decimate(processed_signal, 4);

time = (1:max(length(processed_signal)))/fs;

figure
plot(time, processed_signal)
title('Filtered, Normalized, and Decimated y_{I} Signal')
xlabel('Time [s]')
ylabel('Magnitude Y_{I}')

sound(processed_signal, 300000/4)
