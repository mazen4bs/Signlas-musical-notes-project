f0 = 440;
alpha = 2^(1/12);

n = 0;
fn = f0 * alpha^n;

n_DO = -9;
fn_DO = f0 * alpha^n_DO;
n_CMAJOR = [-9, -7, -5, -4];

F = [];

fprintf('C-Major Scale Frequencies:\n');
for i = 1:length(n_CMAJOR)
    note = '';
    switch n_CMAJOR(i)
        case -9
            note = 'DO';
        case -7
            note = 'RE';
        case -5
            note = 'MI';
        case -4
            note = 'FA';
        otherwise
            note = 'Unknown';
    end

    frequency = f0 * alpha^n_CMAJOR(i);

    f(i) = frequency;

    fprintf('%s: %.4f Hz\n', note, frequency);
end


f_Do = f(1);
f_RE = f(2);
f_MI = f(3);
f_FA = f(4);

fs = 10 * f_FA;

duration = 0.5;
t=0:1/fs:duration-1/fs;

x1 = cos(2*pi*fn_DO*t);
x2 = cos(2*pi*f0*alpha^(-7)*t);
x3 = cos(2*pi*f0*alpha^(-5)*t);
x4 = cos(2*pi*f0*alpha^(-4)*t);

xt = [x1, x2, x3, x4];

audiowrite('DoRiMeFa.wav', xt, fs);

t_seq = linspace(0, 2, numel(xt));

figure;
plot(t_seq, xt)
title('Signal x(t)')
xlabel('Time')
ylabel('Amplitude')

energy_xt = sum(abs(xt).^2/fs);

fprintf('Energy: %.4f\n', energy_xt);

Xf = fft(xt);

L = length(xt);
F = (-fs/2:fs/L:fs/2-fs/L);
Xf_ = fftshift(Xf);

figure;
plot(F, abs(Xf_)/L);
title('Frequency Spectrum of X(f)')
xlabel('Frequency (Hz)')
ylabel('Magnitude')

energy_Xf = (sum(abs(Xf).^2) / length(Xf))/fs;

fprintf('Energy in Frequency Domain: %.4f\n', energy_Xf);


pkg load signal;

cutoff_Freq = (f_RE + f_MI)/2;

[b, a] = butter(20, cutoff_Freq / (fs / 2), 'low');

figure;
freqz(b, a, [], fs);

y1_t = filter(b, a, xt);

audiowrite('DoRiMeFa_LPF.wav', y1_t, fs);

figure;
subplot(2,1,1);
plot(t_seq, xt);
title('Original Signal x(t)');
xlabel('Time');
ylabel('Amplitude');

subplot(2,1,2);
plot(t_seq, y1_t);
title('Filtered Signal y1(t)');
xlabel('Time');
ylabel('Amplitude');

energy_y1_t = sum(abs(y1_t).^2)/fs;
fprintf('Energy of y1_t: %.4f\n', energy_y1_t);

Y1_f = fft(y1_t);
L_y1 = length(y1_t);
F_y1 = (-fs/2:fs/L_y1:fs/2-fs/L_y1);
Y1f_ = fftshift(Y1_f);

figure;
plot(F_y1, abs(Y1f_)/L_y1);
title('Frequency Spectrum of y1(t)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

energy_Y1_f = (sum(abs(Y1_f).^2) / length(Y1_f))/fs;
fprintf('Energy in Frequency Domain (y1_t): %.4f\n', energy_Y1_f);


[d , c] = butter(20, cutoff_Freq / (fs / 2), 'high');

figure;
freqz(d, c, [], fs);

y2_t = filter(d, c, xt);

audiowrite('DoRiMeFa_HPF.wav', y2_t, fs);

figure;
subplot(2,1,1);
plot(t_seq, xt);
title('Original Signal x(t)');
xlabel('Time');
ylabel('Amplitude');

subplot(2,1,2);
plot(t_seq, y2_t);
title('High-pass Filtered Signal y2(t)');
xlabel('Time');
ylabel('Amplitude');

energy_y2_t = sum(abs(y2_t).^2) / fs;

Y2_f = fft(y2_t);

L_highpass = length(y2_t);
F_highpass = (-fs/2:fs/L_highpass:fs/2-fs/L_highpass);
Y2f_ = fftshift(Y2_f);

figure;
plot(F_highpass, abs(Y2f_)/L_highpass);
title('Frequency Spectrum of y2(t)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

fprintf('Energy of y2_t: %.4f\n', energy_y2_t);

energy_Y2_f = (sum(abs(Y2_f).^2) / length(Y2_f)) / fs;

fprintf('Energy in Frequency Domain (y2_t): %.4f\n', energy_Y2_f);
