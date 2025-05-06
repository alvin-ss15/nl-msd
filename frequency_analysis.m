% Perform FFT analysis
Fs = 1/mean(diff(t));  % Average sampling frequency
L = length(t);         % Length of signal

% Prepare for FFT
% Zero-pad for better frequency resolution
n_fft = 2^nextpow2(L);  % Next power of 2 from length of signal

% Calculate FFT for each mass position
X1_fft = fft(x(:,1), n_fft);
X2_fft = fft(x(:,2), n_fft);
X3_fft = fft(x(:,3), n_fft);

% Get single-sided spectrum
P1_X1 = abs(X1_fft/L);
P1_X1 = P1_X1(1:n_fft/2+1);
P1_X1(2:end-1) = 2*P1_X1(2:end-1);

P1_X2 = abs(X2_fft/L);
P1_X2 = P1_X2(1:n_fft/2+1);
P1_X2(2:end-1) = 2*P1_X2(2:end-1);

P1_X3 = abs(X3_fft/L);
P1_X3 = P1_X3(1:n_fft/2+1);
P1_X3(2:end-1) = 2*P1_X3(2:end-1);

% Define frequency vector
f = Fs*(0:(n_fft/2))/n_fft;

% Find peak frequencies using findpeaks (for non-zero frequencies)
[pks1, locs1] = findpeaks(P1_X1, 'MinPeakHeight', max(P1_X1)*0.1, 'SortStr', 'descend');
[pks2, locs2] = findpeaks(P1_X2, 'MinPeakHeight', max(P1_X2)*0.1, 'SortStr', 'descend');
[pks3, locs3] = findpeaks(P1_X3, 'MinPeakHeight', max(P1_X3)*0.1, 'SortStr', 'descend');

% Check DC component (frequency = 0) and add it if significant
dc_threshold = max(P1_X1)*0.1; % Same threshold as for other peaks

% For Mass 1
if P1_X1(1) >= dc_threshold
    pks1 = [P1_X1(1); pks1];
    locs1 = [1; locs1];
end

% For Mass 2
if P1_X2(1) >= dc_threshold
    pks2 = [P1_X2(1); pks2];
    locs2 = [1; locs2];
end

% For Mass 3
if P1_X3(1) >= dc_threshold
    pks3 = [P1_X3(1); pks3];
    locs3 = [1; locs3];
end
% Calculate natural frequencies (linear approximation)
% Mass matrix
M = eye(3)*m;

% Stiffness matrix
K = zeros(3);
K(1,1) = 2*k1; K(1,2) = -k1;
K(2,1) = -k1; K(2,2) = 2*k1; K(2,3) = -k1;
K(3,2) = -k1; K(3,3) = k1;

% Damping matrix (for reference only)
C = zeros(3);
C(1,1) = 2*c1; C(1,2) = -c1;
C(2,1) = -c1; C(2,2) = 2*c1; C(2,3) = -c1;
C(3,2) = -c1; C(3,3) = c1;

% Calculate eigenvalues for undamped system
[V, D] = eig(K, M);
natural_frequencies = sqrt(diag(D))/(2*pi);  % Convert to Hz

% Plot frequency spectra with peak markers and natural frequencies
figure;
subplot(3,1,1);
plot(f, P1_X1);
hold on;
% Mark detected peaks
plot(f(locs1), pks1, 'r.', 'MarkerSize', 8);
% Mark natural frequencies
for i = 1:length(natural_frequencies)
    xline(natural_frequencies(i), '--g', ['f_n' num2str(i)], 'LineWidth', 1.5);
    text(natural_frequencies(i), 0, [num2str(natural_frequencies(i), '%.2f')], ...
         'Color', 'b', 'HorizontalAlignment', 'center');
end
title('Mass 1 Frequency Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;
xlim([0, 5]);

subplot(3,1,2);
plot(f, P1_X2);
hold on;
% Mark detected peaks
plot(f(locs2), pks2, 'r.', 'MarkerSize', 8);
% Mark natural frequencies
for i = 1:length(natural_frequencies)
    xline(natural_frequencies(i), '--g', ['f_n' num2str(i)], 'LineWidth', 1.5);
    text(natural_frequencies(i), 0, [num2str(natural_frequencies(i), '%.2f')], ...
         'Color', 'b', 'HorizontalAlignment', 'center');
end
title('Mass 2 Frequency Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;
xlim([0, 5]);

subplot(3,1,3);
plot(f, P1_X3);
hold on;
% Mark detected peaks
plot(f(locs3), pks3, 'r.', 'MarkerSize', 8);
% Mark natural frequencies
for i = 1:length(natural_frequencies)
    xline(natural_frequencies(i), '--g', ['f_n' num2str(i)], 'LineWidth', 1.5);
    text(natural_frequencies(i), 0, [num2str(natural_frequencies(i), '%.2f')], ...
         'Color', 'b', 'HorizontalAlignment', 'center');
end
title('Mass 3 Frequency Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;
xlim([0, 5]);

% Plot spectrogram to see frequency evolution over time
figure;
window_size = 256;
overlap = 250;
nfft = 1024;

subplot(3,1,1);
spectrogram(x(:,1), hamming(window_size), overlap, nfft, Fs, 'yaxis');
title('Mass 1 Spectrogram');
colorbar;

subplot(3,1,2);
spectrogram(x(:,2), hamming(window_size), overlap, nfft, Fs, 'yaxis');
title('Mass 2 Spectrogram');
colorbar;

subplot(3,1,3);
spectrogram(x(:,3), hamming(window_size), overlap, nfft, Fs, 'yaxis');
title('Mass 3 Spectrogram');
colorbar;

% For linear case
% Calculate the system matrix for the linearized system
M = eye(3)*m;  % Mass matrix
K = zeros(3);  % Stiffness matrix
C = zeros(3);  % Damping matrix

% Fill stiffness matrix (spring connections)
K(1,1) = 2*k1; K(1,2) = -k1; 
K(2,1) = -k1; K(2,2) = 2*k1; K(2,3) = -k1;
K(3,2) = -k1; K(3,3) = k1;

% Fill damping matrix (damper connections)
C(1,1) = 2*c1; C(1,2) = -c1; 
C(2,1) = -c1; C(2,2) = 2*c1; C(2,3) = -c1;
C(3,2) = -c1; C(3,3) = c1;

% Calculate eigenvalues of the system - undamped
[V, D] = eig(K, M);
natural_frequencies = sqrt(diag(D))/(2*pi);  % Convert to Hz

% Display theoretical natural frequencies
disp('Theoretical natural frequencies (Hz):');
disp(natural_frequencies);

% Compare with peaks in FFT plots
% In nonlinear system, these will shift based on amplitude