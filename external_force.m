% Create a cell array to store all force functions
F_external = cell(30, 1);

% Step Inputs (Various amplitudes to test nonlinearity)
F_external{1} = @(t) 0.1 * (t >= 10); % Small step (linear regime)
F_external{2} = @(t) 1.0 * (t >= 10); % Medium step
F_external{3} = @(t) 5.0 * (t >= 10); % Large step (nonlinear regime)
F_external{4} = @(t) -0.1 * (t >= 10); % Small negative step
F_external{5} = @(t) -5.0 * (t >= 10); % Large negative step

% Impulse Inputs (for transient response analysis)
F_external{6} = @(t) 1.0 * (t >= 5 & t <= 5.1); % Short impulse (0.1s)
F_external{7} = @(t) 5.0 * (t >= 5 & t <= 5.1); % Strong short impulse
F_external{8} = @(t) 1.0 * (t >= 5 & t <= 5.5); % Long impulse (0.5s)
F_external{9} = @(t) -1.0 * (t >= 5 & t <= 5.1); % Negative impulse
F_external{10} = @(t) 1.0 * (t >= 5 & t <= 5.1) + -1.0 * (t >= 15 & t <= 15.1); % Double impulse

% Sinusoidal at Natural Frequencies (for resonance testing)
F_external{11} = @(t) 0.5 * sin(2*pi*0.2240*t); % First natural frequency (0.22 Hz) - small amplitude
F_external{12} = @(t) 2.0 * sin(2*pi*0.2240*t); % First natural frequency - large amplitude
F_external{13} = @(t) 0.5 * sin(2*pi*0.6276*t); % Second natural frequency (0.63 Hz) - small amplitude
F_external{14} = @(t) 2.0 * sin(2*pi*0.6276*t); % Second natural frequency - large amplitude
F_external{15} = @(t) 0.5 * sin(2*pi*0.9069*t); % Third natural frequency (0.91 Hz) - small amplitude
F_external{16} = @(t) 2.0 * sin(2*pi*0.9069*t); % Third natural frequency - large amplitude

% Frequency Sweep (for frequency response analysis)
F_external{17} = @(t) 1.0 * sin(2*pi * (0.1 + 1.5*t/60) * t); % Sweep from 0.1 Hz to 1.6 Hz
F_external{18} = @(t) 3.0 * sin(2*pi * (0.1 + 1.5*t/60) * t); % Stronger sweep (nonlinear effects)

% Near-Natural Frequencies (to test bandwidth)
F_external{19} = @(t) 1.0 * sin(2*pi*0.20*t); % Just below first natural frequency
F_external{20} = @(t) 1.0 * sin(2*pi*0.24*t); % Just above first natural frequency
F_external{21} = @(t) 1.0 * sin(2*pi*0.60*t); % Just below second natural frequency
F_external{22} = @(t) 1.0 * sin(2*pi*0.66*t); % Just above second natural frequency

% Multi-Frequency Inputs (for cross-coupling analysis)
F_external{23} = @(t) sin(2*pi*0.22*t) + sin(2*pi*0.63*t); % First + Second modes
F_external{24} = @(t) sin(2*pi*0.22*t) + sin(2*pi*0.91*t); % First + Third modes
F_external{25} = @(t) sin(2*pi*0.63*t) + sin(2*pi*0.91*t); % Second + Third modes
F_external{26} = @(t) sin(2*pi*0.22*t) + sin(2*pi*0.63*t) + sin(2*pi*0.91*t); % All three modes

% Specialized Signals
F_external{27} = @(t) (1 + 0.5*sin(2*pi*0.05*t)) .* sin(2*pi*0.63*t); % Amplitude-modulated at 2nd mode
F_external{28} = @(t) 1.0 * (t >= 10) + 0.5 * sin(2*pi*0.91*t) .* (t >= 10); % Step with oscillation

% Random Signals (for robustness testing)
rng(1); % For reproducibility
noise_values = 2*rand(10000,1) - 1; % Random noise between -1 and 1
F_external{29} = @(t) interp1(linspace(0,60,10000), noise_values, t, 'linear', 0); % Random noise

% Filtered Noise (around 2nd natural frequency)
filtered_noise = filter([0.2 0.2], [1 -0.6], noise_values); % Bandpass-like filter
F_external{30} = @(t) 2.0 * interp1(linspace(0,60,10000), filtered_noise, t, 'linear', 0); % Filtered noise

% Save to a .mat file
save('F_external.mat', 'F_external');