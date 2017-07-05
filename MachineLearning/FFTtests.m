
A = 1;
hz = 2;

Fs = 30;            % Sampling frequency
T = 1/Fs;             % Sampling period
L = 100;             % Length of signal
t = (0:L-1)*T;        % Time vector

S = A*sin(2*pi*hz*t);% + 0.01*randn(size(X));

figure
plot(t,S)

Y = fft(S);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure
f = Fs*(0:(L/2))/L;
plot(f,P1)
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')