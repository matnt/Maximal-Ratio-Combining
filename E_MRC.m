% 4. SIMO Rayleigh fading channel employing MRC with 2, 3, 4 antenna perfect CSIR 
clear
clc

N = 100000; % number of symbols
Es_N0 = 0: 1: 30; % SNR

% generate transmitted signal
% 1 + i, 1 - i, -1 + i, -1 - i
Tx = (2*(rand(1, N) > 0.5) - 1) + j * (2*(rand(1, N) > 0.5) - 1);
Tx_normal = Tx/sqrt(2); % normalized signal
n = (randn(4, N) + j * randn(4, N))/sqrt(2); % noise
h = (randn(4, N) + j * randn(4, N))/sqrt(2); % rayleigh channel

Rx = zeros(3, N); % received signal before detection
Rx_detection = zeros(3, N); % received signal after detection
Rx_output = zeros(3, N); % output signal

count_error = zeros(3, length(Es_N0)); % count number of error symbols
P_s = zeros(3, length(Es_N0)); % probability of symbol error
for i = 1: length(Es_N0)
    
    Rx = h .* (10^(Es_N0(i)/20) * Tx_normal) + n;
    
    for k = 1: 3
       Rx_detection(k, :) = sum(conj(h(1: (k + 1), :)) .* Rx(1: (k + 1), :), 1) ./ sqrt(sum(conj(h(1: (k + 1), :)) .* h(1: (k + 1), :), 1));
       
       Rx_output(k, find(real(Rx_detection(k, :)) > 0 & imag(Rx_detection(k, :)) > 0 )) = 1 + j; 
       Rx_output(k, find(real(Rx_detection(k, :)) > 0 & imag(Rx_detection(k, :)) < 0 )) = 1 - j; 
       Rx_output(k, find(real(Rx_detection(k, :)) < 0 & imag(Rx_detection(k, :)) > 0 )) = -1 + j; 
       Rx_output(k, find(real(Rx_detection(k, :)) < 0 & imag(Rx_detection(k, :)) < 0 )) = -1 - j; 

       count_error(k, i) = size(find(Tx - Rx_output(k,:)), 2);
       P_s(k, i) = count_error(k, i)/N;
    end
    
end
figure(5)

semilogy(Es_N0,P_s(1,:),'-c+')
hold on
semilogy(Es_N0,P_s(2,:),'-co')
semilogy(Es_N0,P_s(3,:),'-cd')
grid on
xlabel('E_s/N_0(dB)')
ylabel('Probability of symbol error (%)')
hold off
legend('MRC 2 antennas','MRC 3 antennas','MRC 4 antennas')
axis([1 30 10^-3 1])
