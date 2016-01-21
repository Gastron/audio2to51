clear
%Load a sample:
[sig,Fs] = audioread('sample.wav');

%Windowing:
winlen = 2^nextpow2(0.1*Fs);
window = sqrt(hann(winlen,'periodic'));
hopsize = winlen/2;
zeroPad = mod(-length(sig),winlen);

%STFT:
x = [zeros(hopsize,2);sig;zeros(hopsize+zeroPad,2)];

Xlen = length(x)/hopsize - 2;
X1 = zeros(winlen,Xlen);
X2 = zeros(winlen,Xlen);
for ind = 1:Xlen
    X1(:,ind) = fft(window.*x(1+hopsize*ind:winlen+hopsize*ind,1),winlen);
    X2(:,ind) = fft(window.*x(1+hopsize*ind:winlen+hopsize*ind,2),winlen);
end

%Short-time coherence function:
lambda = 0.2;
phi_11 = zeros(size(X1));
phi_22 = zeros(size(X2));
phi_12 = zeros(size(X1));
phi_11(:,1) = lambda.*(X1(:,1).*conj(X1(:,1)));
phi_22(:,1) = lambda.*(X2(:,1).*conj(X2(:,1)));
phi_12(:,1) = lambda.*(X1(:,1).*conj(X2(:,1)));
for ind = 2:Xlen
    phi_11(:,ind) = (1-lambda).*phi_11(:,ind-1) + lambda.*(X1(:,ind).*conj(X1(:,ind)));
    phi_22(:,ind) = (1-lambda).*phi_22(:,ind-1) + lambda.*(X2(:,ind).*conj(X2(:,ind)));
    phi_12(:,ind) = (1-lambda).*phi_12(:,ind-1) + lambda.*(X1(:,ind).*conj(X2(:,ind)));
end
phi = abs(phi_12)./(sqrt(phi_11.*phi_22));

%Similarity function:
psi_11 = zeros(size(X1));
psi_22 = zeros(size(X2));
psi_12 = zeros(size(X1));
for ind = 1:Xlen
    psi_11(:,ind) = X1(:,ind).*conj(X1(:,ind));
    psi_22(:,ind) = X2(:,ind).*conj(X2(:,ind));
    psi_12(:,ind) = X1(:,ind).*conj(X2(:,ind));
end
psi = 2*abs(psi_12)./(psi_11+psi_22);


%Ambience extraction:
AmbIndex = 1 - phi;
mu1=1;
mu0=0.1;
thresh = 0.5;
sigma = 4;
AmbWeights = (mu1-mu0)/2*tanh(sigma*pi*(AmbIndex-thresh))+(mu1+mu0)/2;
A1 = AmbWeights.*X1;
A2 = AmbWeights.*X2;


%ISTFT:
a = zeros(length(x),2);
for ind = 1:Xlen
    a(1+hopsize*ind:winlen+hopsize*ind,1) = a(1+hopsize*ind:winlen+hopsize*ind,1) + ...
        window.*ifft(A1(:,ind));
    a(1+hopsize*ind:winlen+hopsize*ind,2) = a(1+hopsize*ind:winlen+hopsize*ind,2) + ...
        window.*ifft(A2(:,ind));
end
a = a(hopsize+1:end-hopsize-zeroPad,:);

%Spectrograms:
figure(1)
imagesc([1 length(x)/Fs],[0 Fs/2],20*log(abs(A1(1:round(winlen/2),:)))); 
set(gca,'YDir','normal')
title('Spectrogram of right ambience signal')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
figure(2)
imagesc([1 length(x)/Fs],[0 Fs/2],20*log(abs(X1(1:round(winlen/2),:)))); 
set(gca,'YDir','normal')
title('Spectrogram of right channel signal')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
