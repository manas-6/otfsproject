clc;
clear all;
%using pulse shape ieee paper and otfs first principles paper
%subcarriers
M=4;

N=4;
channeltaps=2;
averaged_channels=1;
tp=15
for tp=100:100:1000
transmitpower=tp;
T=1;
TotalTime=N*T;

totalsymbols=M*N;
numBits = M*N*2; % Change this value based on your requirements

bits = randi([0, 1], 1, numBits);
%t=0:(1/((numBits/20)+channeltaps-1)):1-(1/((numBits/20)+channeltaps-1));
% Step 2: Map binary data to QPSK symbols
symbols = zeros(1, numBits/2); % Since QPSK maps 2 bits to each symbol
for i = 1:2:numBits
    symbols((i+1)/2) = 2*bits(i) + bits(i+1);
end
 
% Define QPSK modulation symbols
qpskSymbols = exp(1j * pi/4 * (1:2:7));

% Map symbols to QPSK modulation symbols
modulatedSymbols = qpskSymbols(symbols + 1);

X=zeros(M,N);
for k=1:M
for l=1:N
    X(k,l)=modulatedSymbols((k-1)*N+l);
end
end

%for 2ghz fc we get doppler shift as 200hz for vel=30;
%and for vel=150 doppler shift is 1khz;
h_=zeros(1,channeltaps);
v=[2/N*T 1/N*T 3/N*T];
tau=[(T/M) 2*(T/M) 3*(T/M)];
k_i=[1 2];
l_i=[2 1];

h_eff=zeros(M,N);

%for ch_av=1:averaged_channels
averaged_BER=0;

for i=1:channeltaps
    h(i)=1/(sqrt(2))*((randn)+1i*(randn));
    h_eff(k_i(i)+1,l_i(i)+1)=h(i);
end


xdd = reshape(X', [], 1);
Hdd=zeros(M*N,M*N);
for k1=0:M-1
    for k=0:M-1
        for l1=0:N-1
            for l=0:N-1
                
                n=floor((k1-k)/M);
                m=floor((l1-l)/N);
                Hdd(k1*N+l1+1,k*N+l+1)=h_eff(k1-k+1-(n*M),l1-l+1-(m*N))*exp(1j * 2 * pi *n*(l/N))*exp(1j * 2 * pi *(l1-l-(m*N))*(k+(n*M))*(1/(N*M)));
          
            end
        end
    end
end



for j=1:(numBits/2)
noise(j)=1/(sqrt(2))*((randn)+1i*(randn));
end

r=(Hdd*xdd)+((1/(sqrt(transmitpower)))*(noise'));

Y=inv(Hdd)*r

y = reshape(Y,[],1);

for i=1:(numBits/2)
error=y(i)-qpskSymbols;
leastmagnitude=min(abs(error));
index=find((abs(error))==leastmagnitude);
qpskoutput(i)=qpskSymbols(index);
qpskoutputindex(i)=index;
bitnumbers(j)=index(1)-1;

end
mismatches = ((round(qpskoutput.',5) ~= round(xdd,5)));

% Count the number of mismatches
numMismatches = sum(mismatches)
% Initialize an empty string to store binary strings
% binaryStrings = '';
% 
% % Convert each number to binary and concatenate to form bitstream
% for k = 1:length(bitnumbers)
%     % Convert current number to binary string
%     binaryString = dec2bin(bitnumbers(k),2);
%     % Concatenate binary string to the existing bitstream
%     binaryStrings = [binaryStrings, binaryString];
% end
% 
% % Convert bitstream string to numeric vector of bits
% bitstream = binaryStrings - '0';
% differentbits=nnz(bitstream~=bits)
% differentelements=nnz(qpskoutput~=modulatedSymbols)
% 
% averaged_BER=averaged_BER+differentbits;
% 
% end
% bitsunequal(tp/100)=averaged_BER/averaged_channels;
 end
SNR_dB = 10 * log10(100:100:1000)
figure;
plot(SNR_dB,bitsunequal)
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('SNR vs BER');
