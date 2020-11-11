function Xest = weiner_deconv(Y,h,sig2)

%These scripts are Written by
% Srikanth Ryali, PhD & Vinod Menon, PhD
%Department of Psychiatry & Behavioral Sciences
%Stanford Cognitive and Systems Neuroscience Laboratory
%Stanford School of Medicine
%Stanford, USA
%Released only for Research purposes


[M,N] = size(Y);
hext = zeros(1,N);
hext(1:length(h)) = h;
Hw = fft(hext);
Xest = zeros(M,N);
for m = 1:M
    y = Y(m,:);
    Yw = fft(y);
    Xest(m,:) = ifft((conj(Hw).*Yw)./(abs(Hw).^2 + sig2(m,m)));
end