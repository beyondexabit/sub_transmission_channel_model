function [r] = NoiseInsertion(x,OSNRdB,P)

%     fprintf('\n -------------- Noise Insertion ---------------\n')
    
    % OSNR in linear scale:
    OSNR = 10^(OSNRdB/10);
    
    % Average signal power:
    if isfield(P,'NTrib')
        temp = sum(abs(x).^2,1)./P.NTrib;
    else
        warning('Assumed data is transmitted in all tributaries, if not please enter P.NTrib')
    temp = mean(abs(x).^2,1); %mean(abs(x).^2,1);
    end
    AvgSigPower = 2*mean(temp); %mean(temp(temp > max(temp)/3));

   
%     AvgSigPower = 2*sum(mean(abs(x).^2,1))/P.NTrib;
    
    % Noise power:
    NoisePower = (P.DAC.Fs/P.Bref)*AvgSigPower/(2*OSNR);

%     OSNR = (p*Rs/2*Bref)*SNR
% 
%     SNR = (2*Bref)/(p*Rs)*OSNR

    
    
    % AWGN standard deviation:
    StdDev = sqrt(NoisePower/2);
    
    % AWGN generation:
    n = StdDev.*(randn(size(x))+1i*randn(size(x)));
    
    % Inserting noise to the signal:    
    r = x + n;
    
%     fprintf('\n Done. ----------------------------------------\n')
    
end