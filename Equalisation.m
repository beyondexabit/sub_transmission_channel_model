function y = Equalisation(x,W,DFTsize,NOverlap,EqType,P)

% Norm. signal to unitary energy:
% x = reshape(x,[],size(x,2));
% x = bsxfun(@rdivide,x,sqrt(mean(abs(x).^2)));

% Extending the input signal so that the blocks can be properly formed:
Nsym   = size(x,1);
AuxLen = Nsym/(DFTsize-NOverlap);
if AuxLen ~= ceil(AuxLen)
    NExtra = ceil(AuxLen)*(DFTsize-NOverlap)-Nsym;
    x      = [x ; x(1:NExtra,:)];
else
    NExtra = NOverlap;
    x      = [zeros(NExtra/2,size(x,2));x;zeros(NExtra/2,size(x,2))];
end

% Input blocks (N. of samples or symols x N. tributaries x N. blocks):
Blocks = permute(reshape(x,DFTsize-NOverlap,size(x,1)/(DFTsize-NOverlap)...
    ,size(x,2)),[1,3,2]);

% Overlap = zeros(NOverlap,size(x,2));
Overlap = x(end-NOverlap+1:end,:);
if strcmp(EqType,'MISO')
    yBFreq  = zeros(DFTsize,1);

    % Preallocating the output blocks and the overlap for the first block:
    y = zeros(size(Blocks,1)*size(Blocks,3),1);
else
    yBFreq  = zeros(DFTsize,size(x,2));

    % Preallocating the output blocks and the overlap for the first block:  
    y = zeros(size(Blocks,1)*size(Blocks,3),size(Blocks,2));
end
if P.GPU
    y      = gpuArray(y);
    yBFreq = gpuArray(yBFreq);
end

% Equalisation:
for j = 1:size(Blocks,3)
    % Input block with overlap:
    xB = [Overlap ; Blocks(:,:,j)];

    % Overlap:
    Overlap = xB(end-NOverlap+1:end,:);

    % FFT of the input block:
    xBFreq = (1/sqrt(DFTsize))*fftshift(fft(ifftshift(xB,1),DFTsize,1),1);

    % Equalisation:
    switch EqType
        case 'MIMO'
            % MIMO equalisation:
            for i = 1:size(x,2)
                yBFreq(:,i) = sum(reshape(W(:,i,:),[],size(x,2)).*xBFreq,2);
            end
        case 'MISO'
            % MISO equalisation:
            yBFreq(:,1) = sum(reshape(W(:,1,:),[],size(x,2)).*xBFreq,2);
        case 'SISO'
            % SISO equalisation:
            for i = 1:size(x,2)
                yBFreq(:,i) = W(:,i,i).*xBFreq(:,i);
            end
    end

    % IFFT of the block after equalisation (considering downsampling):
    yFDE = sqrt(DFTsize)*fftshift(ifft(ifftshift(yBFreq,1),DFTsize,1),1);

    % Output block:
    yB = yFDE(NOverlap/2+1:end-NOverlap/2,:);

    % Assigning the samples to the output signal:
    y((j-1)*(DFTsize-NOverlap)+1:j*(DFTsize-NOverlap),:) = yB;
end

% Quantity of samples to discard:
DInit = 1+(NOverlap)/2 ;  DFin = NExtra-NOverlap/2;

% Removing the overlapped zeros:
y = y(DInit:end-DFin,:);
y = bsxfun(@rdivide,y,sqrt(mean(abs(y).^2)));

end