function Westim = EqualiserCalc(Hres_estimated,SNR,P)
% V*H*U:
% MMSE equaliser:
if P.GPU
    Westim = zeros(size(Hres_estimated),'gpuArray');
else
    Westim = zeros(size(Hres_estimated));
end
for i = 1:size(Hres_estimated,3)
    % MIMO matrix (P.NTrib x P.NTrib) for the i-th subcarrier:
    Hi = Hres_estimated(:,:,i);

    % MMSE equaliser computation:
    Vi  = eye(P.NTrib)*1/SNR + Hi'*Hi;
    Wi  = Vi\Hi';

    % Assigning the coeffs. to a (P.NDsc + P.NPsc) x P.NTrib x P.NTrib
    % matrix:
    Westim(:,:,i)  = Wi;
end
end