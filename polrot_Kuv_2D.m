
function [kuvX_modified,rot_matrix_nModes_multiplicationfactor]=polrot_Kuv_2D(kuvX_2D,nModes,theta_modes,polrotmethod)

if nargin<4
polrotmethod='add';
end

v = randn(1,3);
v = v./sqrt(v*v');
theta_rot=rand(1)*pi;
rot_matrix_2x2 = pauli_matrix(theta_rot*v);

switch theta_modes
    case 'zero'
        theta_rot=0;%;
        v=[1 0 0];
        rot_matrix_2x2 = pauli_matrix(theta_rot*v);

        rot_matrix_nModes = repmat({rot_matrix_2x2}, 1, nModes);
        polrot= blkdiag(rot_matrix_nModes{:});
        rot_matrix_nModes_multiplicationfactor=[];
    case 'same'

        rot_matrix_nModes = repmat({rot_matrix_2x2}, 1, nModes);
        polrot= blkdiag(rot_matrix_nModes{:});
        rot_matrix_nModes_multiplicationfactor2=polrot;
        rot_matrix_nModes_multiplicationfactor2(rot_matrix_nModes_multiplicationfactor2==0)=1;
rot_matrix_nModes_multiplicationfactor=repmat(rot_matrix_2x2, nModes, nModes);
    case 'different'
        
        rot_matrix_nModes_multiplicationfactor=zeros(2*nModes,2*nModes);
        for jj=1:nModes
        for ii=1:nModes
            theta_rot=rand(1)*pi;
            v = randn(1,3);
            v = v./sqrt(v*v');
            rot_matrix_2x2 = pauli_matrix(theta_rot*v);
            if ii==1
                rot_matrix_nModes= rot_matrix_2x2;
            else
            rot_matrix_nModes = [rot_matrix_nModes,{rot_matrix_2x2}];
            end
rot_matrix_nModes_multiplicationfactor(2*(ii-1)+1:2*(ii-1)+2,2*(jj-1)+1:2*(jj-1)+2)=rot_matrix_2x2;
        end
        end
        polrot= blkdiag(rot_matrix_nModes{:});
         rot_matrix_nModes_multiplicationfactor2=polrot;
        rot_matrix_nModes_multiplicationfactor2(rot_matrix_nModes_multiplicationfactor2==0)=1;

end

switch polrotmethod
    case 'add' % addition of 2x2 matrix in diagonal
kuvX_modified=kuvX_2D + polrot;

    case 'mul_diag' % dot product with 2x2 in diagonal, keeping other elements 1
kuvX_modified=kuvX_2D.*rot_matrix_nModes_multiplicationfactor2;
 
    case 'mul' % dot product with 2x2 in all modes for all modes where theta can be same or different
kuvX_modified=kuvX_2D.*rot_matrix_nModes_multiplicationfactor;

    case 'add_oldkuv'
      kuvX_2D(1:2:end,2:2:end)=0;
            kuvX_2D(2:2:end,1:2:end)=0;
       kuvX_modified=kuvX_2D + polrot;     

end