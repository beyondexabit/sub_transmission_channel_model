function [k]=pauli_matrix(m)

k = m(1) * [0 1; 1 0] + m(2) * [0 -1i; 1i 0] + m(3) * [1 0; 0 -1];

end