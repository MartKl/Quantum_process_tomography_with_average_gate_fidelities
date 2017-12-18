function [ tau_matrix ] = tau_transform( C,h,d,input_tau )
%tau_transform transforms an inputs tau vector into an output tau matrix,
%transformed by the Clifford described by C, h, and d.
%   Detailed explanation goes here


numqubits=length(h)/2;

tau_matrix=eye(2^numqubits,2^numqubits);

for qubit_idx=1:numqubits
    if input_tau(qubit_idx)==0
        if input_tau(qubit_idx+numqubits)==1
            tau_matrix=tau_matrix*vec2tau(C(:,qubit_idx+numqubits))*(1i)^d(qubit_idx+numqubits)*(-1)^h(qubit_idx+numqubits);
        end
    else
        if input_tau(qubit_idx+numqubits)==0
            tau_matrix=tau_matrix*vec2tau(C(:,qubit_idx))*(1i)^d(qubit_idx)*(-1)^h(qubit_idx);
        else
            tau_matrix=tau_matrix*vec2tau(C(:,qubit_idx))*(1i)^d(qubit_idx)*(-1)^h(qubit_idx)...
                *vec2tau(C(:,qubit_idx+numqubits))*(1i)^d(qubit_idx+numqubits)*(-1)^h(qubit_idx+numqubits);
        end
    end
end




end

