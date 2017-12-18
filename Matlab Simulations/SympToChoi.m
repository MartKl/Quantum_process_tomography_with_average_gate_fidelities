function [ choi_mat ] = SympToChoi( P,C,h,d)
% SYMPTOCHOI Changes a Sympletic representation of a Clifford to a Choi representation
%   P,C,h,d are the output from Magesan's "random_Clifford." 
%   The columns of C correspond to (in order) how the Paulis
%   Z_1,Z_2,...,Z_2,X_1,X_2,...,X_n are transformed by the Clifford
%
%
%   author: Shelby Kimmel
%
%   Distributed under Creative Commons Attribution 4.0 International (CC BY 4.0) 
%            See: https://creativecommons.org/licenses/by/4.0/
%
%   Created: 2017


numqubits=size(C,1)/2;
choi_mat=zeros(4^numqubits);

%iterate over all Paulis
for idx1=0:2^numqubits-1
    for idx2=0:2^numqubits-1
        xs=de2bi(idx1,numqubits);
        zs=de2bi(idx2,numqubits);
        vec=[xs,zs]; %symplectic representation of tau vector
        numys=sum(and(xs,zs)); %number of y's in the tensor product
        origpauli=vec2tau(vec)*(-1i)^numys; %convert to pauli by first changing to tau, then adding i's
        transformed_pauli=tau_transform(C,h,d,vec)*(-1i)^numys;
        choi_mat=choi_mat+kron(transformed_pauli,origpauli)/4^numqubits*(-1)^numys;
    end
end

end

