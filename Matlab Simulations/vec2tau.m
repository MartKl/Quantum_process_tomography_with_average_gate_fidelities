function [tau,tau_vec] = vec2tau( vec )
% VEC2TAU  converts a symplectic vector representation of a paul to a tau representation of a pauli
%
% For info about tau vector and symplectic vector, see: https://arxiv.org/pdf/quant-ph/0304125.pdf
%   
%   
%
%   author: Shelby Kimmel
%
%   Distributed under Creative Commons Attribution 4.0 International (CC BY 4.0) 
%            See: https://creativecommons.org/licenses/by/4.0/
%
%   Created: 2017


d=length(vec)/2;
tau_vec=zeros(d,1);

for qubit_idx=1:d
    if vec(qubit_idx)==0
        if vec(qubit_idx+d)==0
            tau_vec(qubit_idx)=0;
        else
            tau_vec(qubit_idx)=1;
        end
    else
        if vec(qubit_idx+d)==0
            tau_vec(qubit_idx)=3;
        else
            tau_vec(qubit_idx)=2;
        end
    end
end

tau=tau_creator(tau_vec,0);



end

