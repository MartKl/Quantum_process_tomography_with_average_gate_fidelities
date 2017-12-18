function [ frob_norm ] = CliffordCompressedSensingNew( numqubits, eta, samples,state)
% CliffordCompressedSensing: runs compressed sensing based on RBT outcomes
%
%   dim: dimension of the Hilbert space that the operation acts
%   eta: error tolerance, eta in our paper
%   samples: number of other Haar random unitaries that are used in the measurement vector
%   state: is a cell that contains maximally entangled states to make making choi matrices easier
% 
%   author: Shelby Kimmel
%
%   Distributed under Creative Commons Attribution 4.0 International (CC BY 4.0) 
%            See: https://creativecommons.org/licenses/by/4.0/
%
%   Created: 2017


% Import packages from +qip library
import qip.random.unitary
import qip.open_systems.*


% Create a random unitary
Utarget=unitary(2^numqubits);
%UtargetChoi=liou2choi(liou(Utarget,Utarget'));
UtargetChoi=make_choi(Utarget,state);
UtargetChoiVec=reshape(UtargetChoi,[16^numqubits,1]);


% Create list of other unitaries to compare to
unitary_list=zeros(samples,16^numqubits);
outcomes=zeros(samples,1);
for idx=1:samples
    [P,C,h,d]=random_Clifford(numqubits);
    utemp_choi=SympToChoi(P,C,h,d);
    %utemp_choi=liou2choi(liou(utemp',utemp));
    %utemp_choi=make_choi(utemp,state);
    unitary_list(idx,:)=reshape(utemp_choi,[1,16^numqubits]);
    outcomes(idx)=unitary_list(idx,:)*UtargetChoiVec; 
end

%Add error to data:
unscaled_error_list=normrnd(0,1,samples,1);
error_list=unscaled_error_list/(norm(unscaled_error_list))*eta;
outcomes=outcomes+error_list;

cvx_begin sdp
    variable rho(4^numqubits,4^numqubits) hermitian
    minimize((unitary_list*reshape(rho,[16^numqubits,1])-outcomes)'*...
        (unitary_list*reshape(rho,[16^numqubits,1])-outcomes));
    TrX(rho,1,[2^numqubits,2^numqubits]) == trace(rho)*eye(2^numqubits)/2^numqubits; 
    TrX(rho,2,[2^numqubits,2^numqubits]) == trace(rho)*eye(2^numqubits)/2^numqubits;
    rho == hermitian_semidefinite(4^numqubits);
cvx_end
frob_norm=norm(rho-UtargetChoi,'fro');


end

