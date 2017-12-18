function output = HaarCompressedSensingNew( numqubits, eta, samples,state)
%HaarCompressedSensing: runs compressed sensing based on RBT outcomes
%   dim: dimension of the Hilbert space that the operation acts
%   eta: error tolerance, eta in our paper
%   samples: number of other Haar random unitaries that are used in the
%   measurement vector
%   state: is a cell that contains maximally entangled states to make making
%   choi matrices easier
% Upgrade: check precision in the end; if not: rerun with slightly relaxed
% equalequality constraints 
% change: measurement unitaries drawn from the Haar measure

% Import packages from Marcus' directory
import qip.random.unitary
import qip.open_systems.*
import qip.partial_trace


% Create a random unitary
Utarget=unitary(2^numqubits);
%UtargetChoi=liou2choi(liou(Utarget,Utarget'));
UtargetChoi=make_choi(Utarget,state);
UtargetChoiVec=reshape(UtargetChoi,[16^numqubits,1]);


% Create list of other unitaries to compare to
unitary_list=zeros(samples,16^numqubits);
outcomes=zeros(samples,1);
for idx=1:samples
%     [P,C,h,d]=random_Clifford(numqubits);
%     utemp_choi=SympToChoi(P,C,h,d);
    %utemp_choi=liou2choi(liou(utemp',utemp));
    %utemp_choi=make_choi(utemp,state);
    utemp_choi=make_choi( unitary(2^numqubits) ,state); % draw Haar random measurements
    unitary_list(idx,:)=reshape(utemp_choi,[1,16^numqubits]);
    outcomes(idx)=unitary_list(idx,:)*UtargetChoiVec; 
end

%Add error to data:
unscaled_error_list= randn(samples,1);
error_list=unscaled_error_list/(norm(unscaled_error_list))*eta;
outcomes=outcomes+error_list;

cvx_begin sdp quiet
    variable rho(4^numqubits,4^numqubits) hermitian semidefinite
    minimize((unitary_list*reshape(rho,[16^numqubits,1])-outcomes)'*...
        (unitary_list*reshape(rho,[16^numqubits,1])-outcomes));
    TrX(rho,1,[2^numqubits,2^numqubits]) == eye(2^numqubits)/2^numqubits; 
    TrX(rho,2,[2^numqubits,2^numqubits]) == eye(2^numqubits)/2^numqubits;
cvx_end

check = strcmp(cvx_status, 'Solved');
j = 0;
while not(check) && j<=6
    myeps = eps*10^j;
    display(['Not solved, set slack: ', num2str(myeps) ]);
    j = j+1;
    
    outcomes = outcomes + eps*randn(samples,1);
    cvx_begin sdp quiet
        variable rho(4^numqubits,4^numqubits) hermitian semidefinite
            minimize((unitary_list*reshape(rho,[16^numqubits,1])-outcomes)'*...
            (unitary_list*reshape(rho,[16^numqubits,1])-outcomes));
        norm( TrX(rho,1,[2^numqubits,2^numqubits]) ...
            -eye(2^numqubits)/2^numqubits, 'fro')<= myeps; 
        norm( TrX(rho,2,[2^numqubits,2^numqubits]) ...
            -eye(2^numqubits)/2^numqubits, 'fro')<= myeps;
    cvx_end
    
    check = strcmp(cvx_status, 'Solved');
end

% output structure containing all produced data
output              = struct;
output.frob_norm    = norm(rho-UtargetChoi,'fro');
output.cvx_cputime  = cvx_cputime;
output.cvx_optbnd   = cvx_optbnd;
output.cvx_optval   = cvx_optval;
output.cvx_slvitr   = cvx_slvitr;
output.cvx_slvtol   = cvx_slvtol;
output.slack            = j;
% output.rho          = rho;
% output.targetChoi   = UtargetChoi;

end

