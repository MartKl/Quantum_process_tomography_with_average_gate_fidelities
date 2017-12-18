% Concatenate different data sets
% dataCliff array is indexed as follows: dataCliff(samples,idx,numQ,etaPower)
% Each array entry is the Frobenius norm of the error of the reconstructed
% Choi matrix

clear

run1 = load('Haar_m01/dataBerlin.mat');
run2 = load('Haar_m02/dataBerlin.mat');
% run3 =

rps = 50;
data = struct;
% dims: (m,reps,numQ=3,eta)
data.frob_norm = run1.data.frob_norm(:,1:rps,:,:);
data.frob_norm(:,rps+1:2*rps,:,:) = run2.data.frob_norm(:,1:rps,:,:);
reps = 2*rps;

data.cvx_optval = run1.data.cvx_optval(:,1:rps,:,:);
data.cvx_optval(:,rps+1:2*rps,:,:) = run2.data.cvx_optval(:,1:rps,:,:);

data.cvxReps = run1.data.cvxReps(:,1:rps,:,:);
data.cvxReps(:,rps+1:2*rps,:,:) = run1.data.cvxReps(:,1:rps,:,:);

display(['Nr. failed reconstructions: ',...
 num2str( sum(isnan(data.cvx_optval(:))) ),...
 ' out of ',...
 num2str( length(data.cvx_optval(:)) )])

m_list = run1.m_list;
eta_list = run1.eta_list;
save('Haar_m01_02')
