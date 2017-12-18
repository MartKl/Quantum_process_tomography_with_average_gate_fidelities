% Concatenate different data sets
% dataCliff array is indexed as follows: dataCliff(samples,idx,numQ,etaPower)
% Each array entry is the Frobenius norm of the error of the reconstructed
% Choi matrix

clear

run1 = load('Haar_eta01/dataBerlin.mat');
run2 = load('Haar_eta02/dataBerlin.mat');
run3 = load('Haar_eta02/dataBerlin.mat');
% run3 =

% rps = 50;
data = struct;
% dims: (m,reps,numQ=3,eta)
data.frob_norm = run1.data.frob_norm(:,1:35,:,:);
data.frob_norm(:,35+1:2*35,:,:) = run2.data.frob_norm(:,1:35,:,:);
data.frob_norm(:,2*35+1:2*35+30,:,:) = run3.data.frob_norm(:,1:30,:,:);
reps = 2*35+30;

data.cvx_optval = run1.data.cvx_optval(:,1:35,:,:);
data.cvx_optval(:,35+1:2*35,:,:) = run2.data.cvx_optval(:,1:35,:,:);
data.cvx_optval(:,2*35+1:2*35+30,:,:) = run3.data.cvx_optval(:,1:30,:,:);

data.cvxReps = run1.data.cvxReps(:,1:35,:,:);
data.cvxReps(:,35+1:2*35,:,:) = run2.data.cvxReps(:,1:35,:,:);
data.cvxReps(:,2*35+1:2*35+30,:,:) = run3.data.cvxReps(:,1:30,:,:);

display(['Nr. failed reconstructions: ',...
 num2str( sum(isnan(data.cvx_optval(:))) ),...
 ' out of ',...
 num2str( length(data.cvx_optval(:)) )])

m_list = run1.m_list;
eta_list = run1.eta_list;
save('Haar_eta01_02_03')
