% Plot the data
% dataCliff array is indexed as follows: dataCliff(samples,idx,numQ,etaPower)
% This stores the Frobenius norm of the error of the reconstructed Choi
% matrix

% samples = # of measurement settings (also called m)
% idx = points to a specific run of the numerical simulation
% numQ = number of qubits
% etaPower = -log_{10}(noise strength)

clear 
% load('run1x')
load('./Haar_eta01_02_03.mat');
% load('./Haar_eta01/dataBerlin.mat');
% reps = jobs(j,2)-1;
dataCliff = squeeze( data.frob_norm(:,1:reps,:,:) );
optval = squeeze( data.cvx_optval(:,1:reps,:,:) );
cvxReps = squeeze( data.cvxReps(:,1:reps,:,:) );


%check precision decreas and failures
n_fails = squeeze(sum(isnan(optval),2));
display([ num2str(sum(n_fails(:))), ' fails out of ', num2str(length(optval(:))) ]);

n_runs = length( cvxReps(:));
n_reps = sum(cvxReps(:));
display([num2str( n_reps ) ' runs to solve ' num2str( n_runs ) ' SDPs, i.e., '... 
    'machine precision ' num2str(n_reps-n_runs) ' times increased']);


% Truncate error values that are > 2
% dims: m, reps, numQ, eta
dataCliffTrunc = min(dataCliff, 2);
numQ=3;

% For j = 1:length(samples), the j'th data point will correspond to the 
% simulation where the number of samples was samples(j)
% Exclude those simulation runs that have fewer than 201 samples
samples = m_list;

% Calculate the averages and standard deviations of the errors
% Call these Delta and Sigma
% Along the way, we calculate the 2nd moments of the errors -- call this Tmp2
% They are indexed like this: Delta(j, numQ, etaPower)
num_trials = size( dataCliffTrunc, 2);
% AVERAGE
Delta = squeeze( mean( dataCliffTrunc, 2, 'omitnan') );
% STANDARD DEVIATION
Sigma = squeeze( std(dataCliffTrunc, 0, 2, 'omitnan') );  

n_fails = squeeze(sum(isnan(optval),2));
display([ num2str(sum(n_fails(:))), ' fails out of ', num2str(length(optval(:))) ]);

% Calculate upper and lower bounds on Delta, for the error bars
Upper = Delta + Sigma;
Lower = max(Delta - Sigma, 10^-15);


% % THIRD SET OF PLOTS:
% % 
% % Fix the number of qubits
% % For each choice of the number of measurement settings m, draw one line:
% % How does the reconstruction error \Delta depend on the noise strength \eta?
n_eta = size(Delta,2);
n_m = size(Delta,1);


numQ = 3;
plt1 = figure(1);
% Draw one line for each j = 1:length(samples)
% errorbar() draws one line for each column of the data matrix, 
EBX = kron(ones([1,n_m]), log10(eta_list)');  % X coordinates
EBY = log10(Delta');         % Y coordinates
EBL = EBY - log10(Lower');   % Deviation below EBY: length of error bar
EBU = log10(Upper') - EBY;   % Deviation above EBY
errorbar(EBX, EBY, EBL, EBU);
xlim([-4.2 0]);
ylim([-4.2,.2]);

FigTitle = ['\Delta(\eta) averaged over ', num2str(reps), ' realizations: '...
     , num2str(numQ), ', qubits, Haar random measurements'];
title(FigTitle);
xlabel(texlabel('log_10(eta)'))
ylabel(texlabel('log_10(Delta)'))
legend([texlabel('m = '), num2str(m_list(1))], ...
       [texlabel('m = '), num2str(m_list(2))], ...
       [texlabel('m = '), num2str(m_list(3))], ...
       [texlabel('m = '), num2str(m_list(4))], ...
       [texlabel('m = '), num2str(m_list(5))], ...
       'Location', 'Best')
name1 = 'Delta(eta)_3qubits_Haar_SDPT3';

saveas(plt1, [name1, '.png']);
savefig(plt1, name1);



plt2 = figure(2);
ebx = kron(ones([1,n_m]), eta_list');
eby = Delta';
eb = Sigma';
errorbar(ebx,eby,eb);
xlim([-.05 1.05]);
ylim([-.05,1.1]);

title(FigTitle)
xlabel(texlabel('eta'))
ylabel(texlabel('Delta'))
legend([texlabel('m = '), num2str(m_list(1))], ...
       [texlabel('m = '), num2str(m_list(2))], ...
       [texlabel('m = '), num2str(m_list(3))], ...
       [texlabel('m = '), num2str(m_list(4))], ...
       [texlabel('m = '), num2str(m_list(5))], ...
       'Location', 'Best');
   
name2 = 'Delta(eta)_3qubits_Haar_SDPT3_non-log';
savefig(plt1, name2);
saveas(plt2, [name2, '.png']);


