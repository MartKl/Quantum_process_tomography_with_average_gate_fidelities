% Plot the data
% dataCliff array is indexed as follows: dataCliff(samples,idx,numQ,etaPower)
% This stores the Frobenius norm of the error of the reconstructed Choi
% matrix

% samples = # of measurement settings (also called m)
% idx = points to a specific run of the numerical simulation
% numQ = number of qubits
% etaPower = -log_{10}(noise strength)
clear

% load('run0x')
% load('./Cliff_m01/dataBerlin.mat');
load('./Haar_m01_02.mat');
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

% COMPUTE AVERAGES AND STANDARD DEVIATIONS OF THE ERRORS:

% set non-converged points to nan:
% dataCliff(isnan(optval)) = nan;


% Truncate error values that are > 1
dataCliffTrunc = min(dataCliff, 1);
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


% Calculate upper and lower bounds on Delta, for the error bars
Upper = Delta + Sigma;
Lower = max(Delta - Sigma, 10^(-15) );

% FIRST SET OF PLOTS:

% Fix the number of qubits
% How does Delta scale with the number of measurement settings, and the noise strength?
% Make a log-log plot

% Note, Matlab's loglog() function cannot draw error bars, so instead we
% rescale the data by hand and then use the errorbar() function to plot it.
% Note, I think a log2 plot is easier to read than log10, since within the 
% range of our data, there are more integer powers of 2 than there are 
% integer powers of 10.
%for numQ = 3

fig1 = figure(1);
% errorbar() draws one line for each column of the data matrix
EBX = kron([1,1,1], log10(samples)');  % X coordinates
EBY = log10(Delta);  % Y coordinates
EBL = EBY - log10(pos(Lower));  % Deviation below EBY
EBU = log10(Upper) - EBY;  % Deviation above EBY
errorbar(EBX, EBY, EBL, EBU)
xlim([2.15 3.35]);
ylim([-4,0.2]);
% plot(log2(samples), log2(Delta(:,numQ,1)), '-o', ...
%      log2(samples), log2(Delta(:,numQ,2)), '-o', ...
%      log2(samples), log2(Delta(:,numQ,3)), '-o');
FigTitle = ['\Delta(m) averaged over ', num2str(reps), ' realizations: '...
     , num2str(numQ), ', qubits, Haar random measurements'];
title(FigTitle);
xlabel(texlabel('log_10(m)'))
ylabel(texlabel('log_10(Delta)'))
legend([texlabel('eta = '), num2str(eta_list(1))], ...
       [texlabel('eta = '), num2str(eta_list(2))], ...
       [texlabel('eta = '), num2str(eta_list(3))], ...
       'Location', 'Best')
name1 = 'Delta(m)_3qubits_Haar_SDPT3';
saveas(fig1, [name1, '.png']);
savefig(fig1, name1);

% % SECOND PLOT: non-log
ebx = kron([1,1,1], samples');
fig2 = figure(2);
errorbar(ebx, Delta, Sigma);
xlim([100 2100]);
% ylim([0,0.1]);

title(FigTitle);
xlabel(texlabel('m'))
ylabel(texlabel('Delta'))
legend([texlabel('eta = '), num2str(eta_list(1))], ...
       [texlabel('eta = '), num2str(eta_list(2))], ...
       [texlabel('eta = '), num2str(eta_list(3))], ...
       'Location', 'Best')

name2 = 'Delta(m)_3qubits_Haar_SDPT3_non-log';
saveas(fig2, [name2, '.png']);
savefig(fig2, name2);

