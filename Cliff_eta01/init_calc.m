% m: measurements
% numQ: number qubits
% eta: noise strength
clear

restoredefaultpath
% run 'C:\cvx\cvx_setup.m';
run /home/martin/.cvx/cvx_setup.m
cvx_solver sdpt3


% ------- parameter ---------
% number of qubits
numQ_list=[3];
% m_list: number of average gate fidelities/ samples
% eta_list: noise strengths

% Delta(m)
% m_list=[160:10:340 360:20:520 550 580 710 750:50:800 900 1000 1100 1200:200:1600 2000];
% eta_list = [.1, .01, .001]; %measurement noise

% Delta(eta)
m_list=[240, 280, 340, 600, 1800];
eta_list = logspace(0,-4,23); %measurement noise

%error bar check
% m_list=[250, 1000];
% eta_list = logspace(-3,-4,4); %measurement noise

%nr. realizations
reps = 100;
%
%max number of comptuation time before restarting matlab
T = 3*60*60;

% ------- test --------
% numQ_list=[1];
% m_list=20:20:200;
% eta_list = [.1, .01]; %measurement noise
% reps = 2;
% T = 20;

% ------- initialize --------
[~, folder] = fileparts(pwd);
job_text = folder;

rng_parameter = str2num( job_text(end-1:end) )-1;

offset = rng_parameter*reps; %to manage rng(j)

mL = length(m_list);
etaL = length(eta_list);
numQL = length(numQ_list);

anz= mL*reps*numQL*etaL;

[w, x, y, z]= ndgrid(1:etaL,1:mL,1:numQL,1:reps);
jobs = horzcat(reshape(w,[anz,1]),reshape(x,[anz,1]),reshape(y,[anz,1]),reshape(z,[anz,1]));
jobs = [jobs(:,2) jobs(:,4) jobs(:,3) jobs(:,1)]; % in order to run over reps at last

j=1;
int2file(j,'j.txt');
int2file(anz,'anz.txt');

% m, k, numQ, eta
data              = struct;
data.frob_norm    = zeros(mL,reps,numQL,etaL);
data.cvx_cputime  = zeros(mL,reps,numQL,etaL);
data.cvx_optbnd   = zeros(mL,reps,numQL,etaL);
data.cvx_optval   = zeros(mL,reps,numQL,etaL);
data.cvx_slvitr   = zeros(mL,reps,numQL,etaL);
data.cvx_slvtol   = zeros(mL,reps,numQL,etaL);
data.cvxReps   = zeros(mL,reps,numQL,etaL);

t0 = clock;

state = importdata('state.mat'); %mk: what is that for?
save('dataBerlin.mat');

exit