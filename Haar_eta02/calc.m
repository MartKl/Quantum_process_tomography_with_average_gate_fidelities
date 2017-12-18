%rng('shuffle')
% run C:\cvx\cvx_setup.m
% run /home/kliesch/.cvx/cvx_setup.m

load('dataBerlin.mat');

p0=genpath('../');
p1=genpath('../QETLAB-master');
p2=genpath('../Matlab Simulations');
addpath(p0,p1,p2,p0);

t0=clock;
while j<=anz && etime(clock,t0)<T
    im = jobs(j,1);
    m = m_list(im); %number measurements
    ieta = jobs(j,4);
    eta = eta_list(ieta); %mesurement noise
    inumQ = jobs(j,3);
    numQ = numQ_list(inumQ); %number qubits
    k = jobs(j,2); %rep
    
    display(['Job ' job_text ...
             ' (m, eta, rep) = (' ...
                num2str(m) ', ' ...
                num2str(eta) ', '...
                num2str(k) '), '...
             datestr(now)]);
    % ---------------------------------------
    % job number as seed 
    rng(j+offset);
    % ---------------------------------------
    output=HaarCompressedSensingNew(numQ, eta , m, state); % Note: numQ, not 2^numQ
    
    data.frob_norm(im, k, inumQ, ieta)      = output.frob_norm;
    data.cvx_cputime(im, k, inumQ, ieta)    = output.cvx_cputime;
    data.cvx_optbnd(im, k, inumQ, ieta)     = output.cvx_optbnd;
    data.cvx_optval(im, k, inumQ, ieta)     = output.cvx_optval;
    data.cvx_slvitr(im, k, inumQ, ieta)     = output.cvx_slvitr;
    data.cvx_slvtol(im, k, inumQ, ieta)     = output.cvx_slvtol;
    data.cvxReps(im, k, inumQ, ieta)     = output.slack+1;
    
    j=j+1;
    save('dataBerlin.mat');
    int2file(j,'j.txt');
end

display(['...exiting...']);
exit