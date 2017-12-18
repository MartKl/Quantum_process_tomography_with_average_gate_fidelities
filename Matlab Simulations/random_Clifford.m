%Easwar Magesan June 2011

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Produces a uniformly random Clifford element
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P,C,h,d]=random_Clifford(n)%Takes in the number of qubits and outputs the symplectic matrix 'C' and vectors 'd' and 'h' which uniquely determine a Clifford element up to phase


U=zeros(2*n,2*n); %define matrix 'U'
U(1:n,1:n)=zeros(n);
U(1:n,n+1:2*n)=eye(n);
U(n+1:2*n,1:n)=zeros(n);
U(n+1:2*n,n+1:2*n)=zeros(n);
P=U+U'; %define matrix with which symplecticity is defined with respect to (note 'P' switches the first n entries and last n entries of a 2*n length vector)
%disp(P)

C=zeros(2*n,2*n); %symplectic matrix for the Clifford operation
h=zeros(2*n,1); %vector of negative signs for the Clifford operation
d=zeros(2*n,1); %vector of phases for the Clifford operation

while C(:,1) == zeros(2*n,1) %ensure the random vector chosen for the first column of 'C' is non-zero
    C(:,1)=randi(2,2*n,1)-1; % define first column of 'C' to be a random bit string
end

for j=2:1:2*n %'j' corresponds to which of the first n columns of 'C' we are currently trying to find
    kit = 0; %label that keeps track of whether two randomly chosen columns of 'C' are linearly dependent (which can't happen)
    while kit == 0
        kit=1; %flip 'kit' initially so that if two columns of 'C' end up being the same the 'while' loop continues to find another random solution
        M=zeros(j-1,2*n); %if we are trying to find the j'th column of 'C' then we have a system of 2*n variables in j-1 equations
        b=zeros(j-1,1); %define the non-homogeneous vector for the system
        if j >= n+1
            b(j-n,1)=1; %in finding the first n columns of 'C', the non-homogeneous vector in the system is the zero vector and in finding the last n columns of 'C', the non-homogeneous vector in the
            %system is equal to '1' at j-n'th entry (due to commutation
            %relations)
        end
        for k=1:1:j-1
            M(k,:) = (P*C(:,k))'; %define the coefficient matrix row by row
        end
%       [RREF,randomsoln]=Gaussianelimold(M,b); %find the random non-zero solution satisfying the current system of equations
%       C(:,j)=randomsoln'; %define the j'th column of 'C' to be the random
%       %non-zero solution found...note the last column of RREF contains the
%       %reduced nonhomogeneous term
        [M_reduced,b_reduced]=gaussianelim(M,b); %Perform Gaussian elimination on the system of equations M*x=b
        randomsoln=randomsolution(M_reduced,b_reduced);
        C(:,j)=randomsoln'; %define the j'th column of 'C' to be the random non-zero solution found
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Test to make sure the new column of 'C' does not form a linearly
        % dependent set with the previously chosen columns
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        R=zeros(j,2*n); %define matrix that represents the new system of equations that includes the new j'th column of 'C'
        
        for k=1:1:j
            R(k,:) = (P*C(:,k))'; %construct the matrix row by row
        end
        
        l=zeros(j,1); %construct a nonhomogeneous vector to input into the function 'Gaussianelim'
        [R_reduced,l_reduced]=gaussianelim(R,l); %Obtain row reduced version of 'R'
        
        if R_reduced(j,:)==zeros(1,2*n) %if the RREF of R contains a row of zero's we keep 'kit' equal to 0 and find another random solution for the j'th column of 'C' (ie. one that forms a 
            %linearly independent set with the previously chosen j-1 columns)
            kit=0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end


C; %symplectic matrix for the random Clifford
h=randi(2,2*n,1)-1; %vector of negative signs for the random Clifford
h; %vector of negative signs
for j=1:1:2*n
    d(j,1)=mod(C(:,j)'*U*C(:,j),2); %define 'd' (this is deterministically chosen from 'C' since the columns of 'C' represent
 %Hermitian elements of the Pauli group)
end
d;