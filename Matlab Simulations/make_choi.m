function rho = make_choi( U, state )
% MAKE_CHOI Creates a Choi representation matrix of the unitary U. 
%    state is cell of state vectors of maximally entangled states of several different
%			dimensions
%
%   author: Shelby Kimmel
%
%   Distributed under Creative Commons Attribution 4.0 International (CC BY 4.0) 
%            See: https://creativecommons.org/licenses/by/4.0/
%
%   Created: 2017

dim=size(U,1);
Ukron=kron(U,eye(dim));
newstate=Ukron*state{dim};
rho=kron(newstate,newstate');



end

