%Easwar Magesan June 2011

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Finds random solution for $y x z$ (y \leq z) matrices A over the finite
%field $Z_2$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [randomsoln]=Randomsolution(RREFcoeff,Rnonhom)%Takes in a system of equations in RREF with coefficient matrix 'RREFcoeff' and non-homogeneous reduced term 'Rnonhom' and outputs a random solution to the system of
%equations. 



 dims=size(RREFcoeff);
 y=dims(1,1); %number of rows in system of equations
 z=dims(1,2); %number of variables in the system of equations
 
 counter = 0; %used for determining whether a solution exists
 
 for j=1:1:y
    if (RREFcoeff(j,:)==zeros(1,z)) & (Rnonhom(j,1)~=0) %if there is a row of zeros in the coefficient matrix with a non-zero value in the nonhomogeneous vector there is no solution and so increase the value of counter
        counter=counter + 1;
    end
 end
        
 if counter > 0
     randomsoln='no solution'; %if counter > 0 then there is no solution
 else %if counter = 0 a solution exists so proceed with finding a random solution

    RREFappend=zeros(y,z+1); %combine coefficient matrix with nonhomogeneous term into one matrix 'RREFappend'
    RREFappend(1:y,1:z)=RREFcoeff;
    RREFappend(:,z+1)=Rnonhom;


   

    holde = zeros(y,1); %variable for recording what column the first non-zero value in the k'th row occurs at (ie. records where the leading 1's are)

    for k=1:1:y %for each row
        for j=1:1:z % for each column not corresponding to the non-homogeneous term
            if RREFappend(k,j)~=0 %find the first non-zero value in the k'th row
                holde(k,1)=j; %records what column the first non-zero value in the k'th row occurs at
            break %break the 'for' loop over columns at this first non-zero entry
            end
        end
    end

    holde;

    temp1=y+1; %variable used for reshaping 'holde' to get rid of any '0' entries that may occur at the end of the vector (initially set it to 'y+1' in case there are no zero's in holde, ie. the coefficient matrix is invertible)

    for j=1:1:y %for each row of RREFappend
        if holde(j,1)==0 %find first '0' entry of holde
            temp1=j;  %set 'temp1' to be the index at which the '0' entry occurs
            break
        end
    end

    holde=holde(1:temp1-1); %redefine 'holde' to only include the non-zero entries
    size1=size(holde); %find the new dimensions of holde
    nonran=size1(1,1); %'nonran' is equal to the number of non-zero entries in holde, ie. the number of leading '1's in RREFappend, which will be the number of non-random components in finding the solution

    r=zeros(5,z); %the first row of 'r' labels the columns that contain the leading ones in RREFappend (ie, these column indices correspond to the entries of holde), the second row labels 
    %the row in RREFappend corresponding to each leading '1',the third row 
    %assigns a uniformly random bit to the columns that don't contain a leading
    %'1', the fourth row gives the solution for the variables with leading '1'
    %in A, and the fifth row contains the uniformly random solution (by just adding
    %the second row and third row together)

    counts=1; %counter for second row of 'r'

    for k=1:1:nonran
        r(1,holde(k,1))=1; %label the columns that contain the leading ones in RREFappend by a '1'
        r(2,holde(k,1))= counts; %label the row corresponding to the leading one in RREFappend
        counts = counts + 1; %increment counter
    end

    while r(5,:) == 0 %ensure that the random solution is non-zero (ie. keep repeating the algorithm to find random solution as long as r(5,:) is the zero vector)

        for j=1:1:z %assign a uniformly random bit to the columns of RREFappend that don't contain a leading
        %'1'
            if r(1,j)==0
                r(3,j)=randi(2,1,1)-1; 
            end
        end

        for j=1:1:z %$for all of the variables
            if r(1,j)==1 %look at only the columns of RREFappend with a leading 1
                for k=1:1:z %for all of the variables
                    if r(1,k)==0 %look at only the columns of RREFappend that don't contain a leading 1 (because the other columns won't contribute to the value of the j'th variable (that corresponds to a 
                        %leading 1))
                        r(4,j)=r(4,j)+r(3,k)*RREFappend(r(2,j),k); %calculate the value of the j'th variable (that corresponds to a leading 1) by taking dot product of the entries of this row in RREFappend that 
                        %aren't leading 1's with the random bits previously
                        %assigned
                        r(4,j)=mod(r(4,j),2); %represent the value by its modular value
                    end
                end
                r(4,j)=RREFappend(r(2,j),z+1)+r(4,j); , %subtract (mod 2) the dot product from the z+1'th column of RREFappend
                r(4,j)=mod(r(4,j),2); %calculate the value of the variables that have leading 1's modulo two
            end
        end
        
        r(5,:)=r(3,:)+r(4,:); %Add the disjoint third and fourth rows of r, each of which contains the random values for the variables corresponding to columns that have (ie. row 4) and don't have (ie. row 3) leading 1's
        for j=1:1:z
            r(5,j)=mod(r(5,j),2); %calculate the value of all variables modulo two
        end
    
    end

    randomsoln=r(5,:);

end