%Easwar Magesan June 2011

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Performs Gaussian elimination for $y x z$ (y \leq z) matrices over the
%finite field $Z_2$ with an extra column for non-homogeneous term (ie. puts into RREF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [RREF,nonhom_reduced]=Gaussianelim(coeff,nonhom)%Takes in a system of equations with coefficient matrix 'coeff' and (possibly) non-homogeneous term 'nonhom'. Outputs the reduced row echelon 
%form of the system (ie. the matrix resulting from performing row
%operations on the system to bring 'coeff' into RREF) as well as the transformed column 'nonhom_reduced')
dims=size(coeff);
y=dims(1,1); %number of rows in system of equations
z=dims(1,2); %number of variables in the system of equations
A=zeros(y,z+1); %Define a matrix that represents the total system of equations (there is an extra column for the non-homogeneous term)
A(1:y,1:z)=coeff; %Construct A
A(:,z+1)=nonhom;
isler = 0; %variable representing how many '1's have been moved to the standard position in REF in the Gaussian elimination process (ie. how many non-zero columns have been encountered)
swit=0;

for j=1:1:z %for each column of coefficient matrix
    for k=isler+1:1:y %for each row with index greater than the largest row index that has been put in REF in Gaussian elimination process
        if A(k,j)~=0 %find the first non-zero value in the j'th column
            swit=1; %set marker equal to '1' if there is a non-zero value in j'th column
            B=A; %define a matrix for swapping the k'th and (isler + 1)'th rows
            A(isler+1,:)=B(k,:); %swap the k'th and (isler + 1)'th rows
            A(k,:)=B(isler+1,:);
        break %break the 'for' loop at this first non-zero entry
        end
    end
    
    if swit==1 %if there was a non-zero value in j'th column then put column in standard form by zeroing all rows other than (isler + 1)'th row
        for l=1:1:isler %still looking at column j, for every row with index \leq isler that has a '1', perform row operation to obtain a '0'
            if A(l,j)==1
                A(l,:)=A(l,:)+A(isler+1,:);
                    for q=j:1:z+1 %set the values that get added to '2' in the row addition equal to '0' (ie. mod 2 addition)
                        if A(l,q)==2
                            A(l,q)=0;
                        end
                    end
            end
        end
        for l=isler+2:1:y %still looking at column j, for every row with index greater than (isler + 1) that has a '1', perform row operation to obtain a '0'
            if A(l,j)==1
                A(l,:)=A(l,:)+A(isler+1,:);
                    for q=j:1:z+1 %set the values that get added to '2' in the row addition equal to '0' (ie. mod 2 addition)
                        if A(l,q)==2
                            A(l,q)=0;
                        end
                    end
            end
        end
        isler = isler + 1; %since there was a non-zero value in j'th column update isler to isler+1
    end
swit=0; %reset the marker 'swit'  
    
end

RREF=A(1:y,1:z);
nonhom_reduced=A(:,z+1);