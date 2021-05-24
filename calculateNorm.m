function Norm = calculateNorm(vec,dim)
%This function gets a 2D vector or an array of 2D vectors and returns the 
%norm of these vectors. vec can be 2xM or Nx2.

%In addition the function can get a 3D matrix of Row x Col x 2 and
%calculate the norm of each vector in the Z dimension
if length(size(vec))==2 %if the matrix is 2D
    if size(vec,1)==2 && size(vec,2) ~= 2
        vec = vec';
        Norm = sqrt(sum(vec.^2,2));
        Norm = Norm';
    elseif size(vec,2)==2
        Norm = sqrt(sum(vec.^2,2));
    else
        %     fprintf('Error Vector is not a 2D vector\n');
        Norm = [];
    end
    
elseif length(size(vec))==3 %if the matrix is 3D
    Norm = sqrt(sum(vec.^2,3));
else
    disp(['operating on dimension: ',num2str(dim)]);
    Norm = sqrt(sum(vec.^2,dim)); 
%     % if neither is size 2 - searach all dims and operate if you find one
%     % that is size 2
%     dims = size(vec);
%     relD = find(dims==2);
%     
%     if ~isempty(relD)
%         disp(['operating on dimension: ',num2str(relD)]);
%         Norm = sqrt(sum(vec.^2,relD));
%     else
%         disp('none of the matrix dimension size 2 - retrun empty')
%         Norm = [];
%     end
end

