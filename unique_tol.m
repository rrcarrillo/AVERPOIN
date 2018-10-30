% unique_A=unique_tol(A) remove the duplicate column vectors from A having
% a level of tolerance when comparing numbers
function unique_A=unique_tol(A)
num_tol=eps(100); % Numeric tolerance
nvectors=size(A,2); % Number of input vectors
unique_A=[];
if numel(A) > 0 % The algorithm only works if the input matrix is empty
   for nvector=1:(nvectors-1) % For each column vector expect the last one
      current_vec=A(:,nvector);
      % Check if the current column vector is equal to any of the remainding
      % ones: A(:,(nvector+1):end)
      if ~any(all(abs(bsxfun(@minus,A(:,(nvector+1):end),current_vec)) < num_tol,1))
         unique_A=[unique_A current_vec]; % Current vector is not duplicated
      end
   end
   unique_A=[unique_A A(:,end)]; % The last column vector is always included
end