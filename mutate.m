function parent_chromosome = mutate(parent_chromosome,rate,mum)
    global nVar numOfObj VarMax VarMin CostFunction;
    [N,~]=size(parent_chromosome);
    index = randperm(N,floor(rate*N));
    for i = index
        % Select at random the parent.
            % Get the chromosome information for the randomnly selected parent.
            child_3 = parent_chromosome(i,:);
            % Perform mutation on eact element of the selected parent.
            for j = 1 :nVar
               r(j) = rand(1);
               if r(j) < 0.5
                   delta(j) = (2*r(j))^(1/(mum+1)) - 1;
               else
                   delta(j) = 1 - (2*(1 - r(j)))^(1/(mum+1));
               end
               % Generate the corresponding child element.
               child_3(j) = child_3(j) + delta(j);
               % Make sure that the generated element is within the decision
               % space.
               if child_3(j) > VarMax(j)
                   child_3(j) = VarMax(j);
               elseif child_3(j) < VarMin(j)
                   child_3(j) = VarMin(j);
               end
            end
            % Evaluate the objective function for the offspring and as before
            % concatenate the offspring chromosome with objective value.    
            child_3(:,nVar + 1: numOfObj + nVar) = evaluate_objective(child_3, numOfObj, nVar);
            parent_chromosome(i,:) = child_3(:);
            % Set the mutation flag
    end
    for i = 1:N
        parent_chromosome(i,nVar + 1: numOfObj+nVar) = CostFunction(parent_chromosome(i,1:nVar));
    end
    clear index i;
end
