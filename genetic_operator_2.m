function f  = genetic_operator_2(parent_chromosome, M, V, mu, mum, l_limit, u_limit, varargin)

%% function f  = genetic_operator(parent_chromosome, M, V, mu, mum, l_limit, u_limit)
% 
% This function is utilized to produce offsprings from parent chromosomes.
% The genetic operators corssover and mutation which are carried out with
% slight modifications from the original design. For more information read
% the document enclosed. 
%
% parent_chromosome - the set of selected chromosomes.
% M - number of objective functions
% V - number of decision varaiables
% mu - distribution index for crossover (read the enlcosed pdf file)
% mum - distribution index for mutation (read the enclosed pdf file)
% l_limit - a vector of lower limit for the corresponding decsion variables
% u_limit - a vector of upper limit for the corresponding decsion variables
%
% The genetic operation is performed only on the decision variables, that
% is the first V elements in the chromosome vector. 

%  Copyright (c) 2009, Aravind Seshadri
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without 
%  modification, are permitted provided that the following conditions are 
%  met:
%
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%      
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%  POSSIBILITY OF SUCH DAMAGE.

[N,m] = size(parent_chromosome);

clear m
p = 1;
% Flags used to set if crossover and mutation were actually performed. 
was_crossover = 0;
was_mutation = 0;
for i = 1 : N
    % With 90 % probability perform crossover
    if rand(1) < 0.9
        % Initialize the children to be null vector.
        % Select the first parent
        parent_1 = round(N*rand(1));
        if parent_1 < 1
            parent_1 = 1;
        end
        % Select the second parent
        parent_2 = parent_1;
        while parent_1 == parent_2
            parent_2 = round(N*rand(1));
            if parent_2 < 1
                parent_2 = 1;
            end
        end
        % Make sure both the parents are not the same. 
%         while isequal(parent_chromosome(parent_1,:),parent_chromosome(parent_2,:))
%             parent_2 = round(N*rand(1));
%             if parent_2 < 1
%                 parent_2 = 1;
%             end
%         end
        % Get the chromosome information for each randomnly selected
        % parents
        parent_1 = parent_chromosome(parent_1,1:V);
        parent_2 = parent_chromosome(parent_2,1:V);
        child_1 = parent_1;
        child_2 = parent_2;
        % Perform corssover for each decision variable in the chromosome.
        u = rand(size(parent_1)) <= 0.5;
        u2 = rand(size(parent_1));
        u3 = rand(size(parent_1));
        y1 = min(parent_1,parent_2);
        y2 = max(parent_1,parent_2);
        beta = 1.0 + (2.0*(y1-l_limit)./(y2-y1));
        alpha = 2.0 +((mu*ones(size(parent_1))+1.0)).^(beta);
        betaq = beta;
        u4 = u2<=(1.0./alpha);
        betaq(u4) = (mu*ones(1,sum(u4))+1.0).^(rand(1,sum(u4)).*alpha(u4));
        betaq(~u4)= (mu*ones(1,sum(~u4))+1.0).^(1.0./(2.0-rand(1,sum(~u4)).*alpha(~u4)));
        child_1(u) = 0.5 * (y1(u)+y2(u)-betaq(u).*(y2(u)-y1(u)));
        
        beta = 1.0 + (2.0*(u_limit-y2)./(y2-y1));
        alpha = 2.0 + ((mu*ones(size(parent_1))+1.0)).^(beta);
        betaq = beta;
        u5 = u3<=(1.0./alpha);
        betaq(u5) = (mu*ones(1,sum(u5))+1.0).^(rand(1,sum(u5)).*alpha(u5));
        betaq(~u5)= (mu*ones(1,sum(~u5))+1.0).^(1.0./(2.0-rand(1,sum(~u5)).*alpha(~u5)));
        child_2(u) = 0.5 * (y1(u)+y2(u)+betaq(u).*(y2(u)-y1(u)));
        
        child_1(child_1 > u_limit) =  u_limit(child_1 > u_limit);
        child_2(child_2 > u_limit) =  u_limit(child_2 > u_limit);
        child_1(child_1 < l_limit) =  l_limit(child_1 < l_limit);
        child_2(child_2 < l_limit) =  l_limit(child_2 < l_limit);
        
        u6 = rand(size(parent_1)) <= 0.5;
        temp = child_1;
        child_1(u6) = child_2(u6);
        child_2(u6) = temp(u6);
        
        % Evaluate the objective function for the offsprings and as before
        % concatenate the offspring chromosome with objective value.
        child_1(:,V + 1: M + V) = evaluate_objective(child_1, M, V);
        child_2(:,V + 1: M + V) = evaluate_objective(child_2, M, V);
        % Set the crossover flag. When crossover is performed two children
        % are generate, while when mutation is performed only only child is
        % generated.
        was_crossover = 1;
        was_mutation = 0;
    % With 10 % probability perform mutation. Mutation is based on
    % polynomial mutation. 
    else
        % Select at random the parent.
        parent_3 = round(N*rand(1));
        if parent_3 < 1
            parent_3 = 1;
        end
        % Get the chromosome information for the randomnly selected parent.
        child_3 = parent_chromosome(parent_3,:);
        % Perform mutation on eact element of the selected parent.
        for j = 1 : V
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
           if child_3(j) > u_limit(j)
               child_3(j) = u_limit(j);
           elseif child_3(j) < l_limit(j)
               child_3(j) = l_limit(j);
           end
        end
        % Evaluate the objective function for the offspring and as before
        % concatenate the offspring chromosome with objective value.    
        child_3(:,V + 1: M + V) = evaluate_objective(child_3, M, V);
        % Set the mutation flag
        was_mutation = 1;
        was_crossover = 0;
    end
    % Keep proper count and appropriately fill the child variable with all
    % the generated children for the particular generation.
    if was_crossover
        child(p,:) = child_1;
        child(p+1,:) = child_2;
        was_cossover = 0;
        p = p + 2;
    elseif was_mutation
        child(p,:) = child_3(1,1 : M + V);
        was_mutation = 0;
        p = p + 1;
    end
end
f = child;
end
