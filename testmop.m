%   Author          : ChenZhe
%   Create time     : Friday, apr.1 2016
% 
%   Get test multi-objective problems from a given name.
%   The method get testing or benchmark problems for Multi-Objective
%   Optimization. (include  dynamic, stationary and Many).
%   The implemented problems included ZDT, DTLZ, ,WFG, UF, CF, UD, 
%   CD, FDA, et..
%   User get the corresponding test problem, which is an instance of class
%   mop, by passing the <problem name> and <optional dimension parameters>.
% 
%   Ref.
% 	[1]:Huband S, Hingston P, Barone L, While L, 2006, A review
%       of multiobjective test problems and a scalable test problem
%       toolkit. IEEE Transactions on Evolutionary Computation,
%       10(5), pp477-506.
% 	[2]:Zitzler, E., Deb, K., & Thiele, L. (2000). Comparison 
%       of multiobjective evolutionary algorithms: Empirical 
%       results. Evolutionary computation, 8(2), 173-195.
% 	[3]:Deb, K., Thiele, L., Laumanns, M., & Zitzler, E. (2005). 
%       Scalable test problems for evolutionary multiobjective 
%       optimization (pp. 105-145). Springer London.
%	[4]:Zhang, Q., Zhou, A., Zhao, S., Suganthan, P. N., Liu, W., &
%       Tiwari, S. (2008). Multiobjective optimization test instances
%       for the CEC 2009 special session and competition. University 
%       of Essex, Colchester, UK and Nanyang technological University,
%       Singapore, special session on performance assessment of 
%       multi-objective optimization algorithms, technical report, 264.
%	[5]:Farina, M., Deb, K., & Amato, P. (2004). Dynamic multiobjective 
%       optimization problems: test cases, approximations, and 
%       applications. Evolutionary Computation, IEEE Transactions on, 
%       8(5), 425-442.
%	[6]:Biswas, S., Das, S., Suganthan, P., & Coello Coello, C. (2014, 
%       July). Evolutionary multiobjective optimization in dynamic 
%       environments: A set of novel benchmark functions. In Evolutionary
%       Computation (CEC), 2014 IEEE Congress on (pp. 3192-3199). IEEE.
% 
%   Input:
%       testname    : (char array)
%       dimension   : (integer) 
%                    In DTLZ, the dim = M - 1 + k;
%                    In WFG, the dim = l + k;
% 
%   GLOBAL Input    : NOTE =.= When select corresponding test problems,
%                     you must assign the below variables, which is
%                     marked by KEYWORD 'global'.
%       dynamic pros: Detail information below
%                     itrCounter - current iter numbers
%                     step       - dynamic step size
%                     window     - dynamic window numbers
%       WFG pros    : Detail information below
%                     k          - no. of position-related parameters
%                     M          - no. of objectives
%                     l          - no. of distance-related parameters
%       DTLZ pros   : Detail information below
%                     M          - no. of objectives
%                     k          - variable which control the dims
% 
%   Output:
%       mop         : (struct) Detail information below
%                     name   - testname
%                     od     - objective dimension
%                     pd     - decision dimension
%                     domain - decision boundary constraints
%                     func   - ref to objective function
function mop = testmop(testname, dimension)
mop         =   struct('name',[],'od',[],'pd',[],'domain',[],'func',[]);
eval(['mop=',testname,'(mop,',num2str(dimension),');']);
end

%% ----------Stationary Multi-Objective Benchmark----------
% ---------------------------------------------------------
% ---------------------------------------------------------

%% KNO1 function generator
function p=kno1(p,dimension)
 p.name='KNO1';
 p.od = 2;
 p.pd = 2;
 p.domain= [0 3;0 3];
 p.func = @evaluate;
 
    %KNO1 evaluation function.
    function y = evaluate(x)
      y=zeros(2,1);
	  c = x(1)+x(2);
	  f = 9-(3*sin(2.5*c^0.5) + 3*sin(4*c) + 5 *sin(2*c+2));
	  g = (pi/2.0)*(x(1)-x(2)+3.0)/6.0;
	  y(1)= 20-(f*cos(g));
	  y(2)= 20-(f*sin(g)); 
    end
end

%% ----------ZDT series. Ref.[2]----------- 
%% ZDT1 function generator
function p=zdt1(p,dim)
 p.name='ZDT1';
 p.pd=dim;
 p.od=2;
 p.domain=[zeros(dim,1) ones(dim,1)];
 p.func=@evaluate;
 
    %zdt1 evaluation function.
    function y=evaluate(x)
        y=zeros(2,1);
        y(1) = x(1);
    	su = sum(x)-x(1);    
		g = 1 + 9 * su / (dim - 1);
		y(2) =g*(1 - sqrt(x(1) / g));
    end
end

%% ZDT2 function generator
function p=zdt2(p,dim)
 p.name='ZDT2';
 p.pd=dim;
 p.od=2;
 p.domain=[zeros(dim,1) ones(dim,1)];
 p.func=@evaluate;
 
    %zdt2 evaluation function.
    function y=evaluate(x)
        y=zeros(2,1);
        g = 1 + 9*sum(x(2:dim))/(dim-1);
        y(1) = x(1);
        y(2) = g*(1-(x(1)/g)^2);
    end
end

%% ZDT3 function generator
function p=zdt3(p,dim)
 p.name='ZDT3';
 p.pd=dim;
 p.od=2;
 p.domain=[zeros(dim,1) ones(dim,1)];
 p.func=@evaluate;
 
    %zdt3 evaluation function.
    function y=evaluate(x)
        y=zeros(2,1);
        g = 1 + 9*sum(x(2:dim))/(dim-1);
        y(1) = x(1);
        y(2) = g*(1-sqrt(x(1)/g) - x(1)/g*sin(10*pi*x(1)));
    end
end

%% ZDT4 function generator
function p=zdt4(p,dim)
 p.name='ZDT4';
 p.pd=dim;
 p.od=2;
 p.domain=[-5*ones(dim,1) 5*ones(dim,1)];
 p.domain(1,1) = 0;
 p.domain(1,2) = 1;
 p.func=@evaluate;
 
    %zdt4 evaluation function.
    function y=evaluate(x)
        y=zeros(2,1);
        g = 1 + 10*(10-1);
        for i = 2:10
            g = g+x(i)^2-10*cos(4*pi*x(i));
        end
        y(1) = x(1);
        y(2) = g*(1-sqrt(x(1)/g));
    end
end

%% ZDT6 function generator
function p=zdt6(p,dim)
 p.name='ZDT6';
 p.pd=dim;
 p.od=2;
 p.domain=[zeros(dim,1) ones(dim,1)];
 p.func=@evaluate;
 
    %zdt6 evaluation function.
    function y=evaluate(x)
        y=zeros(2,1);
        g = 1 + 9*(sum(x(2:dim))/(dim-1))^0.25;
        y(1) = 1 -exp(-4*x(1))*sin(6*pi*x(1))^6;
        y(2) = g*(1-(y(1)/g)^2);
    end
end

%% --------------DTLZ benchmark-------Ref.[3]----

%% DTLZ1 function generator
% as suggested by Deb
% k = 5 
function p=DTLZ1(p,dim)
 global M k;
 p.name='DTLZ1';
 p.pd=dim;
 p.od=M;
 p.domain=[zeros(dim,1) ones(dim,1)];
 p.func=@evaluate;
 
    %DTLZ evaluation function.
    function y=evaluate(x)
        n = (M-1) + k; %this is the default
        if size(x,1) ~= n
           error(['Using k = 5, it is required that the number of dimensions be '...
           'n = (M - 1) + k = %d in this case.'], n)
        end

        xm = x(n-k+1:end,:); %xm contains the last k variables
        g = 100*(k + sum((xm - 0.5).^2 - cos(20*pi*(xm - 0.5)),1));

        % Now, computes the functions. The first and the last will be
        % written separately to facilitate things
        f(1,:) = 1/2*prod(x(1:M-1,:),1).*(1 + g);
        for ii = 2:M-1
           f(ii,:) = 1/2*prod(x(1:M-ii,:),1).*(1 - x(M-ii+1,:)).*(1 + g);
        end
        f(M,:) = 1/2*(1 - x(1,:)).*(1 + g);
        y = f;
    end
end

%% DTLZ2 function generator
% as suggested by Deb
% k = 10;
function p=DTLZ2(p,dim)
 global M k;
 p.name='DTLZ2';
 p.pd=dim;
 p.od=M;
 p.domain=[zeros(dim,1) ones(dim,1)];
 p.func=@evaluate;
 
    %DTLZ evaluation function.
    function y=evaluate(x)
        n = (M-1) + k; %this is the default
        if size(x,1) ~= n
           error(['Using k = 10, it is required that the number of dimensions be'...
           ' n = (M - 1) + k = %d in this case.'], n)
        end

        xm = x(n-k+1:end,:); %xm contains the last k variables
        g = sum((xm - 0.5).^2, 1);

        % Computes the functions
        f(1,:) = (1 + g).*prod(cos(pi/2*x(1:M-1,:)),1);
        for ii = 2:M-1
           f(ii,:) = (1 + g) .* prod(cos(pi/2*x(1:M-ii,:)),1) .* ...
              sin(pi/2*x(M-ii+1,:));
        end
        f(M,:) = (1 + g).*sin(pi/2*x(1,:));
        y = f;
    end
end

%% DTLZ3 function generator
% k = 10
function p=DTLZ3(p,dim)
 global M k;
 p.name='DTLZ3';
 p.pd=dim;
 p.od=M;
 p.domain=[zeros(dim,1) ones(dim,1)];
 p.func=@evaluate;
 
    %DTLZ3 evaluation function.
    function y=evaluate(x)
        % Error check: the number of dimensions must be M-1+k
        n = (M-1) + k; %this is the default
        if size(x,1) ~= n
           error(['Using k = 10, it is required that the number of dimensions be'...
           ' n = (M - 1) + k = %d in this case.'], n)
        end

        xm = x(n-k+1:end,:); %xm contains the last k variables
        g = 100*(k + sum((xm - 0.5).^2 - cos(20*pi*(xm - 0.5)),1));

        % Computes the functions
        f(1,:) = (1 + g).*prod(cos(pi/2*x(1:M-1,:)),1);
        for ii = 2:M-1
           f(ii,:) = (1 + g) .* prod(cos(pi/2*x(1:M-ii,:)),1) .* ...
              sin(pi/2*x(M-ii+1,:));
        end
        f(M,:) = (1 + g).*sin(pi/2*x(1,:));
        y = f;
    end
end


%% DTLZ4 function generator
% k = 10;
function p=DTLZ4(p,dim)
 global M k;
 p.name='DTLZ4';
 p.pd=dim;
 p.od=M;
 p.domain=[zeros(dim,1) ones(dim,1)];
 p.func=@evaluate;
 
    %DTLZ evaluation function.
    function y=evaluate(x)
        alpha = 100;
        % Error check: the number of dimensions must be M-1+k
        n = (M-1) + k; %this is the default
        if size(x,1) ~= n
           error(['Using k = 10, it is required that the number of dimensions be'...
           ' n = (M - 1) + k = %d in this case.'], n)
        end

        xm = x(n-k+1:end,:); %xm contains the last k variables
        g = sum((xm - 0.5).^2, 1);

        % Computes the functions
        f(1,:) = (1 + g).*prod(cos(pi/2*x(1:M-1,:).^alpha),1);
        for ii = 2:M-1
           f(ii,:) = (1 + g) .* prod(cos(pi/2*x(1:M-ii,:).^alpha),1) .* ...
              sin(pi/2*x(M-ii+1,:).^alpha);
        end
        f(M,:) = (1 + g).*sin(pi/2*x(1,:).^alpha);
        y = f;
    end
end

%% DTLZ5 function generator
% k = 10;
function p=DTLZ5(p,dim)
 global M k;
 p.name='DTLZ5';
 p.pd=dim;
 p.od=M;
 p.domain=[zeros(dim,1) ones(dim,1)];
 p.func=@evaluate;
 
    %DTLZ5 evaluation function.
    function y=evaluate(x)
        % Error check: the number of dimensions must be M-1+k
        n = (M-1) + k; %this is the default
        if size(x,1) ~= n
           error(['Using k = 10, it is required that the number of dimensions be'...
           ' n = (M - 1) + k = %d in this case.'], n)
        end

        % There is a gr in the article. But, as used in the file from the authors,
        % gr = g 
        xm = x(n-k+1:end,:); %xm contains the last k variables
        g = sum((xm - 0.5).^2, 1); 

        theta(1,:) = pi/2*x(1,:);
        gr = g(ones(M-2,1),:); %replicates gr for the multiplication below
        theta(2:M-1,:) = pi./(4*(1+gr)) .* (1 + 2*gr.*x(2:M-1,:));

        % Finally, writes down the functions (there was a mistake in the article.
        % There is no pi/2 multiplication inside the cosine and sine functions)
        f(1,:) = (1 + g).*prod(cos(theta(1:M-1,:)),1);
        for ii = 2:M-1
           f(ii,:) = (1 + g) .* prod(cos(theta(1:M-ii,:)),1) .* ...
              sin(theta(M-ii+1,:));
        end
        f(M,:) = (1 + g).*sin(theta(1,:));
        y = f;
    end
end

%% DTLZ6 function generator
% k = 10;
function p=DTLZ6(p,dim)
 global M k;
 p.name='DTLZ6';
 p.pd=dim;
 p.od=M;
 p.domain=[zeros(dim,1) ones(dim,1)];
 p.func=@evaluate;
 
    %DTLZ6 evaluation function.
    function y=evaluate(x)
        % Error check: the number of dimensions must be M-1+k
        n = (M-1) + k; %this is the default
        if size(x,1) ~= n
           error(['Using k = 10, it is required that the number of dimensions be'...
           ' n = (M - 1) + k = %d in this case.'], n)
        end

        % There is a gr in the article. But, as used in the file from the authors,
        % gr = g 
        xm = x(n-k+1:end,:); %xm contains the last k variables
        g = sum(xm.^0.1,1); 

        theta(1,:) = pi/2*x(1,:);
        gr = g(ones(M-2,1),:); %replicates gr for the multiplication below
        theta(2:M-1,:) = pi./(4*(1+gr)) .* (1 + 2*gr.*x(2:M-1,:));

        % Finally, writes down the functions (there was a mistake in the article.
        % There is no pi/2 multiplication inside the cosine and sine functions)
        f(1,:) = (1 + g).*prod(cos(theta(1:M-1,:)),1);
        for ii = 2:M-1
           f(ii,:) = (1 + g) .* prod(cos(theta(1:M-ii,:)),1) .* ...
              sin(theta(M-ii+1,:));
        end
        f(M,:) = (1 + g).*sin(theta(1,:));
        y = f;
    end
end

%% DTLZ7 function generator
% k = 20;
function p=DTLZ7(p,dim)
 global M k;
 p.name='DTLZ7';
 p.pd=dim;
 p.od=M;
 p.domain=[zeros(dim,1) ones(dim,1)];
 p.func=@evaluate;
 
    %DTLZ evaluation function.
    function y=evaluate(x)
        % Error check: the number of dimensions must be M-1+k
        n = (M-1) + k; %this is the default
        if size(x,1) ~= n
           error(['Using k = 20, it is required that the number of dimensions be'...
           ' n = (M - 1) + k = %d in this case.'], n)
        end

        % Writes down the auxiliar function g
        xm = x(n-k+1:end,:); %xm contains the last k variables
        g = 1 + 9/k*sum(xm,1);

        % Now, computes the first M-1 objective functions
        f(1:M-1,:) = x(1:M-1,:);

        % The last function requires another auxiliar variable
        gaux = g(ones(M-1,1),:); %replicates the g function
        h = M - sum(f./(1+gaux).*(1 + sin(3*pi*f)),1);
        f(M,:) = (1 + g).*h;
        y = f;
    end
end

%% --------------WFG benchmark--------Ref.[1]----

%% WFG1 function generator
% dim = k + l;
function p=wfg1(p,dim)
 global k l M;
 p.name='WFG1';
 p.pd=dim;
 p.od=M;
 p.domain=[zeros(dim,1) 2*[1:dim]'];
 p.func=@evaluate;
 
    %WFG1 evaluation function.
    function y=evaluate(x)
        % Initialize
        [noSols, n, S, D, A, Y] = wfg_initialize(x, p.od, k, l, 1);

        % Apply first transformation.
        Ybar = Y;
        lLoop = k + 1;
        shiftA = 0.35;
        Ybar(:, lLoop:n) = s_linear(Ybar(:, lLoop:n), shiftA);
        
        % Apply second transformation.
        Ybarbar = Ybar;
        biasA = 0.8;
        biasB = 0.75;
        biasC = 0.85;
        Ybarbar(:, lLoop:n) = b_flat(Ybarbar(:, lLoop:n), biasA, biasB, ...
            biasC);
        
        % Apply third transformation.
        Ybarbarbar = Ybarbar;
        biasA = 0.02;
        Ybarbarbar = b_poly(Ybarbarbar, biasA);
        
        % Apply fourth transformation.
        T = NaN * ones(noSols, M);
        uLoop = M - 1;
        for i = 1:uLoop
            lBnd = 1+(i-1)*k/(M-1);
            uBnd = i*k/(M-1);
            weights = 2*(lBnd:uBnd);
            T(:,i) = r_sum(Ybarbarbar(:,lBnd:uBnd), weights);
        end
        T(:,M) = r_sum(Ybarbarbar(:,lLoop:n), 2*(lLoop:n));
        
        % Apply degeneracy constants.
        X = T;
        for i = 1:M-1
            X(:,i) = max(T(:,i), A(1, i)) .* (T(:,i) - 0.5) + 0.5;
        end
        
        % Generate objective values.
        fM = h_mixed(X(:,1), 1, 5);
        F = h_convex(X(:,1:uLoop));
        F(:,M) = fM;
        F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* F;
        y = F';
    end
end

%% WFG2 function generator
function p=wfg2(p,dim)
 global k l M;
 p.name='WFG2';
 p.pd=dim;
 p.od=M;
 p.domain=[zeros(dim,1) 2*[1:dim]'];
 p.func=@evaluate;
 
    %WFG2 evaluation function.
    function y=evaluate(x)
        % Initialize
        [noSols, n, S, D, A, Y] = wfg_initialize(x, p.od, k, l, 2);

        % Apply first transformation.
        Ybar = Y;
        lLoop = k + 1;
        shiftA = 0.35;
        Ybar(:, lLoop:n) = s_linear(Ybar(:, lLoop:n), shiftA);
        
        % Apply second transformation.
        Ybarbar = Ybar;
        uLoop = k + l/2;
        for i = lLoop:uLoop
            lBnd = k+2*(i-k)-1;
            uBnd = k+2*(i-k);
            Ybarbar(:,i) = r_nonsep(Ybar(:,lBnd:uBnd), 2);
        end
        
        % Apply third transformation.
        T = NaN * ones(noSols, M);
        uLoop = M - 1;
        weights = ones(1, k/(M-1));
        for i = 1:uLoop
            lBnd = 1+(i-1)*k/(M-1);
            uBnd = i*k/(M-1);
            T(:,i) = r_sum(Ybarbar(:,lBnd:uBnd), weights);
        end
        T(:,M) = r_sum(Ybarbar(:,lLoop:k+l/2), ones(1, (k+l/2)-lLoop+1));
        
        % Apply degeneracy constants.
        X = T;
        for i = 1:M-1
            X(:,i) = max(T(:,i), A(2, i)) .* (T(:,i) - 0.5) + 0.5;
        end
        
        % Generate objective values.
        fM = h_disc(X(:,1), 1, 1, 5);
        F = h_convex(X(:,1:uLoop));
        F(:,M) = fM;
        F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* F;
        y = F';
    end
end

%% WFG3 function generator
function p=wfg3(p,dim)
 global k l M;
 p.name='WFG3';
 p.pd=dim;
 p.od=M;
 p.domain=[zeros(dim,1) 2*[1:dim]'];
 p.func=@evaluate;
 
    %WFG3 evaluation function.
    function y=evaluate(x)
        % Initialize
        [noSols, n, S, D, A, Y] = wfg_initialize(x, p.od, k, l, 3);

        % Apply first transformation.
        Ybar = Y;
        lLoop = k + 1;
        shiftA = 0.35;
        Ybar(:, lLoop:n) = s_linear(Ybar(:, lLoop:n), shiftA);
        
        % Apply second transformation.
        Ybarbar = Ybar;
        uLoop = k + l/2;
        for i = lLoop:uLoop
            lBnd = k+2*(i-k)-1;
            uBnd = k+2*(i-k);
            Ybarbar(:,i) = r_nonsep(Ybar(:,lBnd:uBnd), 2);
        end
        
        % Apply third transformation.
        T = NaN * ones(noSols, M);
        uLoop = M - 1;
        weights = ones(1, k/(M-1));
        for i = 1:uLoop
            lBnd = 1+(i-1)*k/(M-1);
            uBnd = i*k/(M-1);
            T(:,i) = r_sum(Ybarbar(:,lBnd:uBnd), weights);
        end
        T(:,M) = r_sum(Ybarbar(:,lLoop:k+l/2), ones(1, (k+l/2)-lLoop+1));
        
        % Apply degeneracy constants.
        X = T;
        for i = 1:M-1
            X(:,i) = max(T(:,i), A(3, i)) .* (T(:,i) - 0.5) + 0.5;
        end
        
        % Generate objective values.
        F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* h_linear(X(:,1:uLoop));
        y = F';
    end
end

%% WFG4 function generator
function p=wfg4(p,dim)
 global k l M;    
 p.name='WFG4';
 p.pd=dim;
 p.od=M;
 p.domain=[zeros(dim,1) 2*[1:dim]'];
 p.func=@evaluate;
 
    %WFG4 evaluation function.
    function y=evaluate(x)
        % Initialize
        [noSols, n, S, D, A, Y] = wfg_initialize(x, p.od, k, l, 4);

        % Apply first transformation.
        if testNo == 4
            shiftA = 30;
            shiftB = 10;
            shiftC = 0.35;
            Ybar = s_multi(Y, shiftA, shiftB, shiftC);
        else
            shiftA = 0.35;
            shiftB = 0.001;
            shiftC = 0.05;
            Ybar = s_decep(Y, shiftA, shiftB, shiftC);
        end
        
        % Apply second transformation.
        T = NaN * ones(noSols, M);
        lLoop = k + 1;
        uLoop = M - 1;
        weights = ones(1, k/(M-1));
        for i = 1:uLoop
            lBnd = 1+(i-1)*k/(M-1);
            uBnd = i*k/(M-1);
            T(:,i) = r_sum(Ybar(:,lBnd:uBnd), weights);
        end
        T(:,M) = r_sum(Ybar(:,lLoop:n), ones(1, n-lLoop+1));
        
        % Apply degeneracy constants.
        X = T;
        for i = 1:M-1
            X(:,i) = max(T(:,i), A(4, i)) .* (T(:,i) - 0.5) + 0.5;
        end
        
        % Generate objective values.
        F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* h_concave(X(:,1:uLoop));
        y = F';
    end
end

%% WFG5 function generator
function p=wfg5(p,dim)
 global k l M;
 p.name='WFG5';
 p.pd=dim;
 p.od=M;
 p.domain=[zeros(dim,1) 2*[1:dim]'];
 p.func=@evaluate;
 
    %WFG5 evaluation function.
    function y=evaluate(x)
        % Initialize
        [noSols, n, S, D, A, Y] = wfg_initialize(x, p.od, k, l, 5);

        % Apply first transformation.
        shiftA = 0.35;
        shiftB = 0.001;
        shiftC = 0.05;
        Ybar = s_decep(Y, shiftA, shiftB, shiftC);
        
        % Apply second transformation.
        T = NaN * ones(noSols, M);
        lLoop = k + 1;
        uLoop = M - 1;
        weights = ones(1, k/(M-1));
        for i = 1:uLoop
            lBnd = 1+(i-1)*k/(M-1);
            uBnd = i*k/(M-1);
            T(:,i) = r_sum(Ybar(:,lBnd:uBnd), weights);
        end
        T(:,M) = r_sum(Ybar(:,lLoop:n), ones(1, n-lLoop+1));
        
        % Apply degeneracy constants.
        X = T;
        for i = 1:M-1
            X(:,i) = max(T(:,i), A(5, i)) .* (T(:,i) - 0.5) + 0.5;
        end
        
        % Generate objective values.
        F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* h_concave(X(:,1:uLoop));
        y = F';
    end
end

%% WFG6 function generator
function p=wfg6(p,dim)
 global k l M;
 p.name='WFG';
 p.pd=dim;
 p.od=M;
 p.domain=[zeros(dim,1) 2*[1:dim]'];
 p.func=@evaluate;
 
    %WFG6 evaluation function.
    function y=evaluate(x)
        % Initialize
        [noSols, n, S, D, A, Y] = wfg_initialize(x, p.od, k, l, 6);

        % Apply first transformation.
        Ybar = Y;
        lLoop = k + 1;
        shiftA = 0.35;
        Ybar(:, lLoop:n) = s_linear(Ybar(:, lLoop:n), shiftA);
        
        % Apply second transformation.
        T = NaN * ones(noSols, M);
        uLoop = M - 1;
        for i = 1:uLoop
            lBnd = 1+(i-1)*k/(M-1);
            uBnd = i*k/(M-1);
            T(:,i) = r_nonsep(Ybar(:,lBnd:uBnd), k/(M-1));
        end
        T(:,M) = r_nonsep(Ybar(:,k+1:k+l), l);
        
        % Apply degeneracy constants.
        X = T;
        for i = 1:M-1
            X(:,i) = max(T(:,i), A(6, i)) .* (T(:,i) - 0.5) + 0.5;
        end
        
        % Generate objective values.
        F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* h_concave(X(:,1:uLoop));
        y = F';
    end
end

%% WFG7 function generator
function p=wfg7(p,dim)
 global k l M;
 p.name='WFG7';
 p.pd=dim;
 p.od=M;
 p.domain=[zeros(dim,1) 2*[1:dim]'];
 p.func=@evaluate;
 
    %WFG7 evaluation function.
    function y=evaluate(x)
        % Initialize
        [noSols, n, S, D, A, Y] = wfg_initialize(x, p.od, k, l, 7);

        % Apply first transformation.
        Ybar = Y;
        biasA = 0.98 / 49.98;
        biasB = 0.02;
        biasC = 50;
        for i = 1:k
            Ybar(:,i) = b_param(Y(:,i), r_sum(Y(:,i+1:n), ones(1, n-i)), ...
                biasA, biasB, biasC);
        end
        
        % Apply second transformation.
        Ybarbar = Ybar;
        lLoop = k + 1;
        shiftA = 0.35;
        Ybarbar(:, lLoop:n) = s_linear(Ybar(:, lLoop:n), shiftA);
        
        % Apply third transformation.
        T = NaN * ones(noSols, M);
        lLoop = k + 1;
        uLoop = M - 1;
        weights = ones(1, k/(M-1));
        for i = 1:uLoop
            lBnd = 1+(i-1)*k/(M-1);
            uBnd = i*k/(M-1);
            T(:,i) = r_sum(Ybarbar(:,lBnd:uBnd), weights);
        end
        T(:,M) = r_sum(Ybarbar(:,lLoop:n), ones(1, n-lLoop+1));
        
        % Apply degeneracy constants.
        X = T;
        for i = 1:M-1
            X(:,i) = max(T(:,i), A(7, i)) .* (T(:,i) - 0.5) + 0.5;
        end
        
        % Generate objective values.
        F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* h_concave(X(:,1:uLoop));
        y = F';
    end
end

%% WFG8 function generator
function p=wfg8(p,dim)
 global k l M;
 p.name='WFG8';
 p.pd=dim;
 p.od=M;
 p.domain=[zeros(dim,1) 2*[1:dim]'];
 p.func=@evaluate;
 
    %WFG8 evaluation function.
    function y=evaluate(x)
        % Initialize
        [noSols, n, S, D, A, Y] = wfg_initialize(x, p.od, k, l, 8);

        % Apply first transformation.
        Ybar = Y;
        lLoop = k + 1;
        biasA = 0.98 / 49.98;
        biasB = 0.02;
        biasC = 50;
        for i = lLoop:n
            Ybar(:,i) = b_param(Y(:,i), r_sum(Y(:,1:i-1), ones(1, i-1)), ...
                biasA, biasB, biasC);
        end
        
        % Apply second transformation.
        Ybarbar = Ybar;
        shiftA = 0.35;
        Ybarbar(:, lLoop:n) = s_linear(Ybar(:, lLoop:n), shiftA);
        
        % Apply third transformation.
        T = NaN * ones(noSols, M);
        lLoop = k + 1;
        uLoop = M - 1;
        weights = ones(1, k/(M-1));
        for i = 1:uLoop
            lBnd = 1+(i-1)*k/(M-1);
            uBnd = i*k/(M-1);
            T(:,i) = r_sum(Ybarbar(:,lBnd:uBnd), weights);
        end
        T(:,M) = r_sum(Ybarbar(:,lLoop:n), ones(1, n-lLoop+1));
        
        % Apply degeneracy constants.
        X = T;
        for i = 1:M-1
            X(:,i) = max(T(:,i), A(8, i)) .* (T(:,i) - 0.5) + 0.5;
        end
        
        % Generate objective values.
        F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* h_concave(X(:,1:uLoop));
        y = F';
    end
end

%% WFG9 function generator
function p=wfg9(p,dim)
 global k l M;
 p.name='WFG';
 p.pd=dim;
 p.od=M;
 p.domain=[zeros(dim,1) 2*[1:dim]'];
 p.func=@evaluate;
 
    %WFG9 evaluation function.
    function y=evaluate(x)
        % Initialize
        [noSols, n, S, D, A, Y] = wfg_initialize(x, p.od, k, l, 9);

        % Apply first transformation.
        Ybar = Y;
        uLoop = n - 1;
        biasA = 0.98 / 49.98;
        biasB = 0.02;
        biasC = 50;
        for i = 1:uLoop
            Ybar(:,i) = b_param(Y(:,i), r_sum(Y(:,i+1:n), ones(1, n-i)), ...
                biasA, biasB, biasC);
        end
        
        % Apply second transformation.
        Ybarbar = Ybar;
        biasA = 0.35;
        biasB = 0.001;
        biasC = 0.05;
        Ybarbar(:,1:k) = s_decep(Ybar(:,1:k), biasA, biasB, biasC);
        biasA = 30;
        biasB = 95;
        biasC = 0.35;
        Ybarbar(:,k+1:n) = s_multi(Ybar(:,k+1:n), biasA, biasB, biasC);
        
        % Apply third transformation.
        T = NaN * ones(noSols, M);
        uLoop = M - 1;
        for i = 1:uLoop
            lBnd = 1+(i-1)*k/(M-1);
            uBnd = i*k/(M-1);
            T(:,i) = r_nonsep(Ybarbar(:,lBnd:uBnd), k/(M-1));
        end
        T(:,M) = r_nonsep(Ybarbar(:,k+1:k+l), l);
        
        % Apply degeneracy constants.
        X = T;
        for i = 1:M-1
            X(:,i) = max(T(:,i), A(9, i)) .* (T(:,i) - 0.5) + 0.5;
        end
        
        % Generate objective values.
        F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* h_concave(X(:,1:uLoop));
        y = F';
    end
end

%% WFG10 function generator
function p=wfg10(p,dim)
 global k l M;
 p.name='WFG10';
 p.pd=dim;
 p.od=M;
 p.domain=[zeros(dim,1) 2*[1:dim]'];
 p.func=@evaluate;
 
    %WFG10 evaluation function.
    function y=evaluate(x)
        % Initialize
        [noSols, n, S, D, A, Y] = wfg_initialize(x, p.od, k, l, 10);

        % Apply first transformation.
        shiftA = 30;
            shiftB = 10;
            shiftC = 0.35;
            Ybar = s_multi(Y, shiftA, shiftB, shiftC);
        
        % Apply second transformation.
        T = NaN * ones(noSols, M);
        lLoop = k + 1;
        uLoop = M - 1;
        weights = ones(1, k/(M-1));
        for i = 1:uLoop
            lBnd = 1+(i-1)*k/(M-1);
            uBnd = i*k/(M-1);
            T(:,i) = r_sum(Ybar(:,lBnd:uBnd), weights);
        end
        T(:,M) = r_sum(Ybar(:,lLoop:n), ones(1, n-lLoop+1));
        
        % Apply degeneracy constants.
        X = T;
        for i = 1:M-1
            X(:,i) = max(T(:,i), A(10, i)) .* (T(:,i) - 0.5) + 0.5;
        end
        % Generate objective values.
        F = rep(D*X(:,M), [1,M]) + rep(S, [noSols 1]) .* h_convex(X(:,1:uLoop));
        y = F';
    end
end

%% ----------%% CEC2009 series without constraints. Ref.[4]---------
%% UF1
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function p=uf1(p,dim)
 p.name='uf1';
 p.pd=dim;
 p.od=2;
 p.domain=[-1*ones(dim,1) ones(dim,1)];
 p.domain(1,1) = 0;
 p.func=@evaluate;
 
    %evaluation function.
    function y=evaluate(x)
        [dim, num]  = size(x);
        tmp         = zeros(dim,num);
        tmp(2:dim,:)= (x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]))).^2;
        tmp1        = sum(tmp(3:2:dim,:));  % odd index
        tmp2        = sum(tmp(2:2:dim,:));  % even index
        y(1,:)      = x(1,:)             + 2.0*tmp1/size(3:2:dim,2);
        y(2,:)      = 1.0 - sqrt(x(1,:)) + 2.0*tmp2/size(2:2:dim,2);
        clear tmp;
    end
end

%% UF2
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function p=uf2(p,dim)
 p.name='uf2';
 p.pd=dim;
 p.od=2;
 p.domain=[-1*ones(dim,1) ones(dim,1)];
 p.domain(1,1) = 0;
 p.func=@evaluate;
 
    %evaluation function.
    function y=evaluate(x)
        [dim, num]  = size(x);
        X1          = repmat(x(1,:),[dim-1,1]);
        A           = 6*pi*X1 + pi/dim*repmat((2:dim)',[1,num]);
        tmp         = zeros(dim,num);    
        tmp(2:dim,:)= (x(2:dim,:) - 0.3*X1.*(X1.*cos(4.0*A)+2.0).*cos(A)).^2;
        tmp1        = sum(tmp(3:2:dim,:));  % odd index
        tmp(2:dim,:)= (x(2:dim,:) - 0.3*X1.*(X1.*cos(4.0*A)+2.0).*sin(A)).^2;
        tmp2        = sum(tmp(2:2:dim,:));  % even index
        y(1,:)      = x(1,:)             + 2.0*tmp1/size(3:2:dim,2);
        y(2,:)      = 1.0 - sqrt(x(1,:)) + 2.0*tmp2/size(2:2:dim,2);
        clear X1 A tmp;
    end
end

%% UF3
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function p=uf3(p,dim)
 p.name='uf3';
 p.pd=dim;
 p.od=2;
 p.domain=[zeros(dim,1) ones(dim,1)];
 p.func=@evaluate;
 
    %evaluation function.
    function y=evaluate(x)
        [dim, num]   = size(x);
        Y            = zeros(dim,num);
        Y(2:dim,:)   = x(2:dim,:) - repmat(x(1,:),[dim-1,1]).^(0.5+1.5*(repmat((2:dim)',[1,num])-2.0)/(dim-2.0));
        tmp1         = zeros(dim,num);
        tmp1(2:dim,:)= Y(2:dim,:).^2;
        tmp2         = zeros(dim,num);
        tmp2(2:dim,:)= cos(20.0*pi*Y(2:dim,:)./sqrt(repmat((2:dim)',[1,num])));
        tmp11        = 4.0*sum(tmp1(3:2:dim,:)) - 2.0*prod(tmp2(3:2:dim,:)) + 2.0;  % odd index
        tmp21        = 4.0*sum(tmp1(2:2:dim,:)) - 2.0*prod(tmp2(2:2:dim,:)) + 2.0;  % even index
        y(1,:)       = x(1,:)             + 2.0*tmp11/size(3:2:dim,2);
        y(2,:)       = 1.0 - sqrt(x(1,:)) + 2.0*tmp21/size(2:2:dim,2);
        clear Y tmp1 tmp2;
    end
end

%% UF4
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function p=uf4(p,dim)
 p.name='uf4';
 p.pd=dim;
 p.od=2;
 p.domain=[-2*ones(dim,1) 2*ones(dim,1)];
 p.domain(1,1) = 0;
 p.domain(1,2) = 1;
 p.func=@evaluate;
 
    %evaluation function.
    function y=evaluate(x)
        [dim, num]  = size(x);
        Y           = zeros(dim,num);
        Y(2:dim,:)  = x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
        H           = zeros(dim,num);
        H(2:dim,:)  = abs(Y(2:dim,:))./(1.0+exp(2.0*abs(Y(2:dim,:))));
        tmp1        = sum(H(3:2:dim,:));  % odd index
        tmp2        = sum(H(2:2:dim,:));  % even index
        y(1,:)      = x(1,:)          + 2.0*tmp1/size(3:2:dim,2);
        y(2,:)      = 1.0 - x(1,:).^2 + 2.0*tmp2/size(2:2:dim,2);
        clear Y H;
    end
end

%% UF5
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function p=uf5(p,dim)
 p.name='uf5';
 p.pd=dim;
 p.od=2;
 p.domain=[-1*ones(dim,1) ones(dim,1)];
 p.domain(1,1) = 0;
 p.func=@evaluate;
 
    %evaluation function.
    function y=evaluate(x)
        N           = 10.0;
        E           = 0.1;
        [dim, num]  = size(x);
        Y           = zeros(dim,num);
        Y(2:dim,:)  = x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
        H           = zeros(dim,num);
        H(2:dim,:)  = 2.0*Y(2:dim,:).^2 - cos(4.0*pi*Y(2:dim,:)) + 1.0;
        tmp1        = sum(H(3:2:dim,:));  % odd index
        tmp2        = sum(H(2:2:dim,:));  % even index
        tmp         = (0.5/N+E)*abs(sin(2.0*N*pi*x(1,:)));
        y(1,:)      = x(1,:)      + tmp + 2.0*tmp1/size(3:2:dim,2);
        y(2,:)      = 1.0 - x(1,:)+ tmp + 2.0*tmp2/size(2:2:dim,2);
        clear Y H;
    end
end

%% UF6
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function p=uf6(p,dim)
 p.name='uf6';
 p.pd=dim;
 p.od=2;
 p.domain=[-1*ones(dim,1) ones(dim,1)];
 p.domain(1,1) = 0;
 p.func=@evaluate;
 
    %evaluation function.
    function y=evaluate(x)
        N            = 2.0;
        E            = 0.1;
        [dim, num]   = size(x);
        Y            = zeros(dim,num);
        Y(2:dim,:)  = x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
        tmp1         = zeros(dim,num);
        tmp1(2:dim,:)= Y(2:dim,:).^2;
        tmp2         = zeros(dim,num);
        tmp2(2:dim,:)= cos(20.0*pi*Y(2:dim,:)./sqrt(repmat((2:dim)',[1,num])));
        tmp11        = 4.0*sum(tmp1(3:2:dim,:)) - 2.0*prod(tmp2(3:2:dim,:)) + 2.0;  % odd index
        tmp21        = 4.0*sum(tmp1(2:2:dim,:)) - 2.0*prod(tmp2(2:2:dim,:)) + 2.0;  % even index
        tmp          = max(0,(1.0/N+2.0*E)*sin(2.0*N*pi*x(1,:)));
        y(1,:)       = x(1,:)       + tmp + 2.0*tmp11/size(3:2:dim,2);
        y(2,:)       = 1.0 - x(1,:) + tmp + 2.0*tmp21/size(2:2:dim,2);
        clear Y tmp1 tmp2;
    end
end

%% UF7
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function p=uf7(p,dim)
 p.name='uf7';
 p.pd=dim;
 p.od=2;
 p.domain=[-1*ones(dim,1) ones(dim,1)];
 p.domain(1,1) = 0;
 p.func=@evaluate;
 
    %evaluation function.
    function y=evaluate(x)
        [dim, num]  = size(x);
        Y           = zeros(dim,num);
        Y(2:dim,:)  = (x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]))).^2;
        tmp1        = sum(Y(3:2:dim,:));  % odd index
        tmp2        = sum(Y(2:2:dim,:));  % even index
        tmp         = (x(1,:)).^0.2;
        y(1,:)      = tmp       + 2.0*tmp1/size(3:2:dim,2);
        y(2,:)      = 1.0 - tmp + 2.0*tmp2/size(2:2:dim,2);
        clear Y;
    end
end

%% UF8
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function p=uf8(p,dim)
 p.name='uf8';
 p.pd=dim;
 p.od=3;
 p.domain=[-2*ones(dim,1) 2*ones(dim,1)];
 p.domain(1:2,1) = 0;
 p.domain(1:2,2) = 1;
 p.func=@evaluate;
 
    %evaluation function.
    function y=evaluate(x)
        [dim, num]  = size(x);
        Y           = zeros(dim,num);
        Y(3:dim,:)  = (x(3:dim,:) - 2.0*repmat(x(2,:),[dim-2,1]).*sin(2.0*pi*repmat(x(1,:),[dim-2,1]) + pi/dim*repmat((3:dim)',[1,num]))).^2;
        tmp1        = sum(Y(4:3:dim,:));  % j-1 = 3*k
        tmp2        = sum(Y(5:3:dim,:));  % j-2 = 3*k
        tmp3        = sum(Y(3:3:dim,:));  % j-0 = 3*k
        y(1,:)      = cos(0.5*pi*x(1,:)).*cos(0.5*pi*x(2,:)) + 2.0*tmp1/size(4:3:dim,2);
        y(2,:)      = cos(0.5*pi*x(1,:)).*sin(0.5*pi*x(2,:)) + 2.0*tmp2/size(5:3:dim,2);
        y(3,:)      = sin(0.5*pi*x(1,:))                     + 2.0*tmp3/size(3:3:dim,2);
        clear Y;
    end
end

%% UF9
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function p=uf9(p,dim)
 p.name='uf9';
 p.pd=dim;
 p.od=3;
 p.domain=[-2*ones(dim,1) 2*ones(dim,1)];
 p.domain(1:2,1) = 0;
 p.domain(1:2,2) = 1;
 p.func=@evaluate;
 
    %evaluation function.
    function y=evaluate(x)
        E           = 0.1;
        [dim, num]  = size(x);
        Y           = zeros(dim,num);
        Y(3:dim,:)  = (x(3:dim,:) - 2.0*repmat(x(2,:),[dim-2,1]).*sin(2.0*pi*repmat(x(1,:),[dim-2,1]) + pi/dim*repmat((3:dim)',[1,num]))).^2;
        tmp1        = sum(Y(4:3:dim,:));  % j-1 = 3*k
        tmp2        = sum(Y(5:3:dim,:));  % j-2 = 3*k
        tmp3        = sum(Y(3:3:dim,:));  % j-0 = 3*k
        tmp         = max(0,(1.0+E)*(1-4.0*(2.0*x(1,:)-1).^2));
        y(1,:)      = 0.5*(tmp+2*x(1,:)).*x(2,:)     + 2.0*tmp1/size(4:3:dim,2);
        y(2,:)      = 0.5*(tmp-2*x(1,:)+2.0).*x(2,:) + 2.0*tmp2/size(5:3:dim,2);
        y(3,:)      = 1-x(2,:)                       + 2.0*tmp3/size(3:3:dim,2);
        clear Y;
    end
end

%% UF10
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function p=uf10(p,dim)
 p.name='uf10';
 p.pd=dim;
 p.od=3;
p.domain=[-2*ones(dim,1) 2*ones(dim,1)];
p.domain(1:2,1) = 0;
p.domain(1:2,2) = 1;
 p.func=@evaluate;
 
    %evaluation function.
    function y=evaluate(x)
        [dim, num]  = size(x);
        Y           = zeros(dim,num);
        Y(3:dim,:)  = x(3:dim,:) - 2.0*repmat(x(2,:),[dim-2,1]).*sin(2.0*pi*repmat(x(1,:),[dim-2,1]) + pi/dim*repmat((3:dim)',[1,num]));
        H           = zeros(dim,num);
        H(3:dim,:)  = 4.0*Y(3:dim,:).^2 - cos(8.0*pi*Y(3:dim,:)) + 1.0;
        tmp1        = sum(H(4:3:dim,:));  % j-1 = 3*k
        tmp2        = sum(H(5:3:dim,:));  % j-2 = 3*k
        tmp3        = sum(H(3:3:dim,:));  % j-0 = 3*k
        y(1,:)      = cos(0.5*pi*x(1,:)).*cos(0.5*pi*x(2,:)) + 2.0*tmp1/size(4:3:dim,2);
        y(2,:)      = cos(0.5*pi*x(1,:)).*sin(0.5*pi*x(2,:)) + 2.0*tmp2/size(5:3:dim,2);
        y(3,:)      = sin(0.5*pi*x(1,:))                     + 2.0*tmp3/size(3:3:dim,2);
        clear Y H;
    end
end

%% ----------%% CEC2009 series with constraints. Ref.[4]---------
%% CF1
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function p=cf1(p,dim)
p.name   = 'cf1';
p.pd     = dim;
p.od     = 2;
p.domain = [zeros(dim,1) ones(dim,1)];
p.func   = @evaluate;

    function [y,c]=evaluate(x)
        a            = 1.0;
        N            = 10.0;
        [dim, num]   = size(x);
        Y            = zeros(dim,num);
        Y(2:dim,:)   = (x(2:dim,:) - repmat(x(1,:),[dim-1,1]).^(0.5+1.5*(repmat((2:dim)',[1,num])-2.0)/(dim-2.0))).^2;
        tmp1         = sum(Y(3:2:dim,:));% odd index
        tmp2         = sum(Y(2:2:dim,:));% even index 
        y(1,:)       = x(1,:)       + 2.0*tmp1/size(3:2:dim,2);
        y(2,:)       = 1.0 - x(1,:) + 2.0*tmp2/size(2:2:dim,2);
        c(1,:)       = y(1,:) + y(2,:) - a*abs(sin(N*pi*(y(1,:)-y(2,:)+1.0))) - 1.0;
        clear Y;
    end
end

%% CF2
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function p=cf2(p,dim)
p.name   = 'cf2';
p.pd     = dim;
p.od     = 2;
p.domain=[-1*ones(dim,1) ones(dim,1)];
p.domain(1,1) = 0;
p.func   = @evaluate;

    function [y,c]=evaluate(x)
        a           = 1.0;
        N           = 2.0;
        [dim, num]  = size(x);
        tmp         = zeros(dim,num);
        tmp(2:dim,:)= (x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]))).^2;
        tmp1        = sum(tmp(3:2:dim,:));  % odd index
        tmp(2:dim,:)= (x(2:dim,:) - cos(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]))).^2;
        tmp2        = sum(tmp(2:2:dim,:));  % even index
        y(1,:)      = x(1,:)             + 2.0*tmp1/size(3:2:dim,2);
        y(2,:)      = 1.0 - sqrt(x(1,:)) + 2.0*tmp2/size(2:2:dim,2);
        t           = y(2,:) + sqrt(y(1,:)) - a*sin(N*pi*(sqrt(y(1,:))-y(2,:)+1.0)) - 1.0;
        c(1,:)      = sign(t).*abs(t)./(1.0+exp(4.0*abs(t)));
        clear tmp;
    end
end

%% CF3
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function p=cf3(p,dim)
p.name   = 'cf3';
p.pd     = dim;
p.od     = 2;
p.domain=[-2*ones(dim,1) 2*ones(dim,1)];
p.domain(1,1) = 0;
p.domain(1,2) = 1;
p.func   = @evaluate;

    function [y,c]=evaluate(x)
        a            = 1.0;
        N            = 2.0;
        [dim, num]   = size(x);
        Y            = zeros(dim,num);
        Y(2:dim,:)   = x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
        tmp1         = zeros(dim,num);
        tmp1(2:dim,:)= Y(2:dim,:).^2;
        tmp2         = zeros(dim,num);
        tmp2(2:dim,:)= cos(20.0*pi*Y(2:dim,:)./sqrt(repmat((2:dim)',[1,num])));
        tmp11        = 4.0*sum(tmp1(3:2:dim,:)) - 2.0*prod(tmp2(3:2:dim,:)) + 2.0;  % odd index
        tmp21        = 4.0*sum(tmp1(2:2:dim,:)) - 2.0*prod(tmp2(2:2:dim,:)) + 2.0;  % even index
        y(1,:)       = x(1,:)          + 2.0*tmp11/size(3:2:dim,2);
        y(2,:)       = 1.0 - x(1,:).^2 + 2.0*tmp21/size(2:2:dim,2);
        c(1,:)       = y(2,:) + y(1,:).^2 - a*sin(N*pi*(y(1,:).^2-y(2,:)+1.0)) - 1.0;   
        clear Y tmp1 tmp2;
    end
end

%% CF4
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function p=cf4(p,dim)
p.name   = 'cf4';
p.pd     = dim;
p.od     = 2;
p.domain=[-2*ones(dim,1) 2*ones(dim,1)];
p.domain(1,1) = 0;
p.domain(1,2) = 1;
p.func   = @evaluate;

    function [y,c]=evaluate(x)
        [dim, num]  = size(x);
        tmp         = zeros(dim,num);
        tmp(2:dim,:)= x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
        tmp1        = sum(tmp(3:2:dim,:).^2);  % odd index
        tmp2        = sum(tmp(4:2:dim,:).^2);  % even index
        index1      = tmp(2,:) < (1.5-0.75*sqrt(2.0));
        index2      = tmp(2,:)>= (1.5-0.75*sqrt(2.0));
        tmp(2,index1) = abs(tmp(2,index1));
        tmp(2,index2) = 0.125 + (tmp(2,index2)-1.0).^2;
        y(1,:)      = x(1,:)                  + tmp1;
        y(2,:)      = 1.0 - x(1,:) + tmp(2,:) + tmp2;
        t           = x(2,:) - sin(6.0*pi*x(1,:)+2.0*pi/dim) - 0.5*x(1,:) + 0.25;
        c(1,:)      = sign(t).*abs(t)./(1.0+exp(4.0*abs(t)));
        clear tmp index1 index2;
    end
end

%% CF5
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function p=cf5(p,dim)
p.name   = 'cf5';
p.pd     = dim;
p.od     = 2;
p.domain=[-2*ones(dim,1) 2*ones(dim,1)];
p.domain(1,1) = 0;
p.domain(1,2) = 1;
p.func   = @evaluate;

    function [y,c]=evaluate(x)
        [dim, num]  = size(x);
        tmp         = zeros(dim,num);
        tmp(2:dim,:)= x(2:dim,:) - 0.8*repmat(x(1,:),[dim-1,1]).*cos(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
        tmp1        = sum(2.0*tmp(3:2:dim,:).^2-cos(4.0*pi*tmp(3:2:dim,:))+1.0);  % odd index
        tmp(2:dim,:)= x(2:dim,:) - 0.8*repmat(x(1,:),[dim-1,1]).*sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));    
        tmp2        = sum(2.0*tmp(4:2:dim,:).^2-cos(4.0*pi*tmp(4:2:dim,:))+1.0);  % even index
        index1      = tmp(2,:) < (1.5-0.75*sqrt(2.0));
        index2      = tmp(2,:)>= (1.5-0.75*sqrt(2.0));
        tmp(2,index1) = abs(tmp(2,index1));
        tmp(2,index2) = 0.125 + (tmp(2,index2)-1.0).^2;
        y(1,:)      = x(1,:)                  + tmp1;
        y(2,:)      = 1.0 - x(1,:) + tmp(2,:) + tmp2;
        c(1,:)      = x(2,:) - 0.8*x(1,:).*sin(6.0*pi*x(1,:)+2.0*pi/dim) - 0.5*x(1,:) + 0.25;
        clear tmp;
    end
end

%% CF6
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function p=cf6(p,dim)
p.name   = 'cf6';
p.pd     = dim;
p.od     = 2;
p.domain=[-2*ones(dim,1) 2*ones(dim,1)];
p.domain(1,1) = 0;
p.domain(1,2) = 1;
p.func   = @evaluate;

    function [y,c]=evaluate(x)
        [dim, num]  = size(x);
        tmp         = zeros(dim,num);
        tmp(2:dim,:)= x(2:dim,:) - 0.8*repmat(x(1,:),[dim-1,1]).*cos(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
        tmp1        = sum(tmp(3:2:dim,:).^2);  % odd index
        tmp(2:dim,:)= x(2:dim,:) - 0.8*repmat(x(1,:),[dim-1,1]).*sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));    
        tmp2        = sum(tmp(2:2:dim,:).^2);  % even index
        y(1,:)      = x(1,:)            + tmp1;
        y(2,:)      = (1.0 - x(1,:)).^2 + tmp2;
        tmp         = 0.5*(1-x(1,:))-(1-x(1,:)).^2;
        c(1,:)      = x(2,:) - 0.8*x(1,:).*sin(6.0*pi*x(1,:)+2*pi/dim) - sign(tmp).*sqrt(abs(tmp));
        tmp         = 0.25*sqrt(1-x(1,:))-0.5*(1-x(1,:));
        c(2,:)      = x(4,:) - 0.8*x(1,:).*sin(6.0*pi*x(1,:)+4*pi/dim) - sign(tmp).*sqrt(abs(tmp));    
        clear tmp;
    end
end

%% CF7
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function p=cf7(p,dim)
p.name   = 'cf7';
p.pd     = dim;
p.od     = 2;
p.domain=[-2*ones(dim,1) 2*ones(dim,1)];
p.domain(1,1) = 0;
p.domain(1,2) = 1;
p.func   = @evaluate;

    function [y,c]=evaluate(x)
        [dim, num]  = size(x);
        tmp         = zeros(dim,num);
        tmp(2:dim,:)= x(2:dim,:) - cos(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
        tmp1        = sum(2.0*tmp(3:2:dim,:).^2-cos(4.0*pi*tmp(3:2:dim,:))+1.0);  % odd index
        tmp(2:dim,:)= x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
        tmp2        = sum(2.0*tmp(6:2:dim,:).^2-cos(4.0*pi*tmp(6:2:dim,:))+1.0);  % even index
        tmp(2,:)    = tmp(2,:).^2;
        tmp(4,:)    = tmp(4,:).^2;
        y(1,:)      = x(1,:)                                  + tmp1;
        y(2,:)      = (1.0 - x(1,:)).^2 + tmp(2,:) + tmp(4,:) + tmp2;
        tmp         = 0.5*(1-x(1,:))-(1-x(1,:)).^2;
        c(1,:)      = x(2,:) - sin(6.0*pi*x(1,:)+2*pi/dim) - sign(tmp).*sqrt(abs(tmp));
        tmp         = 0.25*sqrt(1-x(1,:))-0.5*(1-x(1,:));
        c(2,:)      = x(4,:) - sin(6.0*pi*x(1,:)+4*pi/dim) - sign(tmp).*sqrt(abs(tmp));    
        clear tmp;
    end
end

%% CF8
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function p=cf8(p,dim)
p.name   = 'cf8';
p.pd     = dim;
p.od     = 3;
p.domain=[-4*ones(dim,1) 4*ones(dim,1)];
p.domain(1:2,1) = 0;
p.domain(1:2,2) = 1;
p.func   = @evaluate;

    function [y,c]=evaluate(x)
        N           = 2.0;
        a           = 4.0;
        [dim, num]  = size(x);
        Y           = zeros(dim,num);
        Y(3:dim,:)  = (x(3:dim,:) - 2.0*repmat(x(2,:),[dim-2,1]).*sin(2.0*pi*repmat(x(1,:),[dim-2,1]) + pi/dim*repmat((3:dim)',[1,num]))).^2;
        tmp1        = sum(Y(4:3:dim,:));  % j-1 = 3*k
        tmp2        = sum(Y(5:3:dim,:));  % j-2 = 3*k
        tmp3        = sum(Y(3:3:dim,:));  % j-0 = 3*k
        y(1,:)      = cos(0.5*pi*x(1,:)).*cos(0.5*pi*x(2,:)) + 2.0*tmp1/size(4:3:dim,2);
        y(2,:)      = cos(0.5*pi*x(1,:)).*sin(0.5*pi*x(2,:)) + 2.0*tmp2/size(5:3:dim,2);
        y(3,:)      = sin(0.5*pi*x(1,:))                     + 2.0*tmp3/size(3:3:dim,2);
        c(1,:)      = (y(1,:).^2+y(2,:).^2)./(1.0-y(3,:).^2) - a*abs(sin(N*pi*((y(1,:).^2-y(2,:).^2)./(1.0-y(3,:).^2)+1.0))) - 1.0;
        clear Y;
    end
end

%% CF9
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function p=cf9(p,dim)
p.name   = 'cf9';
p.pd     = dim;
p.od     = 3;
p.domain=[-2*ones(dim,1) 2*ones(dim,1)];
p.domain(1:2,1) = 0;
p.domain(1:2,2) = 1;
p.func   = @evaluate;

    function [y,c]=evaluate(x)
        N           = 2.0;
        a           = 3.0;
        [dim, num]  = size(x);
        Y           = zeros(dim,num);
        Y(3:dim,:)  = (x(3:dim,:) - 2.0*repmat(x(2,:),[dim-2,1]).*sin(2.0*pi*repmat(x(1,:),[dim-2,1]) + pi/dim*repmat((3:dim)',[1,num]))).^2;
        tmp1        = sum(Y(4:3:dim,:));  % j-1 = 3*k
        tmp2        = sum(Y(5:3:dim,:));  % j-2 = 3*k
        tmp3        = sum(Y(3:3:dim,:));  % j-0 = 3*k
        y(1,:)      = cos(0.5*pi*x(1,:)).*cos(0.5*pi*x(2,:)) + 2.0*tmp1/size(4:3:dim,2);
        y(2,:)      = cos(0.5*pi*x(1,:)).*sin(0.5*pi*x(2,:)) + 2.0*tmp2/size(5:3:dim,2);
        y(3,:)      = sin(0.5*pi*x(1,:))                     + 2.0*tmp3/size(3:3:dim,2);
        c(1,:)      = (y(1,:).^2+y(2,:).^2)./(1.0-y(3,:).^2) - a*sin(N*pi*((y(1,:).^2-y(2,:).^2)./(1.0-y(3,:).^2)+1.0)) - 1.0;
        clear Y;
    end
end

%% CF10
% x and y are columnwise, the imput x must be inside the search space and
% it could be a matrix
function p=cf10(p,dim)
p.name   = 'cf10';
p.pd     = dim;
p.od     = 3;
p.domain=[-2*ones(dim,1) 2*ones(dim,1)];
p.domain(1:2,1) = 0;
p.domain(1:2,2) = 1;
p.func   = @evaluate;

    function [y,c]=evaluate(x)
        a           = 1.0;
        N           = 2.0;
        [dim, num]  = size(x);
        Y           = zeros(dim,num);
        Y(3:dim,:)  = x(3:dim,:) - 2.0*repmat(x(2,:),[dim-2,1]).*sin(2.0*pi*repmat(x(1,:),[dim-2,1]) + pi/dim*repmat((3:dim)',[1,num]));
        H           = zeros(dim,num);
        H(3:dim,:)  = 4.0*Y(3:dim,:).^2 - cos(8.0*pi*Y(3:dim,:)) + 1.0;
        tmp1        = sum(H(4:3:dim,:));  % j-1 = 3*k
        tmp2        = sum(H(5:3:dim,:));  % j-2 = 3*k
        tmp3        = sum(H(3:3:dim,:));  % j-0 = 3*k
        y(1,:)      = cos(0.5*pi*x(1,:)).*cos(0.5*pi*x(2,:)) + 2.0*tmp1/size(4:3:dim,2);
        y(2,:)      = cos(0.5*pi*x(1,:)).*sin(0.5*pi*x(2,:)) + 2.0*tmp2/size(5:3:dim,2);
        y(3,:)      = sin(0.5*pi*x(1,:))                     + 2.0*tmp3/size(3:3:dim,2);
        c(1,:)      = (y(1,:).^2+y(2,:).^2)./(1.0-y(3,:).^2) - a*sin(N*pi*((y(1,:).^2-y(2,:).^2)./(1.0-y(3,:).^2)+1.0)) - 1.0;
        clear Y H;
    end
end

%% -----------Dynamic Multi-Objective Benchmark----------
% -------------------------------------------------------
% -------------------------------------------------------

%% ---------FDA series. Ref.[5]----------
%% FDA1
% x are columnwise and y are rowwise, the input x must be inside the search space and
% it could be a matrxix 
% time window for which the process remains constant, is recommended as 5
% step controls the distance between 2 paretos nT, is recommended as 10
% n is recommended as 20
function p=fda1(p,dim)
p.name   = 'fda1';
p.pd     = dim;
p.od     = 2;
p.domain = [-1*ones(dim,1) ones(dim,1)];
p.domain(1,1) = 0;
p.func   = @evaluate;

    function y=evaluate(x)
        global itrCounter step window;
        x            =x';
        n            =length(x);
        f1           =x(1);
        t            =(floor(itrCounter/window))/step;
        G            =sin(0.5*pi*t);
        temp         =x(2:n);
        Gtemp        =G*ones(n-1,1);
        k            =(temp-Gtemp).^2;
        arbit        =sum(k);
        g            =1+arbit;
        f2           =g*(1-(f1/g)^0.5);
        y            =[f1,f2];
        clear temp Gtemp k;
    end
end 

%% FDA2
% x are columnwise and y are rowwise, the input x must be inside the search space and
% it could be a matrxix 
% time window for which the process remains constant
% step controls the distance between 2 paretos nT
function p=fda2(p,dim)
p.name   = 'fda2';
p.pd     = dim;
p.od     = 2;
p.domain = [-1*ones(dim,1) ones(dim,1)];
p.domain(1,1) = 0;
p.func   = @evaluate;

    function y=evaluate(x)
    	global itrCounter params step window;
        x            =x';
        n            =length(x);
        f1           =x(1);
        t            =(floor(itrCounter/window))/step;
        H            =0.75+0.7*sin(0.5*pi*t);
        tempsize     =(n-1)/2;
        temp         =x(2:tempsize+1);
        temp2        =x(end-tempsize+1:end);
        g            =1+sum(temp.^2);
        Htemp        =H*ones(tempsize,1);
        k            =(temp2-Htemp).^2;
        arbit        =(sum(k)+H)^-1;
        f2           =g*((1-f1/g)^arbit);
        y            =[f1,f2];
        clear temp temp2 Htemp k;
    end
end 

%% FDA2-nsga2
% x are columnwise and y are rowwise, the input x must be inside the search space and
% it could be a matrxix 
% time window for which the process remains constant
% step controls the distance between 2 paretos nT
function p=fda2_nsga2(p,dim)
p.name   = 'fda2_nsga2';
p.pd     = dim;
p.od     = 2;
p.domain = [-1*ones(dim,1) ones(dim,1)];
p.domain(1,1) = 0;
p.func   = @evaluate;

    function y=evaluate(x)
    	global itrCounter params step window;
        x            =x';
        n            =length(x);
        f1           =x(1);
        t            =2*floor(itrCounter/window)*(window/(step*window-window));
%         t            =(floor(itrCounter/window))/step;
        H            =2*sin(0.5*pi*(t-1));
        temp         =x(2:6);
        temp2        =x(7:13);
        g            =1+sum(temp.^2);   
        Htemp        =sum((temp2-H/4).^2);
        h            =1-(f1/g)^(2^(H+Htemp));
        f2           =g*h;
        y            =[f1,f2];
        clear temp temp2 Htemp k h;
    end
end 

%% FDA3
% x are columnwise and y are rowwise, the input x must be inside the search space and
% it could be a matrxix 
% time window for which the process remains constant
% step controls the distance between 2 paretos nT
function p=fda3(p,dim)
p.name   = 'fda3';
p.pd     = dim;
p.od     = 2;
p.domain = [-1*ones(dim,1) ones(dim,1)];
p.domain(1,1) = 0;
p.func   = @evaluate;

    function y=evaluate(x)
        global itrCounter step window;
        x            =x';
        n            =length(x);
        t            =(floor(itrCounter/window))/step;
        G            =abs(sin(0.5*3.14*t));
        F            =2*G;  %10^(2*sin(0.5*pi*t))
        f1           =x(1)^F;
        Gtemp        =G*ones(n-1,1);
        temp         =x(2:n);
        g            =1+G+sum((temp-Gtemp).^2);
        f2           =g*(1-(f1/g)^0.5);
        y            =[f1,f2];
        clear Gtemp temp;
    end
end 

%% FDA4
% x are columnwise and y are rowwise, the input x must be inside the search space and
% it could be a matrxix 
% time window for which the process remains constant
% step controls the distance between 2 paretos nT
function p=fda4(p,dim)
p.name   = 'fda4';
p.pd     = dim;
p.od     = dim-9;
p.domain = [zeros(dim,1) ones(dim,1)];
p.func   = @evaluate;

    function y=evaluate(x)
        global itrCounter step window;
        x            =x';
        n            =length(x);
        t            =(floor(itrCounter/window))/step;
        temp         =x(p.od:end);
        G            =abs(sin(0.5*pi*t));
        gtemp        =(temp-G).^2;
        g            =sum(gtemp);
        f1temp       =cos((x(1:p.od-1)*pi)/2);
        f1           =(1+g)*prod(f1temp);
        f            =[f1];
        for k = 2:p.od-1
            fxtemp      =cos((x(1:p.od-k)*pi)/2);
            fx          =(1+g)*prod(fxtemp)*sin((x(p.od-k+1)*pi)/2);
            f           =[f,fx];
            clear fxtemp fx;
        end
        fxtemp       =sin((x(1)*pi)/2);
        fx           =(1+g)*fxtemp;
        f            =[f,fx];
        y            =f;
        clear temp gtemp f1temp fxtemp fx f;
    end
end 

%% FDA5
% x agre columnwise and y are rowwise, the input x must be inside the search space and
% it could be a matrxix 
% time window for which the process remains constant
% step controls the distance between 2 paretos nT
function p=fda5(p,dim)
p.name   = 'fda5';
p.pd     = dim;
p.od     = 2;
p.domain = [zeros(dim,1) ones(dim,1)];
p.func   = @evaluate;

    function y=evaluate(x)
        global itrCounter step window;
        n            =length(x);
        t            =(floor(itrCounter/window))/step;
        temp         =x(p.od:end);
        G            =abs(sin(0.5*pi*t));
        gtemp        =(temp-G).^2;
        g            =G+sum(gtemp);
        F            =1+100*(sin(0.5*pi*t).^4);
        y            =x(1:p.od-1).^F;
        f1temp       =cos((y(1:p.od-1)*pi)/2);
        f1           =(1+g)*prod(f1temp);
        f            =[f1];
        for k = 2:p.od-1
          fxtemp      =cos((y(1:p.od-k)*pi)/2);
          fx          =(1+g)*prod(fxtemp)*sin((y(p.od-k+1)*pi)/2);
          f           =[f;fx]
          clear fxtemp fx;
        end
        fxtemp       =sin((y(1)*pi)/2);
        fx           =(1+g)*fxtemp;
        f            =[f;fx];
        y            =f;
        clear temp gtemp f1temp fxtemp fx f;
    end
end 

%% ---------CEC2014 Benchmark without constraintrs. Ref.[6]----------
%% ud1
% x agre columnwise and y are rowwise, the input x must be inside the search space and
% it could be a matrxix 
% time window for which the process remains constant
% step controls the distance between 2 paretos nT
function p=ud1(p,dim)
p.name   = 'ud1';
p.pd     = dim;
p.od     = 2;
p.domain = [-1*ones(dim,1) ones(dim,1)];
p.domain(1,1) = 0;
p.func   = @evaluate;

    function y=evaluate(x)
        global itrCounter step window;
        [dim, num]  = size(x);
        t           =(floor(itrCounter/window))/step;
        G           =sin(0.5*pi*t);
        Gtemp       =G*ones(1,num);
        F           =ceil(dim*G);
        Ftemp       =F*ones(dim-1,num);
        tmp         = zeros(dim,num);
        tmp(2:dim,:)= (x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(Ftemp+repmat((2:dim)',[1,num])))).^2;
        tmp1        = sum(tmp(3:2:dim,:));  % odd index
        tmp(2:dim,:)= (x(2:dim,:)- sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(Ftemp+repmat((2:dim)',[1,num])))).^2;
        tmp2        = sum(tmp(2:2:dim,:));  % even index
        y(1,:)      = x(1,:)+abs(Gtemp) + 2.0*tmp1/size(3:2:dim,2);
        y(2,:)      = 1.0 - (x(1,:)).^1+abs(Gtemp) + 2.0*tmp2/size(2:2:dim,2);
        clear tmp tmp1 tmp2;
        f=[y(1,:);y(2,:)];
        y = f;
        clear G Gtemp F Ftemp f;
    end
end

%% ud2
% x agre columnwise and y are rowwise, the input x must be inside the search space and
% it could be a matrxix 
% time window for which the process remains constant
% step controls the distance between 2 paretos nT
function p=ud2(p,dim)
p.name   = 'ud2';
p.pd     = dim;
p.od     = 2;
p.domain = [-2*ones(dim,1) 2*ones(dim,1)];
p.domain(1,1) = 0;
p.domain(1,2) = 1;
p.func   = @evaluate;

    function y=evaluate(x)
        global itrCounter step window;
        [dim, num]  = size(x);
        t           =(floor(itrCounter/window))/step;
        G           =sin(0.5*pi*t);
        Gtemp       =G*ones(1,num);
        tmp         = zeros(dim,num);
        tmp(2:dim,:)= (x(2:dim,:)-Gtemp - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(repmat((2:dim)',[1,num])))).^2;
        tmp1        = sum(tmp(3:2:dim,:));  % odd index
        tmp(2:dim,:)= (x(2:dim,:)-Gtemp - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(repmat((2:dim)',[1,num])))).^2;
        tmp2        = sum(tmp(2:2:dim,:));  % even index
        y(1,:)      = x(1,:)+abs(Gtemp) + 2.0*tmp1/size(3:2:dim,2);
        y(2,:)      = 1.0 - (x(1,:).^1)+abs(Gtemp) + 2.0*tmp2/size(2:2:dim,2);
        clear tmp tmp1 tmp2;
        f           =[y(1,:);y(2,:)];
        y           =f;
        clear f G Gtemp F Ftemp; 
    end
end

%% ud3
% x agre columnwise and y are rowwise, the input x must be inside the search space and
% it could be a matrxix 
% time window for which the process remains constant
% step controls the distance between 2 paretos nT
function p=ud3(p,dim)
p.name   = 'ud3';
p.pd     = dim;
p.od     = 2;
p.domain = [-1*ones(dim,1) 2*ones(dim,1)];
p.domain(1,1) = 0;
p.domain(1,2) = 1;
p.func   = @evaluate;

    function y=evaluate(x)
        global itrCounter step window;
        [dim, num]  = size(x);
        t           =(floor(itrCounter/window))/step;
        G           =sin(0.5*pi*t);
        Gtemp       =G*ones(1,num);
        Y           = zeros(dim,num);
        Y(2:dim,:)  = (x(2:dim,:) - repmat(x(1,:),[dim-1,1]).^(1+Gtemp+1.5*(repmat((2:dim)',[1,num])-2.0)/(dim-2.0))-Gtemp).^2;
        tmp1        = sum(Y(3:2:dim,:));% odd index
        tmp2        = sum(Y(2:2:dim,:));% even index 
        y(1,:)      = x(1,:)+abs(Gtemp) + 2.0*tmp1/size(3:2:dim,2);
        y(2,:)      = 1.0 - (x(1,:)).^1+abs(Gtemp) + 2.0*tmp2/size(2:2:dim,2);
        clear tmp tmp1 tmp2 Y G Gtemp;
    end
end

%% ud4
% x agre columnwise and y are rowwise, the input x must be inside the search space and
% it could be a matrxix 
% time window for which the process remains constant
% step controls the distance between 2 paretos nT
function p=ud4(p,dim)
p.name   = 'ud4';
p.pd     = dim;
p.od     = 2;
p.domain = [-1*ones(dim,1) ones(dim,1)];
p.domain(1,1) = 0;
p.func   = @evaluate;

    function y=evaluate(x)
        global itrCounter step window;
        [dim, num]  = size(x);
        t           =(floor(itrCounter/window))/step;
        G           =sin(0.5*pi*t);
        Gtemp       =G*ones(1,num);
        N           = 10.0;
        E           = 0.1;
        Y           = zeros(dim,num);
        Y(2:dim,:)  = x(2:dim,:)- sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
        H           = zeros(dim,num);
        H(2:dim,:)  = 2.0*Y(2:dim,:).^2 - cos(4.0*pi*Y(2:dim,:)) + 1.0;
        tmp1        = sum(H(3:2:dim,:));  % odd index
        tmp2        = sum(H(2:2:dim,:));  % even index
        tmp         = (0.5/N+E)*abs(sin(2.0*N*pi*x(1,:))-abs(2*N*Gtemp));
        y(1,:)      = x(1,:)      + tmp + 2.0*tmp1/size(3:2:dim,2);
        y(2,:)      = 1.0 - x(1,:)+ tmp + 2.0*tmp2/size(2:2:dim,2);
        clear G Gtemp Y H tmp tmp1 tmp2;
    end
end

%% ud5
% x agre columnwise and y are rowwise, the input x must be inside the search space and
% it could be a matrxix 
% time window for which the process remains constant
% step controls the distance between 2 paretos nT
function p=ud5(p,dim)
p.name   = 'ud5';
p.pd     = dim;
p.od     = 2;
p.domain = [-1*ones(dim,1) ones(dim,1)];
p.domain(1,1) = 0;
p.func   = @evaluate;

    function y=evaluate(x)
        global itrCounter step window;
        [dim, num]   = size(x);
        t            =(floor(itrCounter/window))/step;
        G            =sin(0.5*pi*t);
        Gtemp        =G*ones(1,num);
        N            = 2.0;
        E            = 0.1;
        Y            = zeros(dim,num);
        Y(2:dim,:)  = x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
        tmp1         = zeros(dim,num);
        tmp1(2:dim,:)= Y(2:dim,:).^2;
        tmp2         = zeros(dim,num);
        tmp2(2:dim,:)= cos(20.0*pi*Y(2:dim,:)./sqrt(repmat((2:dim)',[1,num])));
        tmp11        = 4.0*sum(tmp1(3:2:dim,:)) - 2.0*prod(tmp2(3:2:dim,:)) + 2.0;  % odd index
        tmp21        = 4.0*sum(tmp1(2:2:dim,:)) - 2.0*prod(tmp2(2:2:dim,:)) + 2.0;  % even index
        tmp          = max(0,(1.0/N+2.0*E)*(sin(2.0*N*pi*x(1,:))-abs(2*N*Gtemp)));
        y(1,:)       = x(1,:)       + tmp + 2.0*tmp11/size(3:2:dim,2);
        y(2,:)       = 1.0 - x(1,:) + tmp + 2.0*tmp21/size(2:2:dim,2);
        clear Y tmp1 tmp2 tmp G Gtemp tmp11 tmp21;
    end
end

%% ud6
% x agre columnwise and y are rowwise, the input x must be inside the search space and
% it could be a matrxix 
% time window for which the process remains constant
% step controls the distance between 2 paretos nT
function p=ud6(p,dim)
p.name   = 'ud6';
p.pd     = dim;
p.od     = 3;
p.domain = [-2*ones(dim,1) 2*ones(dim,1)];
p.domain(1:2,1) = 0;
p.domain(1:2,2) = 1;
p.func   = @evaluate;

    function y=evaluate(x)
        global itrCounter step window;
        [dim, num]  = size(x);
        t           =(floor(itrCounter/window))/step;
        G           =sin(0.5*pi*t);
        F           =1+abs(G);
        Gtemp       =G*ones(1,num);
        Ftemp       =F*ones(1,num);
        Y           = zeros(dim,num);
        Y(3:dim,:)  = (x(3:dim,:)-Gtemp - 2.0*repmat(x(2,:),[dim-2,1]).*sin(2.0*pi*repmat(x(1,:),[dim-2,1]) + pi/dim*repmat((3:dim)',[1,num]))).^2;
        tmp1        = sum(Y(4:3:dim,:));  % j-1 = 3*k
        tmp2        = sum(Y(5:3:dim,:));  % j-2 = 3*k
        tmp3        = sum(Y(3:3:dim,:));  % j-0 = 3*k
        y(1,:)      = Ftemp.*(cos(0.5*pi*x(1,:)).*cos(0.5*pi*x(2,:)) + 2.0*tmp1/size(4:3:dim,2));
        y(2,:)      = Ftemp.*(cos(0.5*pi*x(1,:)).*sin(0.5*pi*x(2,:)) +2.0*tmp2/size(5:3:dim,2));
        y(3,:)      = Ftemp.*(sin(0.5*pi*x(1,:))                     + 2.0*tmp3/size(3:3:dim,2));
        clear Y G Gtemp Ftemp tmp1 tmp2 tmp3;
    end
end

%% ud7
% x agre columnwise and y are rowwise, the input x must be inside the search space and
% it could be a matrxix 
% time window for which the process remains constant
% step controls the distance between 2 paretos nT
function p=ud7(p,dim)
p.name   = 'ud7';
p.pd     = dim;
p.od     = 2;
p.domain = [-1*ones(dim,1) ones(dim,1)];
p.domain(1,1) = 0;
p.func   = @evaluate;

    function y=evaluate(x)
        global itrCounter step window;
        [dim, num]  = size(x);
        %     tmp         = zeros(dim,num);
        t           =(floor(itrCounter/window))/step;
        %     N=2;
        G           =sin(0.5*pi*t);
        %    Gtemp=G*ones(1,num);
        F           =ceil(dim*G);
        Ftemp       =F*ones(dim-1,num);
        %     K=0.5+abs(G);
        %     Ktemp=K*ones(1,num);
        %     H=1/N+(N-1/N)*abs(G);
        H           =0.5+abs(G);
        %     Htemp=H*ones(1,num);
        %     tmp(2:dim,:)= (x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(F+repmat((2:dim)',[1,num])))-G).^2;
        %     tmp1        = sum(tmp(3:2:dim,:));  % odd index
        %     tmp2        = sum(tmp(2:2:dim,:));  % even index
        tmp         = zeros(dim,num);
        tmp(2:dim,:)= (x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(Ftemp+repmat((2:dim)',[1,num])))).^2;
        tmp1        = sum(tmp(3:2:dim,:));  % odd index
        tmp(2:dim,:)= (x(2:dim,:)- sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(Ftemp+repmat((2:dim)',[1,num])))).^2;
        tmp2        = sum(tmp(2:2:dim,:));  % even index
        y(1,:)      = x(1,:)+ + 2.0*tmp1/size(3:2:dim,2);
        y(2,:)      = 1.0 - H.*(x(1,:)).^H + 2.0*tmp2/size(2:2:dim,2);
        clear tmp;
    end
end

%% ud8
% x agre columnwise and y are rowwise, the input x must be inside the search space and
% it could be a matrxix 
% time window for which the process remains constant
% step controls the distance between 2 paretos nT
function p=ud8(p,dim)
p.name   = 'ud8';
p.pd     = dim;
p.od     = 2;
p.domain = [-2*ones(dim,1) 2*ones(dim,1)];
p.domain(1,1) = 0;
p.domain(1,2) = 1;
p.func   = @evaluate;

    function y=evaluate(x)
        global itrCounter step window;
        [dim, num]  = size(x);
        t=(floor(itrCounter/window))/step;
        G=sin(0.5*pi*t);
        Gtemp=G*ones(1,num);
        %F=ceil(dim*G);
        %Ftemp=F*ones(dim-1,num);
        K=0.5+abs(G);
        %Ktemp=K*ones(1,num);
        %     H=1/N+(N-1/N)*abs(G);
        %H=0.5+abs(G);
        %Htemp=H*ones(1,num);
        %     tmp(2:dim,:)= (x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(F+repmat((2:dim)',[1,num])))-G).^2;
        %     tmp1        = sum(tmp(3:2:dim,:));  % odd index
        %     tmp2        = sum(tmp(2:2:dim,:));  % even index
        tmp         = zeros(dim,num);
        tmp(2:dim,:)= (x(2:dim,:)-Gtemp - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(repmat((2:dim)',[1,num])))).^2;
        tmp1        = sum(tmp(3:2:dim,:));  % odd index
        tmp(2:dim,:)= (x(2:dim,:)-Gtemp - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(repmat((2:dim)',[1,num])))).^2;
        tmp2        = sum(tmp(2:2:dim,:));  % even index
        y(1,:)      = x(1,:)+ 2.0*tmp1/size(3:2:dim,2);
        y(2,:)      = 1.0 - K.*(x(1,:)).^K + 2.0*tmp2/size(2:2:dim,2);
        clear tmp;
    end
end

%% ud9
% x agre columnwise and y are rowwise, the input x must be inside the search space and
% it could be a matrxix 
% time window for which the process remains constant
% step controls the distance between 2 paretos nT
function p=ud9(p,dim)
p.name   = 'ud9';
p.pd     = dim;
p.od     = 2;
p.domain = [-1*ones(dim,1) 2*ones(dim,1)];
p.domain(1,1) = 0;
p.domain(1,2) = 1;
p.func   = @evaluate;

    function y=evaluate(x)
        global itrCounter step window;
        [dim, num]  = size(x);
        %     tmp         = zeros(dim,num);
        t=(floor(itrCounter/window))/step;
        %     N=2;
        G=sin(0.5*pi*t);
        Gtemp=G*ones(1,num);
        %     F=ceil(dim*G);
        %     Ftemp=F*ones(dim-1,num);
        K=0.5+abs(G);
        %     Ktemp=K*ones(1,num);
        %     H=1/N+(N-1/N)*abs(G);
        %     H=0.5+abs(G);
        %     Htemp=H*ones(1,num);
        %     tmp(2:dim,:)= (x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(F+repmat((2:dim)',[1,num])))-G).^2;
        %     tmp1        = sum(tmp(3:2:dim,:));  % odd index
        %     tmp2        = sum(tmp(2:2:dim,:));  % even index
        %     tmp         = zeros(dim,num);
        %     tmp(2:dim,:)= (x(2:dim,:)-Gtemp - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(Ftemp+repmat((2:dim)',[1,num])))).^2;
        %     tmp1        = sum(tmp(3:2:dim,:));  % odd index
        %     tmp(2:dim,:)= (x(2:dim,:)-Gtemp - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(Ftemp+repmat((2:dim)',[1,num])))).^2;
        %     tmp2        = sum(tmp(2:2:dim,:));  % even index
        Y            = zeros(dim,num);
        Y(2:dim,:)   = (x(2:dim,:) - repmat(x(1,:),[dim-1,1]).^(1+Gtemp+1.5*(repmat((2:dim)',[1,num])-2.0)/(dim-2.0))-Gtemp).^2;
        tmp1         = sum(Y(3:2:dim,:));% odd index
        tmp2         = sum(Y(2:2:dim,:));% even index 
        y(1,:)      = x(1,:) + 2.0*tmp1/size(3:2:dim,2);
        y(2,:)      = 1.0 - K.*(x(1,:)).^K+ 2.0*tmp2/size(2:2:dim,2);
        clear tmp;
    end
end

%% ud10
% x agre columnwise and y are rowwise, the input x must be inside the search space and
% it could be a matrxix 
% time window for which the process remains constant
% step controls the distance between 2 paretos nT
function p=ud10(p,dim)
p.name   = 'ud10';
p.pd     = dim;
p.od     = 2;
p.domain = [-1*ones(dim,1) ones(dim,1)];
p.domain(1,1) = 0;
p.func   = @evaluate;

    function y=evaluate(x)
        global itrCounter step window;
        [dim, num]  = size(x);
        t=(floor(itrCounter/window))/step;
        G=sin(0.5*pi*t);
        K=0.5+abs(G);
        Gtemp=G*ones(1,num);
        N           = 10.0;
        E           = 0.1;
        [dim, num]  = size(x);
        Y           = zeros(dim,num);
        Y(2:dim,:)  = x(2:dim,:)- sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
        H           = zeros(dim,num);
        H(2:dim,:)  = 2.0*Y(2:dim,:).^2 - cos(4.0*pi*Y(2:dim,:)) + 1.0;
        tmp1        = sum(H(3:2:dim,:));  % odd index
        tmp2        = sum(H(2:2:dim,:));  % even index
        tmp         = (0.5/N+E)*abs(sin(2.0*N*pi*x(1,:))-abs(2*N*Gtemp));
        y(1,:)      = x(1,:)      + tmp + 2.0*tmp1/size(3:2:dim,2);
        y(2,:)      = 1- x(1,:).*K+ tmp + 2.0*tmp2/size(2:2:dim,2);
    end
end

%% ud11
% x agre columnwise and y are rowwise, the input x must be inside the search space and
% it could be a matrxix 
% time window for which the process remains constant
% step controls the distance between 2 paretos nT
function p=ud11(p,dim)
p.name   = 'ud11';
p.pd     = dim;
p.od     = 2;
p.domain = [-1*ones(dim,1) ones(dim,1)];
p.domain(1,1) = 0;
p.func   = @evaluate;

    function y=evaluate(x)
        global itrCounter step window;
        [dim, num]   = size(x);
        t=(floor(itrCounter/window))/step;
        G=sin(0.5*pi*t);
        K=0.5+abs(G);
        Gtemp=G*ones(1,num);
        N            = 2.0;
        E            = 0.1;
        Y            = zeros(dim,num);
        Y(2:dim,:)  = x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
        tmp1         = zeros(dim,num);
        tmp1(2:dim,:)= Y(2:dim,:).^2;
        tmp2         = zeros(dim,num);
        tmp2(2:dim,:)= cos(20.0*pi*Y(2:dim,:)./sqrt(repmat((2:dim)',[1,num])));
        tmp11        = 4.0*sum(tmp1(3:2:dim,:)) - 2.0*prod(tmp2(3:2:dim,:)) + 2.0;  % odd index
        tmp21        = 4.0*sum(tmp1(2:2:dim,:)) - 2.0*prod(tmp2(2:2:dim,:)) + 2.0;  % even index
        tmp          = max(0,(1.0/N+2.0*E)*(sin(2.0*N*pi*x(1,:))-abs(2*N*Gtemp)));
        y(1,:)       = x(1,:)       + tmp + 2.0*tmp11/size(3:2:dim,2);
        y(2,:)       = 1.0 - x(1,:).*K + tmp + 2.0*tmp21/size(2:2:dim,2);
        clear Y tmp1 tmp2;
    end
end

%% ud12
% x agre columnwise and y are rowwise, the input x must be inside the search space and
% it could be a matrxix 
% time window for which the process remains constant
% step controls the distance between 2 paretos nT
function p=ud12(p,dim)
p.name   = 'ud12';
p.pd     = dim;
p.od     = 3;
p.domain = [-2*ones(dim,1) 2*ones(dim,1)];
p.domain(1:2,1) = 0;
p.domain(1:2,2) = 1;
p.func   = @evaluate;

    function y=evaluate(x)
        global itrCounter step window;
        [dim, num]  = size(x);
        t           =(floor(itrCounter/window))/step;
        G           =sin(0.5*pi*t);
        F           =1+abs(G);
        Gtemp       =G*ones(1,num);
        Ftemp       =F*ones(1,num);
        Y           = zeros(dim,num);
        Y(3:dim,:)  = (x(3:dim,:)-Gtemp - 2.0*repmat(x(2,:),[dim-2,1]).*sin(2.0*pi*repmat(x(1,:),[dim-2,1]) + pi/dim*repmat((3:dim)',[1,num]))).^2;
        tmp1        = sum(Y(4:3:dim,:));  % j-1 = 3*k
        tmp2        = sum(Y(5:3:dim,:));  % j-2 = 3*k
        tmp3        = sum(Y(3:3:dim,:));  % j-0 = 3*k
        y(1,:)      = Ftemp.*(cos(0.5*pi*x(1,:)).*cos(0.5*pi*x(2,:)) +abs(Gtemp)+ 2.0*tmp1/size(4:3:dim,2));
        y(2,:)      = Ftemp.*(cos(0.5*pi*x(1,:)).*sin(0.5*pi*x(2,:)) + abs(Gtemp)+2.0*tmp2/size(5:3:dim,2));
        y(3,:)      = Ftemp.*(sin(0.5*pi*x(1,:))                     +abs(Gtemp)+ 2.0*tmp3/size(3:3:dim,2));
        clear Y;
    end
end

%% ud13
% x agre columnwise and y are rowwise, the input x must be inside the search space and
% it could be a matrxix 
% time window for which the process remains constant
% step controls the distance between 2 paretos nT
function p=ud13(p,dim)
p.name   = 'ud13';
p.pd     = dim;
p.od     = 2;
p.domain = [-1*ones(dim,1) ones(dim,1)];
p.domain(1,1) = 0;
p.func   = @evaluate;

    function y=evaluate(x)
        global itrCounter step window;
        [dim, num]  = size(x);
        if mod(iter,window1)==0 
            Hpf;
        end
        tmp         = zeros(dim,num);
        tmp(2:dim,:)= (x(2:dim,:) -Gps- sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(Fps+repmat((2:dim)',[1,num])))).^2;
        tmp1        = sum(tmp(3:2:dim,:));  % odd index
        tmp(2:dim,:)= (x(2:dim,:)- Gps-sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(Fps+repmat((2:dim)',[1,num])))).^2;
        tmp2        = sum(tmp(2:2:dim,:));  % even index
        y(1,:)      = x(1,:)+abs(Gpf) + 2.0*tmp1/size(3:2:dim,2);
        y(2,:)      = 1.0 - Mpf*(x(1,:)).^Hpf+abs(Gpf) + 2.0*tmp2/size(2:2:dim,2);
        clear tmp;
    end
end

%% ud14
% x agre columnwise and y are rowwise, the input x must be inside the search space and
% it could be a matrxix 
% time window for which the process remains constant
% step controls the distance between 2 paretos nT
function p=ud14(p,dim)
p.name   = 'ud14';
p.pd     = dim;
p.od     = 2;
p.domain = [-2*ones(dim,1) 2*ones(dim,1)];
p.domain(1,1) = 0;
p.domain(1,2) = 1;
p.func   = @evaluate;

    function y=evaluate(x)
        global itrCounter step window;
        [dim, num]  = size(x);
        if mod(iter,window1)==0
            Hpf;
        end
        Y           = zeros(dim,num);
        Y(2:dim,:)  = (x(2:dim,:) - repmat(x(1,:),[dim-1,1]).^(1+Gps+1.5*(repmat((2:dim)',[1,num])-2.0)/(dim-2.0))-Gps).^2;
        tmp1        = sum(Y(3:2:dim,:));% odd index
        tmp2        = sum(Y(2:2:dim,:));% even index 
        y(1,:)      = x(1,:)+abs(Gpf) + 2.0*tmp1/size(3:2:dim,2);
        y(2,:)      = 1.0 - Mpf*(x(1,:)).^Hpf+abs(Gpf) + 2.0*tmp2/size(2:2:dim,2);
        clear tmp;
    end
end

%% ---------CEC2014 Benchmark with constraintrs. Ref.[]----------
%% cd1
% x agre columnwise and y are rowwise, the input x must be inside the search space and
% it could be a matrxix 
% time window for which the process remains constant
% step controls the distance between 2 paretos nT
% c: constrains value
function p=cd1(p,dim)
p.name   = 'cd1';
p.pd     = dim;
p.od     = 2;
p.domain = [zeros(dim,1) ones(dim,1)];
p.func   = @evaluate;

    function [y,c]=evaluate(x)
        global itrCounter step window;
        [dim, num]   = size(x);
        t=(floor(itrCounter/window))/step;
        G=sin(0.5*3.14*t);
        Gtemp=G*ones(dim-1,num);
        Ftemp=1.5*ones(1,num)+G*ones(1,num);
        a            = 1.0;
        N            = 10.0;
        [dim, num]   = size(x);
        Y            = zeros(dim,num);
        Y(2:dim,:)   = (x(2:dim,:) - repmat(x(1,:),[dim-1,1]).^(0.5+1.5*(repmat((2:dim)',[1,num])-2.0)/(dim-2.0))-Gtemp).^2;
        tmp1         = sum(Y(3:2:dim,:));% odd index
        tmp2         = sum(Y(2:2:dim,:));% even index 
        y(1,:)       = x(1,:)       + 2.0*tmp1/size(3:2:dim,2);
        y(2,:)       = 1.0 - x(1,:).^Ftemp + 2.0*tmp2/size(2:2:dim,2);
        c(1,:)       = y(1,:).^Ftemp + y(2,:) - a*abs(sin(N*pi*(y(1,:).^Ftemp-y(2,:)+1.0))) - 1.0;
    end
end

%% cd2
% x agre columnwise and y are rowwise, the input x must be inside the search space and
% it could be a matrxix 
% time window for which the process remains constant
% step controls the distance between 2 paretos nT
% c: constrains value
function p=cd2(p,dim)
p.name   = 'cd2';
p.pd     = dim;
p.od     = 2;
p.domain = [zeros(dim,1) ones(dim,1)];
p.func   = @evaluate;

    function [y,c]=evaluate(x)
        global itrCounter step window;
        [dim, num]   = size(x);
        t=(floor(itrCounter/window))/step;
        G=sin(0.5*3.14*t);
        Gtemp=G*ones(dim-1,num);
        Ftemp=1.5*ones(1,num)+G*ones(1,num);
        a            = 1.0;
        N            = 10.0;
        [dim, num]   = size(x);
        Y            = zeros(dim,num);
        Y(2:dim,:)   = (x(2:dim,:) - repmat(x(1,:),[dim-1,1]).^(0.5+1.5*(repmat((2:dim)',[1,num])-2.0)/(dim-2.0))-Gtemp).^2;
        tmp1         = sum(Y(3:2:dim,:));% odd index
        tmp2         = sum(Y(2:2:dim,:));% even index 
        y(1,:)       = x(1,:)       + 2.0*tmp1/size(3:2:dim,2);
        y(2,:)       = 1.0+abs(Gtemp(1,:)) - x(1,:) + 2.0*tmp2/size(2:2:dim,2);
        c(1,:)       = y(1,:) + y(2,:) - a*abs(sin(N*pi*(y(1,:).^Ftemp-y(2,:)+1.0))) - 1.0-abs(Gtemp(1,:));
        f=[y(1,:);y(2,:)];
    end
end


%% ------------------Transformation functions------------------
% ---------Warning: below are not test problem ----------------
% -------------------------------------------------------------

% WFG benchmark initialize
function [noSols, n, S, D, A, Y] = wfg_initialize(Z, M, k, l, testNo)

    Z = Z';
    % Check for correct number of inputs.
    if nargin ~= 5
        error('Five inputs are required.');
    end

    % Get total number of decision variables and no. of candidates.
    [noSols, n] = size(Z);

    % Data input checks.
    if n ~= (k + l)
        error('Inconsistent number of variabes.');
    end
    if rem(k,M-1) ~= 0
        error('k must be divisible by M-1.');
    end
    if (testNo == 2 || testNo == 3) && (rem(l,2) ~= 0)
        error('l must be a multiple of 2 for WFG2 and WFG3.');
    end

    % Could also check data input range z_i in [0, 2i].

    % Initialise function-wide constants.
    NO_TESTS = 10;
    S = NaN * ones(1, M);
    for i = 1:M
        S(i) = 2*i;
    end
    D = 1;
    A = ones(NO_TESTS, M-1);
    A(3,2:M-1) = 0;

    % Transform all variable ranges to [0 1].
    x = Z;
    for i = 1:n
        Y(:,i) = Z(:,i) ./ (2*i);
    end
end

% Reduction: weighted sum.
function ybar = r_sum(y, weights)

[noSols noY] = size(y);
wgtMatrix=rep(weights,[noSols 1]);
ybar = y .* wgtMatrix;
ybar = sum(ybar, 2) ./ sum(wgtMatrix, 2);

end

% Reduction: non-separable.
function y_bar = r_nonsep(y, A)

[noSols noY] = size(y);

y_bar = 0;
for j = 1:noY
    innerSum = 0;
    for k = 0:(A-2)
        innerSum = innerSum + abs(y(:,j) - y(:,1+mod(j+k,noY)));
    end
    y_bar = y_bar + y(:,j) + innerSum;
end
y_bar = y_bar / ( (noY/A) * ceil(A/2) * (1+2*A-2*ceil(A/2)) );

end

% Bias: polynomial.
function y_bar = b_poly(y, alpha)

y_bar = y.^alpha;

end

% Bias: flat region
function y_bar = b_flat(y, A, B, C)

[noSols noY] = size(y);
min1 = min(0, floor(y - B));
min2 = min(0, floor(C - y));
y_bar = A + min1*A.*(B-y)/B - min2*(1-A).*(y-C)/(1-C);

% Machine precision problems can cause y_bar to go slightly negative so
% force >=0 condition.
y_bar=max(0,y_bar);

end

% Bias: parameter dependent.
function ybar = b_param(y, uy, A, B, C)

[noSols noY] = size(y);
v = A - (1 - 2*uy) .* abs(floor(0.5 - uy) + A);
v = rep(v, [1 noY]);
ybar = y.^(B + (C-B)*v);

end

% Shift: linear.
function ybar = s_linear(y, A)

ybar = abs(y - A) ./ abs(floor(A - y) + A);

end

% Shift: deceptive.
function ybar = s_decep(y, A, B, C)

y1 = floor(y - A + B) * (1 - C + (A - B)/B) / (A - B);
y2 = floor(A + B - y) * (1 - C + (1 - A - B)/B) / (1 - A - B);
ybar = 1 + (abs(y - A) - B) .* (y1 + y2 + 1/B);

end

% Shift: multi-modal.
function ybar = s_multi(y, A, B, C)

y1 = abs(y-C) ./ (2*(floor(C-y)+C));
ybar = (1 + cos((4*A+2)*pi*(0.5 - y1)) + 4*B*y1.^2) / (B+2);

end

% Shape functions.
% Linear.
function f = h_linear(x)

[noSols mMinusOne] = size(x);

M = mMinusOne + 1;
f = NaN * ones(noSols, M);

f(:,1) = prod(x, 2);
for i = 2:mMinusOne
    f(:,i) = prod(x(:,1:M-i), 2) .* (1 - x(:,M-i+1));
end
f(:,M) = 1 - x(:,1);
end

% Convex.
function f = h_convex(x)

[noSols mMinusOne] = size(x);

M = mMinusOne + 1;
f = NaN * ones(noSols, M);

f(:,1) = prod(1-cos(x*pi/2), 2);
for i = 2:mMinusOne
    f(:,i) = prod(1-cos(x(:,1:M-i)*pi/2), 2) .* (1-sin(x(:,M-i+1)*pi/2));
end
f(:,M) = 1-sin(x(:,1)*pi/2);

end

% Concave.
function f = h_concave(x)

[noSols mMinusOne] = size(x);

M = mMinusOne + 1;
f = NaN * ones(noSols, M);

f(:,1) = prod(sin(x*pi/2), 2);
for i = 2:mMinusOne
    f(:,i) = prod(sin(x(:,1:M-i)*pi/2), 2) .* cos(x(:,M-i+1)*pi/2);
end
f(:,M) = cos(x(:,1)*pi/2);

end

% Mixed.
function f = h_mixed(x,alpha,A)

f = (1 - x(:,1) - cos(2*A*pi*x(:,1) + pi/2) / (2*A*pi)).^alpha;

end

% Disconnected.
function f = h_disc(x,alpha,beta,A)

f = 1 - x(:,1).^alpha .* cos(A * x(:,1).^beta * pi).^2;

end

function MatOut = rep(MatIn,REPN)

% Get size of input matrix
   [N_D,N_L] = size(MatIn);

% Calculate
   Ind_D = rem(0:REPN(1)*N_D-1,N_D) + 1;
   Ind_L = rem(0:REPN(2)*N_L-1,N_L) + 1;

% Create output matrix
   MatOut = MatIn(Ind_D,Ind_L);
end