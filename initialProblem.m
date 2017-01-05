function initialProblem()
    global CostFunction nVar VarMin VarMax numOfObj TestProblem dynamic;
    switch TestProblem
        case 1
            CostFunction=@(x) MyCost1(x); 
            nVar=5;
            VarMin=-4*ones(1,nVar);
            VarMax=4*ones(1,nVar);
            numOfObj = 2;
        case 2
            CostFunction=@(x) MyCost2(x);
            nVar=3;
            VarMin=-5*ones(1,nVar);
            VarMax=5*ones(1,nVar);
            numOfObj = 2;
        case 3
            CostFunction=@(x) MyCost3(x); 
            nVar=2;
            VarMin=0;
            VarMax=1;
            numOfObj = 2;
        case 4
            global k M;
            k = 5;
            M = 3;
            mop = testmop('DTLZ1',7);
            CostFunction=@(x) mop.func(x);
            nVar = 7;
            VarMin = mop.domain(:,1)';
            VarMax = mop.domain(:,2)';
            numOfObj = 3;
            dynamic = 0;
        case 5
            CostFunction =@(x) MyCost5(x);%zdt6
            nVar = 30;
            numOfObj = 2;
            VarMin=zeros(1,nVar);
            VarMax=ones(1,nVar);
        case 6
            CostFunction = @(x) MyCost6(x);%zdt3
            nVar = 30;
            VarMin=zeros(1,nVar);
            VarMax=ones(1,nVar);
            numOfObj = 2;
        case 7
            CostFunction=@(x) MyCost7(x); %zdt1
            nVar=30;
            VarMin=zeros(1,nVar);
            VarMax=ones(1,nVar);
            numOfObj = 2;

        case 8
            CostFunction=@(x) MyCost8(x); %zdt2
            nVar=30;
            VarMin=zeros(1,nVar);
            VarMax=ones(1,nVar);
            numOfObj = 2;
        case 9
            CostFunction=@(x) MyCost9(x); %zdt4
            nVar=10;
            VarMin=-5*ones(1,nVar);
            VarMax=5*ones(1,nVar);
            VarMin(1) = 0;
            VarMax(1) = 1;
            numOfObj = 2;
        case 10
            global k M;
            k = 10;
            M = 3;
            mop = testmop('DTLZ2',12);
            CostFunction=@(x) mop.func(x);
            nVar = 12;
            VarMin = mop.domain(:,1)';
            VarMax = mop.domain(:,2)';
            numOfObj = 3;
            dynamic = 0;
        case 12
            CostFunction=@(x) MyCost12(x); %cec09-1
            nVar=30;
            VarMin=-1*ones(1,nVar);
            VarMax=ones(1,nVar);
            VarMin(1) = 0;
            VarMax(1) = 1;
            numOfObj = 2;
        case 13
            CostFunction=@(x) MyCost13(x); %cec09-2
            nVar=30;
            VarMin=-1*ones(1,nVar);
            VarMax=ones(1,nVar);
            VarMin(1) = 0;
            VarMax(1) = 1;
            numOfObj = 2;
         case 14
            CostFunction=@(x) MyCost14(x); %cec09-3
            nVar=30;
            VarMin=zeros(1,nVar);
            VarMax=ones(1,nVar);
            numOfObj = 2;
         case 15
            CostFunction=@(x) MyCost15(x); %cec09-7
            nVar=30;
            VarMin=-1*ones(1,nVar);
            VarMax=ones(1,nVar);
            VarMin(1) = 0;
            VarMax(1) = 1;
            numOfObj = 2;
        case 16
            CostFunction=@(x) MyCost16(x); %cec09-4
        case 17
            CostFunction=@(x) MyCost17(x); %cec09-5
        case 18
            CostFunction=@(x) MyCost18(x); %cec09-6
        case 19
            CostFunction=@(x) MyCost19(x); %cec09-8
        case 20
            CostFunction=@(x) MyCost20(x); %cec09-9
        case 21
            CostFunction=@(x) MyCost21(x); %cec09-10
        case 22
            global k l M;
            k = 2;
            l = 4;
            M = 2;
            mop = testmop('wfg1',6);
            CostFunction=@(x) mop.func(x);
            nVar = 6;
            VarMin = mop.domain(:,1)';
            VarMax = mop.domain(:,2)';
            numOfObj = 2;
            dynamic = 0;
        case 23
            global k l M;
            k = 2;
            l = 4;
            M = 2;
            mop = testmop('wfg2',6);
            CostFunction=@(x) mop.func(x);
            nVar = 6;
            VarMin = mop.domain(:,1)';
            VarMax = mop.domain(:,2)';
            numOfObj = 2;
            dynamic = 0;
        case 24
            CostFunction=@(x) wfg(x, 2, 2, 4, 3); %wfg3
            nVar = 6;
            VarMin=zeros(1,nVar);
            VarMax=ones(1,nVar);
            numOfObj = 2;
            for i = 1:nVar
                VarMax(i) = VarMax(i).*(i*2);
            end
            clear i;
        case 25
            CostFunction=@(x) wfg(x, 2, 2, 4, 4); %wfg4
            nVar = 6;
            VarMin=zeros(1,nVar);
            VarMax=ones(1,nVar);
            numOfObj = 2;
            for i = 1:nVar
                VarMax(i) = VarMax(i).*(i*2);
            end
            clear i;
        case 26
            CostFunction=@(x) wfg(x, 2, 2, 4, 5); %wfg5
            nVar = 6;
            VarMin=zeros(1,nVar);
            VarMax=ones(1,nVar);
            numOfObj = 2;
            for i = 1:nVar
                VarMax(i) = VarMax(i).*(i*2);
            end
            clear i;
        case 27
            CostFunction=@(x) wfg(x, 2, 2, 4, 6); %wfg6
            nVar = 6;
            VarMin=zeros(1,nVar);
            VarMax=ones(1,nVar);
            numOfObj = 2;
            for i = 1:nVar
                VarMax(i) = VarMax(i).*(i*2);
            end
            clear i;
        case 28
            CostFunction=@(x) wfg(x, 2, 2, 4, 7); %wfg7
            nVar = 6;
            VarMin=zeros(1,nVar);
            VarMax=ones(1,nVar);
            numOfObj = 2;
            for i = 1:nVar
                VarMax(i) = VarMax(i).*(i*2);
            end
            clear i;
        case 29
            CostFunction=@(x) wfg(x, 2, 2, 4, 8); %wfg8
            nVar = 6;
            VarMin=zeros(1,nVar);
            VarMax=ones(1,nVar);
            numOfObj = 2;
            for i = 1:nVar
                VarMax(i) = VarMax(i).*(i*2);
            end
            clear i;
        case 30
            CostFunction=@(x) wfg(x, 2, 2, 4, 9); %wfg9
            nVar = 6;
            VarMin=zeros(1,nVar);
            VarMax=ones(1,nVar);
            numOfObj = 2;
            for i = 1:nVar
                VarMax(i) = VarMax(i).*(i*2);
            end
            clear i;
        case 31
            mop = testmop('fda1',20);
            CostFunction=@(x) mop.func(x);
            nVar = 20;
            VarMin = mop.domain(:,1)';
            VarMax = mop.domain(:,2)';
            numOfObj = 2;
            dynamic = 1;
        case 32
            mop = testmop('fda2',31);
            CostFunction=@(x) mop.func(x);
            nVar = 31;
            VarMin = mop.domain(:,1)';
            VarMax = mop.domain(:,2)';
            numOfObj = 2;
            dynamic = 1;
         case 33
            mop = testmop('fda3',30);
            CostFunction=@(x) mop.func(x);
            nVar = 30;
            VarMin = mop.domain(:,1)';
            VarMax = mop.domain(:,2)';
            numOfObj = 2;
            dynamic = 1;
        case 34
            mop = testmop('fda4',12);
            CostFunction=@(x) mop.func(x);
            nVar = 12;
            VarMin = mop.domain(:,1)';
            VarMax = mop.domain(:,2)';
            numOfObj = 3;
            dynamic = 1;
        case 35
            global k M;
            k = 10;
            M = 3;
            mop = testmop('DTLZ5',12);
            CostFunction=@(x) mop.func(x);
            nVar = 12;
            VarMin = mop.domain(:,1)';
            VarMax = mop.domain(:,2)';
            numOfObj = 3;
            dynamic = 0;
        case 36
            global k M;
            k = 20;
            M = 3;
            mop = testmop('DTLZ7',22);
            CostFunction=@(x) mop.func(x);
            nVar = 22;
            VarMin = mop.domain(:,1)';
            VarMax = mop.domain(:,2)';
            numOfObj = 3;
            dynamic = 0;
        case 37
            mop = testmop('fda2_nsga2',13);
            CostFunction=@(x) mop.func(x);
            nVar = 13;
            VarMin = mop.domain(:,1)';
            VarMax = mop.domain(:,2)';
            numOfObj = 2;
            dynamic = 1;
    end
end