%Group Teaching Optimization Algorithm with Information sharing (ISGTOA) Source Code %

function [BestCost,BestValue]=ISGTOA(fhd,nPop,nVar,VarMin,VarMax,MaxIt,X)
%fhd--------objective function
%nPop-------population size
%nVar-------dimension
%VarMin-----lower limits of variables
%VarMax-----upper limits of variables
%MaxIt------the maximum number of function evaluations
%X----------initilization population
%BestCost---convergence curve
%BestValue--the optimal solution
FES=0; % FES---the current number of function evaluations
for i=1:nPop
    cost(i)=fhd(X(i,:));
    FES=FES+1;
end
[~,lb]=sort(cost);%sort for population
for i=1:nPop
    if i<=nPop/2
        X1(i,:)=X(lb(i),:);% X1-----outstanding group
        cost1(i)=cost(lb(i));
    else
        X2(i-nPop/2,:)=X(lb(i),:);% X2-----common group
        cost2(i-nPop/2)=cost(lb(i));
    end
end
TEACHER=X(lb(1),:);%the current optimal solution
COST=cost(lb(1));
t=0;
NUM=2;
POT(1).Pos=TEACHER;% initilization archive B
POT(1).Cost=COST;
POT(2).Pos=X(lb(2),:);
POT(2).Cost=cost(lb(2));
while FES<MaxIt
    t=t+1;
    [X1,cost1,POT,NUM]=TLBO1(fhd,nPop/2,nVar,VarMin,VarMax,X1,cost1,TEACHER,COST,POT,NUM,X);% optimize outstanding group 
    [X2,cost2,POT,NUM]=TLBO2(fhd,nPop/2,nVar,VarMin,VarMax,X2,cost2,TEACHER,COST,POT,NUM,X);%optimize common group
    FES=FES+2*nPop;%update the  current number of function evaluations
    for i=1:nPop
        if i<=nPop/2
            X(i,:)=X1(i,:);
            cost(i)=cost1(i);
        else
            X(i,:)=X2(i-nPop/2,:);
            cost(i)=cost2(i-nPop/2);
        end
    end
    [~,lb]=sort(cost);%sort for population
    for i=1:nPop
        if i<=nPop/2
            X1(i,:)=X(lb(i),:);
            cost1(i)=cost(lb(i));
        else
            X2(i-nPop/2,:)=X(lb(i),:);
            cost2(i-nPop/2)=cost(lb(i));
        end
    end
    if cost(lb(1))<COST
        TEACHER=X(lb(1),:);
        COST=cost(lb(1));
    end
    MEAN=(X(lb(1),:)+X(lb(2),:)+X(lb(3),:))/3;%teacher allocation phase
    MCOST =fhd(MEAN);
    FES=FES+1;
    if MCOST<COST
        TEACHER=MEAN;
        COST=MCOST;
    end
    BestCost(t)=COST;
    BestValue=COST;
    
end

end

function [pop,cost,POT,NUM]=TLBO1(fhd,nPop,nVar,VarMin,VarMax,pop,cost,TEACHER,COST,POT,NUM,GOT)
VarSize = [1 nVar]; 
popx=pop;
costx=cost;
GOT=GOT(randperm(nPop*2),:); %generate archive A
for it=1:1
    Mean = 0;
    for i=1:nPop
        Mean = Mean + pop(i,:);
    end
    Mean = Mean/nPop;
    Teacher = TEACHER;
    BEST=COST;
    XBEST=Teacher;
    M=mean(cost);
    for i=1:nPop
        ki=round(1+rand);
        fi=rand;
        gi=1-fi;
        if cost(i)>M
            newsol= pop(i,:)+ rand(1,nVar).*(Teacher-ki*(fi*Mean+gi*pop(i,:)));
        else
            newsol= pop(i,:)+ceil(1+abs(randn)).*randn.*(GOT(i,:)-pop(i,:));
        end
        newsol = max( newsol, VarMin);
        newsol = min( newsol, VarMax);
        newsol_cost =fhd(newsol);
        if  newsol_cost<cost(i)
            pop(i,:) = newsol;
            cost(i)= newsol_cost;
            NUM=NUM+1;
            if NUM<=nPop*2
                POT(NUM).Cost=cost(i);
                POT(NUM).Pos=pop(i,:) ;
            else
                a=randperm(nPop*2,1);
                b=randperm(nPop*2,1);
                while a==b
                    a=randperm(nPop*2,1);
                    b=randperm(nPop*2,1);
                end
                if POT(b).Cost<POT(a).Cost
                    a=b;
                end
                POT(a).Cost=cost(i);
                POT(a).Pos=pop(i,:) ;
            end
            if newsol_cost < BEST
                BEST=newsol_cost;
                XBEST=newsol;
            end
        end
        
    end
    for i=1:nPop
        A = 1:nPop;
        A(i)=[];
        j = A(randi(nPop-1));
        
        Step = pop(i,:) - pop(j,:);
        if cost(j) < cost(i)
            Step = -Step;
        end
        a=randperm(length(POT),1);
        b=randperm(length(POT),1);
        if POT(b).Cost<POT(a).Cost
            a=b;
        end
        if POT(a).Cost<cost(i)
            newsol = pop(i,:) +rand(VarSize).*Step+rand(VarSize).*(POT(a).Pos-popx(i,:));
        else
            newsol = pop(i,:) +rand(VarSize).*Step-rand(VarSize).*(POT(a).Pos-popx(i,:));
        end   
        newsol = max( newsol, VarMin);
        newsol = min( newsol, VarMax);
        newsol_cost =fhd(newsol);
        if newsol_cost<cost(i)
            pop(i,:) = newsol;
            cost(i)=newsol_cost;
            if NUM<=nPop*2
                POT(NUM).Cost=cost(i);
                POT(NUM).Pos=pop(i,:) ;
            else
                a=randperm(nPop*2,1);
                b=randperm(nPop*2,1);
                while a==b
                    a=randperm(nPop*2,1);
                    b=randperm(nPop*2,1);
                end
                if POT(b).Cost<POT(a).Cost
                    a=b;
                end
                POT(a).Cost=cost(i);
                POT(a).Pos=pop(i,:) ;
            end
            if newsol_cost < BEST
                BEST = newsol_cost;
                XBEST=newsol;
            end
        end
    end
end

end

function [pop,cost,POT,NUM]=TLBO2(fhd,nPop,nVar,VarMin,VarMax,pop,cost,TEACHER,COST,POT,NUM,GOT)
popx=pop;
GOT=GOT(randperm(nPop*2),:);%generate archive A
VarSize = [1 nVar]; 
for it=1:1
    Teacher = TEACHER;
    BEST=COST;
    M=mean(cost);
    for i=1:nPop
        if cost(i)>M
            newsol= pop(i,:)+ 2.*rand(1,nVar).*(Teacher-pop(i,:));
        else
            newsol= pop(i,:)+ceil(1+abs(randn)).*randn.*(GOT(i,:)-pop(i,:));
        end
        newsol = max( newsol, VarMin);
        newsol = min( newsol, VarMax);
        newsol_cost =fhd(newsol);
        if  newsol_cost<cost(i)
            pop(i,:) = newsol;
            cost(i)= newsol_cost;
            if NUM<=nPop*2
                POT(NUM).Cost=cost(i);
                POT(NUM).Pos=pop(i,:) ;
            else
                a=randperm(nPop*2,1);
                b=randperm(nPop*2,1);
                while a==b
                    a=randperm(nPop*2,1);
                    b=randperm(nPop*2,1);
                end
                if POT(b).Cost<POT(a).Cost
                    a=b;
                end
                POT(a).Cost=cost(i);
                POT(a).Pos=pop(i,:) ;
            end
            if newsol_cost < BEST
                BEST=newsol_cost;
                XBEST=newsol;
            end
        end
    end
    for i=1:nPop
        A = 1:nPop;
        A(i)=[];
        j = A(randi(nPop-1));
        Step = pop(i,:) - pop(j,:);
        if cost(j) < cost(i)
            Step = -Step;
        end
        a=randperm(length(POT),1);
        b=randperm(length(POT),1);
        if POT(b).Cost<POT(a).Cost
            a=b;
        end
        if POT(a).Cost<cost(i)
            newsol = pop(i,:) +rand(VarSize).*Step+rand(VarSize).*(POT(a).Pos-popx(i,:));
        else
            newsol = pop(i,:) +rand(VarSize).*Step-rand(VarSize).*(POT(a).Pos-popx(i,:));
        end
        newsol = max( newsol, VarMin);
        newsol = min( newsol, VarMax);
        newsol_cost =fhd(newsol);
        if newsol_cost<cost(i)
            pop(i,:) = newsol;
            cost(i)=newsol_cost;
            if NUM<=nPop*2
                POT(NUM).Cost=cost(i);
                POT(NUM).Pos=pop(i,:) ;
            else
                a=randperm(nPop*2,1);
                b=randperm(nPop*2,1);
                while a==b
                    a=randperm(nPop*2,1);
                    b=randperm(nPop*2,1);
                end
                if POT(b).Cost<POT(a).Cost
                    a=b;
                end
                POT(a).Cost=cost(i);
                POT(a).Pos=pop(i,:) ;
            end
            if newsol_cost < BEST
                BEST = newsol_cost;
                XBEST=newsol;
            end
        end
    end
end

end




