clc
clear all
fhd=@sphere;
nPop=50;
nVar=50;
VarMin=-100.*ones(1,nVar);
VarMax=100.*ones(1,nVar);
MaxIt=50000;
for i=1:nPop
    X(i,:)=VarMin+(VarMax-VarMin).*rand(1,nVar);
end
[BestCost,BestValue]=ISGTOA(fhd,nPop,nVar,VarMin,VarMax,MaxIt,X);
plot(1:495,BestCost,'r')
title('Objective space')
xlabel('Iteration');
ylabel('Best score obtained so far');