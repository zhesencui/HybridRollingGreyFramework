%_________________________________________________________________________%
%  RGM-MOA framework with artificial bee colony algorithm (RGM-ABC) source codes             %
%                                                                         %
%  Developed in MATLAB                                      %
%                                                                         %
%  Author and programmer: Zhesen Cui                              %
%                                                                         %
%         E-Mail: cuizhesen@gmail.com
%_________________________________________________________________________%
clc;
clear;
clear all;
close all;
load raw_data
%Problem Definition
CostFunction=@(x) Sphere(x);        % Cost Function
nVar=1;             % Number of Decision Variables
VarSize=[1 nVar];   % Decision Variables Matrix Size
VarMin=0;         % Decision Variables Lower Bound
VarMax=1;         % Decision Variables Upper Bound
%ABC Settings
MaxIt=100;              % Maximum Number of Iterations
%MaxIt
nPop=50;               % Population Size (Colony Size)
nOnlooker=nPop;         % Number of Onlooker Bees
L=round(0.6*nVar*nPop); % Abandonment Limit Parameter (Trial Limit)
a=1;                    % Acceleration Coefficient Upper Bound
[M,N] = size(raw_data);
result_rgm=zeros(1,2);
PRE_all=zeros(6,1);

for sample=1:M
    input=raw_data(sample,1:5);
    output=raw_data(sample,6);
    % Empty Bee Structure
    empty_bee.Position=[];
    empty_bee.Cost=[];
    % Initialize Population Array
    pop=repmat(empty_bee,nPop,1);
    % Initialize Best Solution Ever Found
    BestSol.Cost=inf;
    % Create Initial Population
    x=initialization(nPop,1,1,0);
    for i=1:nPop
        pop(i).Position=x(1,i);
        [PRE MRE pop(i).Position] = RollingGM(pop(i).Position,input,output);
        pop(i).Cost=MRE;
        if pop(i).Cost<=BestSol.Cost
            BestSol=pop(i);
            BestSol.Cost=pop(i).Cost;
        end
    end
    % Abandonment Counter
    C=zeros(nPop,1);
    % Array to Hold Best Cost Values
    
    BestCost=zeros(MaxIt,1);
    %% ABC Main Loop
    %MaxIt
    for it=1:MaxIt
        % Recruited Bees
        for i=1:nPop
            % Choose k randomly, not equal to i
            K=[1:i-1 i+1:nPop];
            k=K(randi([1 numel(K)]));
            % Define Acceleration Coeff.
            phi=a*unifrnd(-1,+1,VarSize);
            % New Bee Position
            newbee.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
            % Evaluation
            %newbee.Cost=CostFunction(newbee.Position);
            [PRE MRE] =RollingGM(newbee.Position,input,output);
            newbee.Cost = MRE;
            % Comparision
            if newbee.Cost<=pop(i).Cost
                pop(i)=newbee;
            else
                C(i)=C(i)+1;
            end
        end
        
        % Calculate Fitness Values and Selection Probabilities
        F=zeros(nPop,1);
        MeanCost = mean([pop.Cost]);
        for i=1:nPop
            F(i) = exp(-pop(i).Cost/MeanCost); % Convert Cost to Fitness
        end
        P=F/sum(F);
        
        % Onlooker Bees
        for m=1:nOnlooker
            
            % Select Source Site
            i=RouletteWheelSelection(P);
            
            % Choose k randomly, not equal to i
            K=[1:i-1 i+1:nPop];
            k=K(randi([1 numel(K)]));
            
            % Define Acceleration Coeff.
            phi=a*unifrnd(-1,+1,VarSize);
            
            % New Bee Position
            newbee.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
            
            % Evaluation
            %newbee.Cost=CostFunction(newbee.Position);
            [PRE MRE] =RollingGM(newbee.Position,input,output);
            newbee.Cost= MRE;
            
            % Comparision
            if newbee.Cost<=pop(i).Cost
                pop(i)=newbee;
            else
                C(i)=C(i)+1;
            end
            
        end
        
        % Scout Bees
        for i=1:nPop
            if C(i)>=L
                pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
                [PRE MRE] =RollingGM(pop(i).Position,input,output);
                pop(i).Cost= MRE;
            end
        end
        
        % Update Best Solution Ever Found
        for i=1:nPop
            if pop(i).Cost<=BestSol.Cost
                BestSol=pop(i);
                BestSol.Cost=pop(i).Cost;
            end
        end
        
        BestCost(it)=BestSol.Cost;
        BestAns = BestSol.Position;
        %the best MRE
        Best_score = BestCost(it) ;
    end
    result_rgm=[result_rgm;BestAns BestCost(it)];
end

result_rgm=result_rgm(2:end,1);
for i=1:M
    input=raw_data(i,1:5);
    output=raw_data(i,6);
    lambda=result_rgm(i);
    [PRE MRE lambda]=RollingGM(lambda,input,output);
    PRE_all=[PRE_all,PRE'];
end

forecasting_results=PRE_all(6,2:end);
fprintf('forecasting values for tourism demands from 2013 to 2018\n')
forecasting_results
