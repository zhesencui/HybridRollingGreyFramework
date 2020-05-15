%_________________________________________________________________________%
%  RGM-MOA framework with grey wolf optimizer (RGM-GWO) source codes             %
%                                                                         %
%  Developed in MATLAB                                      %
%                                                                         %
%  Author and programmer: Zhesen Cui                              %
%                                                                         %
%         E-Mail: cuizhesen@gmail.com
%_________________________________________________________________________%
close all
clear
clc
format compact

load raw_data
SearchAgents_no=50; % Number of search agents
Max_iteration=100; % Maximum numbef of iterations
dim=1; % number of your variables
lb=0;
ub=1;
% v = 5;

[M,N] = size(raw_data);
result_rgm=zeros(1,1);
PRE_all=zeros(6,1);
for sample=1:M
    input=raw_data(sample,1:5);
    output=raw_data(sample,6);
    % initialize alpha, beta, and delta_pos
    Alpha_pos=zeros(1,dim);
    Alpha_score=inf; % change this to -inf for maximization problems
    
    Beta_pos=zeros(1,dim);
    Beta_score=inf; % change this to -inf for maximization problems
    
    Delta_pos=zeros(1,dim); %
    Delta_score=inf; % change this to -inf for maximization problems
    
    %Initialize the positions of search agents
    Positions=initialization(SearchAgents_no,dim,ub,lb);
    
    %在这里要用的是一个20*1的矩阵，所以需要再转置一下
    Positions = Positions';
    
    Convergence_curve=zeros(1,Max_iteration);
    
    l=0; % Loop counter
    
    % Main loop
    
    while l<100
        for i=1:size(Positions,1)
            % Return back the search agents that go beyond the boundaries of the search space
            
            Flag4ub=Positions(i,:)>ub;
            Flag4lb=Positions(i,:)<lb;
            
            Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb; % ~表示取反
            
            [PRE MRE Positions(i,:)]=RollingGM(Positions(i,:),input,output);
            
            fitness=MRE;
            
            % Update Alpha, Beta, and Delta
            if fitness<Alpha_score
                Alpha_score=fitness; % Update alpha
                Alpha_pos=Positions(i,:);
            end
            
            if fitness>Alpha_score && fitness<Beta_score
                Beta_score=fitness; % Update beta
                Beta_pos=Positions(i,:);
            end
            
            if fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score
                Delta_score=fitness; % Update delta
                Delta_pos=Positions(i,:);
            end
        end
        
        a=2-l*((2)/Max_iteration); % a decreases linearly fron 2 to 0
        
        % Update the Position of search agents including omegas
        for i=1:size(Positions,1)
            for j=1:size(Positions,2)
                
                r1=rand(); % r1 is a random number in [0,1]
                r2=rand(); % r2 is a random number in [0,1]
                
                A1=2*a*r1-a;
                C1=2*r2;
                
                D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j));
                X1=Alpha_pos(j)-A1*D_alpha;
                
                r1=rand();
                r2=rand();
                
                A2=2*a*r1-a;
                C2=2*r2;
                
                D_beta=abs(C2*Beta_pos(j)-Positions(i,j));
                X2=Beta_pos(j)-A2*D_beta;
                
                r1=rand();
                r2=rand();
                
                A3=2*a*r1-a;
                C3=2*r2;
                
                D_delta=abs(C3*Delta_pos(j)-Positions(i,j));
                X3=Delta_pos(j)-A3*D_delta;
                
                Positions(i)=(X1+X2+X3)/3;
                
            end
        end
        l=l+1;
        Convergence_curve(l)=Alpha_score;
        Best_score = Alpha_score;
    end
    
    bestc=Alpha_pos(1,1);
    bestGWOaccuarcy=Alpha_score;
    result_rgm=[result_rgm;bestc];
    
end

result_rgm=result_rgm(2:end,:);
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

