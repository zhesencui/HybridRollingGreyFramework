%_________________________________________________________________________%
%  RGM-MOA framework with Whale Optimization Algorithm (RGM-WOA) source codes             %
%                                                                         %
%  Developed in MATLAB                                      %
%                                                                         %
%  Author and programmer: Zhesen Cui                              %
%                                                                         %
%         E-Mail: cuizhesen@gmail.com
%_________________________________________________________________________%

function main()
clc;clear all; close all;

load raw_data
E0 = 0.001;
Max_iter = 100;
narvs = 1;
SearchAgents_no = 50;
dim=1;
lb=0;
ub=1;
[M,N] = size(raw_data);
result_rgm=zeros(1,2);
PRE_all=zeros(6,1);
for sample=1: M
    input=raw_data(sample,1:5);
    output=raw_data(sample,6);
    % initialize position vector and score for the leader
    Leader_pos=zeros(1,dim);
    Leader_score=inf; %change this to -inf for maximization problems
    
    %Initialize the positions of search agents
    Positions=initialization(SearchAgents_no,dim,ub,lb);
    Positions=Positions';
    
    Convergence_curve=zeros(1,Max_iter);
    
    t=0;% Loop counter
    
    % Main loop
    while t<Max_iter
        for i=1:size(Positions,1)
            % Return back the search agents that go beyond the boundaries of the search space
            Flag4ub=Positions(i,:)>ub;
            Flag4lb=Positions(i,:)<lb;
            Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
            
            [PRE MRE Positions(i,:)]=RollingGM(Positions(i,:),input,output);
            
            % Calculate objective function for each search agent
            %fitness=fobj(Positions(i,:));
            fitness=MRE;
            % Update the leader
            if fitness<Leader_score % Change this to > for maximization problem
                Leader_score=fitness; % Update alpha
                Leader_pos=Positions(i,:);
            end
        end
        
        a=2-t*((2)/Max_iter); % a decreases linearly fron 2 to 0 in Eq. (2.3)
        
        % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
        a2=-1+t*((-1)/Max_iter);
        
        % Update the Position of search agents
        for i=1:size(Positions,1)
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]
            
            A=2*a*r1-a;  % Eq. (2.3) in the paper
            C=2*r2;      % Eq. (2.4) in the paper
            
            b=1;               %  parameters in Eq. (2.5)
            l=(a2-1)*rand+1;   %  parameters in Eq. (2.5)
            
            p = rand();        % p in Eq. (2.6)
            
            for j=1:size(Positions,2)
                
                if p<0.5
                    if abs(A)>=1
                        rand_leader_index = floor(SearchAgents_no*rand()+1);
                        X_rand = Positions(rand_leader_index, :);
                        D_X_rand=abs(C*X_rand(j)-Positions(i,j)); % Eq. (2.7)
                        Positions(i,j)=X_rand(j)-A*D_X_rand;      % Eq. (2.8)
                        
                    elseif abs(A)<1
                        D_Leader=abs(C*Leader_pos(j)-Positions(i,j)); % Eq. (2.1)
                        Positions(i,j)=Leader_pos(j)-A*D_Leader;      % Eq. (2.2)
                    end
                    
                elseif p>=0.5
                    
                    distance2Leader=abs(Leader_pos(j)-Positions(i,j));
                    % Eq. (2.5)
                    Positions(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos(j);
                    
                end
                
            end
        end
        t=t+1;
        Convergence_curve(t)=Leader_score;
        [globalbest_faval,i] = min(Convergence_curve);
        %;Leader_score represent minimum MRE
        Best_score = Leader_score;
        Best_score = globalbest_faval;
        Best_pos = Leader_pos;
    end
    result_rgm=[result_rgm;Best_pos Best_score];
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

