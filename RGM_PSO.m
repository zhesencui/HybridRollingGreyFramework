%_________________________________________________________________________%
%  RGM-MOA framework with particle swarm optimization (RGM-PSO) source codes             %
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
MaxNum = 100;
narvs = 1;
particlesize = 50;
dim=1;
lb=0;
ub=1;
c1 = 2;
c2 = 2;
w =0.6;
vmax = 3*(ub-lb)/10;

[M,N] = size(raw_data);
result_rgm=zeros(1,2);
PRE_all=zeros(6,1);
for sample=1:M
    input=raw_data(sample,1:5);
    output=raw_data(sample,6);
    x = initialization(particlesize,dim,ub,lb);
    v = initialization(particlesize,dim,ub,lb);
    x = x';
    v = v';
    for i= 1:particlesize
        f(i) =inf;
    end
    personalbest_x=x;
    personalbest_faval=f;
    [globalbest_faval,i] = min(personalbest_faval);
    globalbest_x = personalbest_x(i,:);
    k = 1;
    while k <= MaxNum
        for i = 1:particlesize
            [PRE MRE x(i)]=RollingGM(x(i),input,output);
            f(i)=MRE;
            if f(i) < personalbest_faval(i)
                personalbest_faval(i) = f(i);
                personalbest_x(i,:) = x(i,:);
            end
        end
        [globalbest_faval,i] = min(personalbest_faval);
        globalbest_x = personalbest_x(i,:);
        for i = 1:particlesize
            v(i,:) = w*v(i,:) + c1*rand*(personalbest_x(i,:) - x(i,:)) + c2*rand*personalbest_x(i,:);
            for j = 1:narvs
                if v(i,j) > vmax
                    v(i,j) = vmax;
                elseif v(i,j) < -vmax
                    v(i,j) = -vmax;
                end
            end
            x(i,:) = x(i,:) + v(i,:);
        end
        k = k + 1;
        Best_score = globalbest_faval;
    end
    Value1 = 1/globalbest_faval - 1;
    Value1 = num2str(Value1);
    %disp(strcat('the maximum value',' = ', Value1));
    Value2 = globalbest_x;
    Value2 = num2str(Value2);
    %disp(strcat('the maximum x = ', Value2));
    result_rgm=[result_rgm;globalbest_x globalbest_faval];
end

result_rgm=result_rgm(2:end,1);
for i=1:M
    input=raw_data(i,1:5);
    output=raw_data(i,6);
    lambda=result_rgm(i);
    [PRE MRE lambda]=RollingGM(lambda,input,output);
    PRE_all=[PRE_all,PRE'];
end
fprintf('forecasting values for tourism demands from 2013 to 2018\n')
forecasting_results=PRE_all(6,2:end);
forecasting_results
