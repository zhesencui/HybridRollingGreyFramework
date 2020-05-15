clc,clear
load compar_raw_data
[M,N] = size(compar_raw_data);
for i=1:(M-5)
    input=compar_raw_data(1:4+i)';
    output=compar_raw_data(5+i);
    lambda=0.5;
    [PRE MRE lambda]=RollingGM(lambda,input,output);
    PRE_all(i)=PRE(end);
end
format shortg
forecasting_results=PRE_all';

fprintf('forecasting values for tourism demands from 2013 to 2018\n')
forecasting_results
