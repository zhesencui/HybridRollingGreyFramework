clear,clc
load raw_data
[M,N] = size(raw_data);
PRE_all=zeros(6,1);
for sample=1:M
    input=raw_data(sample,1:5);
    output=raw_data(sample,6);
    [PRE MRE lambda]=RollingGM(0.5,input,output);
    PRE_all=[PRE_all,PRE'];
end
  format shortg
  forecasting_results=PRE_all(6,2:end);
  fprintf('forecasting values for tourism demands from 2013 to 2018\n')
  forecasting_results
 
  