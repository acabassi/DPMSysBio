close all
hypers = importdata('data-30-2_Results_Chain1.csv');
hypers = hypers.data(:,62:end);
%clear
%close all
%load('hypers.mat')
for j = 1:6
    subplot(2,3,j)
    hist(hypers(:,j))
end


