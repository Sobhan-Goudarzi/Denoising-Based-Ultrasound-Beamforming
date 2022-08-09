clc
clear all
close all
warning off
%%
count = 0;
Phi = sparse(195456,484137);% Please update the matrix size based on your selections
for i=1:8
    load(['B_k_',num2str(i),'.mat'])
    Phi = Phi+B_k;
end
save('Phi','Phi')
