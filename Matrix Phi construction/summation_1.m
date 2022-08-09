clc
clear all
close all
warning off
%%
count = 0;
for i=1:8
    B_k = sparse(195456,484137);% Please update the matrix size based on your selections
    for j=1:16
        count = count+1;
        load(['A_k_',num2str(count),'.mat'])
        B_k = B_k+A_k;
    end
    save(['B_k_',num2str(i)],'B_k')
end