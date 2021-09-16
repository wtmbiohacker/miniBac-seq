function ribo_ratio = plotRiboRatio(counts, gene_name, Variable_name, ribofile)
% This function is to calculate the ratio of ribosomal RNA counts in all 
% varibales and display the ratio with bar plot.
%
% Usage: ribo_ratio = plotRiboRatio(counts, gene_name, Variable_name, ribofile)
%
% Input:
% The input counts is a N*M matrix, storing the read counts of genes. The N
% represents the number of genes. The M represents the number of variables.
%
% The input gene_name is name of genes. The order is same as rows of counts.
%
% The input Variable_name is name of variables.The order is same as columns
% of counts.
%
% The input ribofile is a text file contained genes coded ribosomal RNA.
%
% Output:
% The output gene_ratio is a 1*M vector contained the the ratio of 
% ribosomal RNA counts in all variables. The M represents the number of 
% variables.
% Written by Ping Shen @ Tsinghua University
% Version 0.0.1. Created on Dec 2, 2020. Last modified on Dec 2, 2020.

% load ribosomal RNA
ribo = readtable(ribofile,'ReadVariableNames',false);
ribo = table2cell(ribo);
order = zeros(1,numel(ribo));
for m = 1:numel(ribo)
    order(m) = find(strcmp(gene_name, ribo{m}));
end

ribo_ratio = zeros(1,numel(Variable_name));
for i = 1:numel(Variable_name)
    ribo_ratio(i) = sum(counts(order,i)) / sum(counts(:,i)) * 100;
end

figure
set(gcf,'position', [100 100 1100 650])
b = bar(ribo_ratio,0.6);
b.EdgeColor = 'white';
b.FaceColor = [0.5 0.5 0.5];
set(gca,'Xtick',[1:numel(Variable_name)],'Xticklabel',Variable_name,'FontSize',18,'YGrid','on')
xtickangle(330)
xlabel('Library','FontSize',20)
ylabel('Percentage of ribosomal reads (%)', 'FontSize',25)
ylim([0 max(ribo_ratio)+1])
box on