function [tpm,gene,variableName,rpkm,rc_table] = rcNormalize(rc_wd, ginfo_wd, varargin)
%
% Transfrom rawcount matrix to TPM and FPKM
% 
% Usage   [tpm,gene,variableName,rpkm,rawcount] = rcNormalize(rc_wd, ginfo_wd)
%         [tpm,gene,variableName,rpkm,rawcount] = rcNormalize(rc_wd, ginfo_wd, 'Remove_ncRNA',path)

% The input rc_wd should be the absolute path of a csv file contained rawcount of
% genes.
% The input ginfo_wd should be the absolute path of a gff3 file contained
% genes' information of certain bacterial strain.
% 
% The optianl input Remove_ncRNA should be a absolute path of a list contained genes
% would be removed. The default of Remove_ncRNA is [] (an empty array);
% 
% The output tpm is a n*m matrix contained normalized reads. The row 
% represents genes. The column represents different samples.
% 
% The output gene is a n*1 cell contained the name of gene which is in the
% same order of row of tpm matrix.
% 
% The output variableName is a m*1 cell contained the name of sample which is in the
% same order of column of tpm matrix.
%
% The output rpkm is a n*m matrix contained normalized reads. The row 
% represents genes. The column represents different samples.
%
% The output rc_table is a n*m table contained rawcount reads from table of rc_wd. The row 
% represents genes. The column represents different samples.
% 
% Written by Ping Shen
% Version 0.0.1. Created on May 5, 2021. Last modified on May 5, 2021.
%

argin = inputParser;
addParamValue(argin,'Remove_ncRNA',[])
parse(argin,varargin{:})
Remove_ncRNA = argin.Results.Remove_ncRNA;

%% load the gene info 
gff3 = GFFAnnotation(ginfo_wd);
gname = {};
idx = strfind(gff3.Attributes, 'ID=gene-');
gidx = [];
for n = 1:length(idx)
    if idx{n}
        gstart = regexp(gff3.Attributes{n},'gene=')+5;
        gend = regexp(gff3.Attributes{n},'gene_biotype')-2;
        gname = [gname;gff3.Attributes{n}(gstart: gend)];
        gidx = [gidx; n];
    end
end
Length = abs(gff3.Start(gidx) - gff3.Stop(gidx))+1;
%% whether remove ncRNA(ssrA & rRNA or undesired genes)
switch Remove_ncRNA
    case []
        % do nothing
    otherwise 
        rg = readtable(Remove_ncRNA,'ReadVariableNames',false);
        removeGene = table2cell(rg);
        for m = 1:numel(removeGene)
            order(m) = find(strcmp(gname, removeGene{m}));
        end
        gname(order,:) = [];
        Length(order,:) = [];
end
%% nomalize the rawcounts to tpm
rc_table = readtable(rc_wd, 'ReadRowNames',true, 'ReadVariableNames',true);
rc_gene = rc_table.Properties.RowNames;
variableName = rc_table.Properties.VariableNames;
rc = table2array(rc_table);

[gene, og, ng] = intersect(rc_gene, gname, 'stable');
L = double(Length(ng));
rawcount = rc(og,:);
sf = sum(rawcount);
rpkm = (rawcount./(L./1000))./(sf./1000000);
tpm = rpkm./(sum(rpkm)./1e6);
