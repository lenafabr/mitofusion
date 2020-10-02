%filename = '../test/results_initrun/initrun.*.snap.out';


function [positions, proteins,timevals,nregions,nclust,ids,states,idtags] = readSnap(filename)
%function snapinfo = readSnap(filename)
%%

%filename = sprintf('../test/results_initrun2/initrun.%d.snap.out',fileno);
%data = dlmread(filename);

of = fopen(filename);
data = zeros(0,6);
ct=0;
while 1
    line = fgetl(of);% read in a line
    if (line==-1)
        break
    end
    vals = strsplit(line);
    goodvals = find(cellfun(@(x) ~isempty(x)>0, vals));    
    vals = vals(goodvals);
    ct = ct+1;
    isnum = cellfun(@(x) ~isnan(str2double(x)), vals);
    n = nnz(isnum);
    data(ct,1:n) = cellfun(@(x) str2num(x), vals(isnum));
end
fclose(of)

%%
curind = 1;
sc = 1; % snapshot line counter

while (curind < size(data,1)) % cycle over snapshot
    %[sc curind]
    dataline = data(curind,:);
    timevals(sc) = dataline(1);
    nregions(sc) = dataline(2); % number of regions, info(2)
    nclust(sc) = dataline(3);

    % current cluster index
    ids{sc} = data(curind+1:curind+nclust(sc),1);
    positions{sc} = data(curind+1:curind+nclust(sc),2);
    % protein{sc,fileno} has the information on protein content of all
    % clusters at sc timestep of fileno fileproteins(
    proteins{sc} = data(curind+1:curind+nclust(sc),3);
    states{sc} = data(curind+1:curind+nclust(sc),4);
    % tag ids, unique to each individual cluster
    idtags{sc} =  data(curind+1:curind+nclust(sc),5);
    
    curind = curind+ nclust(sc) +1;
    sc=sc+1;
    
end

end
