function [final_combis] = get_constrained_combis(inputlist, k, n_constrain, printXLS)
% Input #1: path to spreadsheet list, cell array or table
% Input #2: number k of elements that are drawn
% Input #3: constrain on occurrence of each element in final combinations
% Input #4: print back to excel

% check input type
if isa(inputlist, 'string')
    % read excel table
    elements=readtable(inputlist);
    elements = table2cell(elements);
elseif isa(inputlist,'table')
    elements = table2cell(inputlist);
    printIt = false;
elseif isa(inputlist,'cell')
    elements = inputlist;
    printIt = false;
end

% Make sure, the list elements are aligned row-wise
elements = reshape(elements,[numel(elements) 1]);

if k == 0 || n_constrain == 0 || k>size(elements,1)
    error('Please enter valid arguments');
end

ids = 1:size(elements,1); % because the entries in your excel table are indexed by 1,...,n 
combis= nchoosek(ids,k); % all possible combinations of k(=2) elements drawn from  available IDs with no regard to order
% since the result of nchoosek is ordered, shuffle it
shuffle = randperm(size(combis,1)); 
combis = combis(shuffle, :); 
final_combis=[]; % for the final combinations
cnt=zeros(numel(ids),1);

for iID = 1:numel(ids) % loop through IDs
    tmpIdx = find(combis(:,:)==iID); % find (linar) indices of the current ID in the combis array
    [r,~]= ind2sub(size(combis),tmpIdx); % convert linear indices into row and col indices
    for iR = 1:numel(r) % loop through rows of combis array in which current ID occurrs
        iRth_row = r(iR);
        takeIt = true; % true, if the iR-th row of combis array should be added to the final array
        for ik = 1:k % loop through k columns of the combis array in the rows indexed by r
            cnt(combis(iRth_row,ik)) = cnt(combis(iRth_row,ik))+1; % counter of the index in iR-th row and ik-th column of combis array
            if takeIt == true && cnt(combis(iRth_row,ik)) <= n_constrain, takeIt = true; else, takeIt = false; end % if one of the indices contained in the iR-th row occurrs already 4 times in the final array, dont take this row
        end
        if takeIt % if none of the indices in iR-th row of combis array occurrs already 4 times in the list of final combinatinos, add it
            final_combis = [final_combis; combis(iRth_row,:)]; %...add the iR-th row
        else % subtract again from counter of the indices contained in the iRth row, if row has not been taken
            for ik = 1:k % loop through 2 (k) columns of the combis array in the rows indexed by r
                cnt(combis(iRth_row,ik)) = cnt(combis(iRth_row,ik))-1; % counter of the index in iR-th row and ik-th column of combis array
            end
        end
        combis(iRth_row,:) = NaN; % ...and delete it
    end
end

% Shuffle the first and second positions
for idx = 1:size(final_combis,1)
    final_combis(idx,:) = final_combis(idx,randperm(k));
end

% Create cell array with image pairs and export to excel
final_list = cell(size(final_combis)); % create empty cell of same size as final combis array
for ik = 1:k % loop through columns/number of paired images (k=2)
    final_list(:,ik) = elements(final_combis(:,ik)); % select image names from list using ids in the final combis array and fill them into the cell array (column-wise, i.e. loop through columns)
end

if printIt
    pathparts = regexp(inputlist, '\', 'split');
    pathname = strcat(fullfile(pathparts{1:end-1}),'\','final_list.xls');
    % Write into to an excel sheet
    writecell(final_list, pathname);
end
end