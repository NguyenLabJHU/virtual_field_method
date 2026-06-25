function model = findSharedSurfaceNode2Node(model, mat1, mat2)
% Find unique mat1-mat2 shared faces and add 1 psurfid and update surfacesp{1,3} with all [Nsurf,7]
% Mat1 = Outer surface; Mat2 = Inner surface

elements = model.elements;    % Ne x 9: [elnum n1 n2 n3 n4 n5 n6 n7 n8]
elemmat = model.elemmat(:);   % Ne x 1

num_elems = size(elements, 1);
hex_faces = [
    1 2 3 4;
    5 6 7 8;
    1 2 6 5;
    2 3 7 6;
    3 4 8 7;
    4 1 5 8
];
nodes = elements(:,2:9);      % skip first col (element number)
elnum  = elements(:,1);       % element numbers

% Step 1: Build global face map (key: sorted nodes, value: element indices)
face_map = containers.Map('KeyType','char','ValueType','any');
for elem_idx = 1:num_elems
    for f = 1:6
        face_nodes = nodes(elem_idx, hex_faces(f,:)); % <--- NO SORT
        face_nodes_sorted = sort(face_nodes);          % <--- sort only for the key!
        key = sprintf('%d_', face_nodes_sorted);
        if isKey(face_map, key)
            face_map(key) = [face_map(key) elem_idx];
        else
            face_map(key) = elem_idx;
        end
    end
end

% Step 2: Gather all unique mat1-mat2 surface pairs
surf_rows = [];
mat1_idx = find(elemmat == mat1);

for j = 1:length(mat1_idx)
    eidx = mat1_idx(j);
    for f = 1:6
        face_nodes = nodes(eidx, hex_faces(f,:));      % <--- NO SORT
        face_nodes_sorted = sort(face_nodes);          % <--- sort only for the key!
        key = sprintf('%d_', face_nodes_sorted);
        elist = face_map(key);

        % Check for mat2-neighbors only
        for k = 1:length(elist)
            nbr_idx = elist(k);
            if nbr_idx ~= eidx && elemmat(nbr_idx) == mat2
                % Only lowest element reports the interface (unique)
                if elnum(eidx) < elnum(nbr_idx)
                    % surf_rows: [count, n1 n2 n3 n4, elemB_number, which_face]
                    surf_rows(end+1,:) = [0, face_nodes, elnum(eidx), elnum(nbr_idx), f]; %mat1, mat2
                end
            end
        end
    end
end

% Set 'count' column
if ~isempty(surf_rows)
    surf_rows(:,1) = (1:size(surf_rows,1))';
end

% Append a string to psurfid
if iscell(model.psurfid)
    npsurfs = numel(model.psurfid);
else
    model.psurfid = cellstr(model.psurfid);
    npsurfs = numel(model.psurfid);
end
model.psurfid{npsurfs+1} = ['ContactNode2Node - ','Mat',num2str(mat1),'&',num2str(mat2)];

% Update surfacesp
npsurfacesp = numel(model.surfacesp);
model.surfacesp{1,npsurfacesp+1} = surf_rows;

end