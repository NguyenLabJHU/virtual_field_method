function model = read_febio_model_test(mydir, mymodel)
% Robust FEBio reader using XML parser
% Works for nested <Material>/<material> tags and MeshDomains
% Returns same fields as your previous function

fullpath = fullfile(mydir, mymodel);

%% --------- Read XML ----------
xDoc = xmlread(fullpath);


%% --------- Nodes ----------
nodeList = xDoc.getElementsByTagName('node');
Nnd = nodeList.getLength;
nodes = zeros(Nnd, 4); % [id x y z]
for k = 0:Nnd-1
    node = nodeList.item(k);
    id = str2double(node.getAttribute('id'));
    txt = char(node.getTextContent);
    nums = sscanf(txt, '%f , %f , %f');
    if numel(nums) < 3
        error('Cannot parse coordinates for node %d', id);
    end
    nodes(k+1, :) = [id, nums(:)'];
end

%% --------- Elements ----------
elBlocks = xDoc.getElementsByTagName('Elements');
Nblocks = elBlocks.getLength;
blockname = strings(1, Nblocks);
Nel_per_block = zeros(1, Nblocks);

totalNel = 0;
% first count total elements
for b = 0:Nblocks-1
    elemBlock = elBlocks.item(b);
    blockname(b+1) = string(elemBlock.getAttribute('name'));
    elems = elemBlock.getElementsByTagName('elem');
    Nel_per_block(b+1) = elems.getLength;
    totalNel = totalNel + Nel_per_block(b+1);
end

nen = 8; % hex8 assumed
elements = zeros(totalNel, nen+1);
elemmat  = zeros(totalNel, 1);
count = 0;
elements_per_block = cell(1, Nblocks);

for b = 0:Nblocks-1
    elemBlock = elBlocks.item(b);
    elems = elemBlock.getElementsByTagName('elem');
    elements_in_block = zeros(Nel_per_block(b+1),nen+1);
    for e = 0:elems.getLength-1
        elem = elems.item(e);
        eid = str2double(elem.getAttribute('id'));
        connStr = char(elem.getTextContent);
        conn = sscanf(connStr, '%d , %d , %d , %d , %d , %d , %d , %d');
        if numel(conn) ~= nen
            error('Element %d in block "%s" does not have %d nodes', eid, char(blockname(b+1)), nen);
        end

        node_ids = conn(:);                                      % [8 x 1] node IDs
        idxs = arrayfun(@(nid) find(nodes(:,1)==nid,1), node_ids);  % indices in node array
        node_xyz = nodes(idxs,2:4);                              % [8 x 3] coordinates
        
        pass = hex8_order_xyz(node_xyz);                    % robust ref order

        count = count + 1;
        elements(count, :) = [eid, conn(:)'];          % store with proper node order

        elements_in_block(e+1,:) = [eid, conn(:)']; 
    end

    elements_per_block{b+1} = elements_in_block;    

end

%% --------- Surfaces ----------
surfNodes = xDoc.getElementsByTagName('Surface');
numSurfs = surfNodes.getLength;
surfaces = struct('name',cell(1,numSurfs),'conn',cell(1,numSurfs));

for s = 0:numSurfs-1
    surf = surfNodes.item(s);
    surfName = string(surf.getAttribute('name')); % or use char() if you want
    quads = surf.getElementsByTagName('quad4');
    quadData = [];  % Collect quad [id, conn] here

    for q = 0:quads.getLength-1
        quad = quads.item(q);
        qid = str2double(quad.getAttribute('id'));
        conn = sscanf(char(quad.getTextContent), '%d , %d , %d , %d');
        quadData = [quadData; qid, conn(:).']; %#ok<AGROW>
    end

    surfaces(s+1).name = surfName;
    surfaces(s+1).conn = quadData;
end

%change from cells to array
for i = 1:numel(surfaces)
     surfaces(i).conn = find_surface_elements(surfaces(i).conn, elements);
end


%% --------- Surface loads (pressure) ----------
surfLoadNodes = xDoc.getElementsByTagName('surface_load');
Nsurfloads = surfLoadNodes.getLength;

% Preallocate
psurfid   = strings(1, Nsurfloads); % surface names
sidesetlc = zeros(1, Nsurfloads);   % load curve IDs
surfacesp = cell(1, Nsurfloads);    % numeric matrices for each surface

for i = 0:Nsurfloads-1
    sl = surfLoadNodes.item(i);

    % Only process pressure loads
    if strcmp(char(sl.getAttribute('type')), 'pressure')
        % Get surface name from surface_load tag
        surfName = char(sl.getAttribute('surface'));
        psurfid(i+1) = string(surfName);

        % Get load curve id from <pressure ...>
        pressureNode = sl.getElementsByTagName('pressure').item(0);
        sidesetlc(i+1) = str2double(pressureNode.getAttribute('lc'));

        % Now find the corresponding <Surface name="...">
        surfaceNodes = xDoc.getElementsByTagName('Surface');
        for j = 0:surfaceNodes.getLength-1
            surfNode = surfaceNodes.item(j);
            if strcmp(char(surfNode.getAttribute('name')), surfName)
                % Collect all child elements (like quad4, tri3, etc.)
                childNodes = surfNode.getChildNodes();
                elems = [];
                for k = 0:childNodes.getLength-1
                    child = childNodes.item(k);
                    if child.getNodeType() == child.ELEMENT_NODE
                        % Split by comma and convert to numbers
                        nums = str2double(strsplit(strtrim(char(child.getTextContent())), ','));
                        elems = [elems; nums]; %#ok<AGROW>
                    end
                end
                nFaces = size(elems,1);
                surfacesp{i+1} = [(1:nFaces)' elems]; % Store numeric matrix
                break;
            end
        end
    end
end

%change from cells to array
for i = 1:numel(surfacesp)
    surfacesp{i} = find_surface_elements(surfacesp{i}, elements);
end

%% --------- Materials ----------
matNodes = xDoc.getElementsByTagName('material');
Nmat = matNodes.getLength;
matid = zeros(1, Nmat);
matname = strings(1, Nmat);
matprop1 = zeros(1, Nmat);  % C1 or equivalent
matprop2 = zeros(1, Nmat);  % C2 or equivalent
matprop3 = zeros(1, Nmat);  % k (bulk/volumetric modulus)
matprop4 = zeros(1, Nmat);  % C3
matprop5 = zeros(1, Nmat);  % C4
matprop6 = zeros(1, Nmat);  % C5
matprop7 = zeros(1, Nmat);  % lam_max

matmodel = strings(1, Nmat);

materials =struct;


for m = 0:Nmat-1
    mat = matNodes.item(m);
    matid(m+1) = str2double(mat.getAttribute('id'));
    materials(m+1).id = matid(m+1);
    matname(m+1) = string(mat.getAttribute('name'));
    materials(m+1).name = matname(m+1);

    % Elastic block (where the material model lives in your XML)
    elastic = mat.getElementsByTagName('elastic');
    if elastic.getLength > 0
        el = elastic.item(0);
    else
        el = mat;
    end

    etype = lower(strtrim(string(el.getAttribute('type')))); % e.g. 'coupled mooney-rivlin'

    % Read parameters present in your file
    c = readScalar(el, {'c','C'});
    C1 = readScalar(el, {'c1','C1'});
    C2 = readScalar(el, {'c2','C2'});
    c3 = readScalar(el, {'c3','C3'});
    c4 = readScalar(el, {'c4','C4'});
    c5 = readScalar(el, {'c5','C5'});
    K  = readScalar(el, {'k','K'});
    k1  = readScalar(el, {'k1','K1'});
    k2  = readScalar(el, {'k2','K2'});
    kappa = readScalar(el, {'Kappa','kappa'});
    gamma = readScalar(el, {'Gamma','gamma'});
    lam_max = readScalar(el, {'lam_max','lambda_max','lamMax'}); % alias support
    E = readScalar(el, {'E'});
    v = readScalar(el, {'v'});


    % Fiber tag: <fiber type="local">1,2</fiber>
    fiber = readFiber(el); % struct with fields: type, indices

    % Fiber tag: <fiber type="local">1,2</fiber>
    mat_axis = readMatAxis(el); % struct with fields: type, indices

    % Decide model
    const_model = "";
    
    etype_lc = lower(strtrim(etype)); % Remove spaces & lowercase

    if strcmpi(etype_lc, 'hgo unconstrained')
        const_model = "hgo-unconstrained";
    elseif strcmpi(etype_lc, 'coupled mooney-rivlin')
        const_model = "coupled-mooney-rivlin";
    elseif strcmpi(etype_lc, 'mooney-rivlin')
        const_model = "mooney-rivlin";    
    elseif strcmpi(etype_lc, 'coupled trans-iso mooney-rivlin')
        const_model = "ti-mooney-rivlin";
    elseif strcmpi(etype_lc, 'neo-hookean')
        const_model = "neo-hookean";
    elseif strcmpi(etype_lc, 'isotropic elastic') || strcmpi(etype_lc, 'isotropic-elastic')
        const_model = "linear-elastic";    
    else
        const_model = "unknown";
    end
    
    matmodel(m+1) = const_model;
    materials(m+1).model = matmodel(m+1);

    % Store material parameters for computeStress and arrays
    params = struct;

    switch const_model
        case {"mooney-rivlin","coupled-mooney-rivlin"}
            params.C1 = defaultZero(C1);
            params.C2 = defaultZero(C2);
            params.K  = defaultZero(K);
            matprop1(m+1) = params.C1;
            matprop2(m+1) = params.C2;
            matprop3(m+1) = params.K;

        case "ti-mooney-rivlin"
            params.C1 = defaultZero(C1);
            params.C2 = defaultZero(C2);
            params.K  = defaultZero(K);
            matprop1(m+1) = params.C1;
            matprop2(m+1) = params.C2;
            matprop3(m+1) = params.K;
            params.c3 = defaultZero(c3);
            params.c4 = defaultZero(c4);
            params.c5 = defaultZero(c5);
            matprop4(m+1) =  params.c3;
            matprop5(m+1) =  params.c4;
            matprop6(m+1) =  params.c5;
            
            if ~isempty(lam_max)
                params.lam_max = lam_max;
                params.flam    = lam_max;
            end

            matprop7(m+1) = params.flam;

            if isfield(fiber, 'type')
                params.fiber_type = fiber.type;
            end
            if isfield(fiber, 'indices')
                params.fiber_indices = fiber.indices;
            end

            % --- Generic child tag extraction ---
            fiberNodes = el.getElementsByTagName('fiber');
            if fiberNodes.getLength > 0
                fiberNode = fiberNodes.item(0);  % Get the first <fiber> under this material's <elastic>
                childNodes = fiberNode.getChildNodes();
                for k = 1:childNodes.getLength
                    child = childNodes.item(k-1);
                    if child.getNodeType == child.ELEMENT_NODE
                        tag = char(child.getNodeName);
                        val = strtrim(char(child.getTextContent));
                        if ~isempty(regexp(val, '^[\d\.\-]+(,[\d\.\-]+)+$', 'once'))
                            val_parsed = str2double(strsplit(val, ','));
                            params.(tag) = val_parsed;
                        else
                            params.(tag) = val;
                        end
                    end
                end
            end

        case "hgo-unconstrained"
            params.c      = defaultZero(c);   % <c>
            params.k1     = defaultZero(k1);   % <k1>
            params.k      = defaultZero(K);    % <k>
            params.k2     = defaultZero(k2);   % <k2>
            params.kappa  = defaultZero(kappa);   % <kappa>
            params.gamma  = defaultZero(gamma);   % <gamma>

            matprop1(m+1) = params.c;
            matprop2(m+1) = params.k1;
            matprop3(m+1) = params.k;
            matprop4(m+1) = params.k2;
            matprop5(m+1) = params.kappa;
            matprop6(m+1) = params.gamma;

            if isfield(mat_axis, 'type')
                params.mat_axis_type = mat_axis.type;
            end
            if isfield(mat_axis, 'indices')
                params.mat_axis_indices = mat_axis.indices;
            end

            % Extract cylindrical axis info if present
            mat_axis_node = el.getElementsByTagName('mat_axis');
            if mat_axis_node.getLength > 0
                axis_node = mat_axis_node.item(0);
                axis_children = axis_node.getChildNodes();
                for k_axis = 1:axis_children.getLength
                    axis_child = axis_children.item(k_axis-1);
                    if axis_child.getNodeType == axis_child.ELEMENT_NODE
                        nodeName = char(axis_child.getNodeName);
                        txt = strtrim(char(axis_child.getTextContent));
                        if ~isempty(regexp(txt,'^[\d\.\-]+(,[\d\.\-]+)+$', 'once'))
                            params.(nodeName) = str2double(strsplit(txt,','));
                        else
                            params.(nodeName) = txt;
                        end
                    end
                end
            end

        case "linear-elastic"
            params.E = defaultZero(E);
            params.v = defaultZero(v);

            matprop1(m+1) = params.E;   % E
            matprop3(m+1) = params.v;   % poisson    

        otherwise % neo-hookean fallback if ever needed
            params.E = defaultZero(E);  
            params.v  = defaultZero(v);
            matprop1(m+1) = params.E;
            matprop3(m+1) = params.v;
    end

    materials(m+1).parameters = params;

end


%% --------- Map element blocks to materials (MeshDomains) ----------
meshDomains = xDoc.getElementsByTagName('MeshDomains');
blockmat = strings(1, Nblocks);
blockmatid = zeros(1, Nblocks);

if meshDomains.getLength > 0
    md = meshDomains.item(0);
    solidDomains = md.getElementsByTagName('SolidDomain');
    for sd = 0:solidDomains.getLength-1
        dom = solidDomains.item(sd);
        dname = string(dom.getAttribute('name'));
        mname = string(dom.getAttribute('mat'));
        % match with blockname
        for b = 1:Nblocks
            if blockname(b) == dname
                blockmat(b) = mname;
                idx = find(matname == mname);
                if ~isempty(idx)
                    blockmatid(b) = idx;
                else
                    blockmatid(b) = 0;
                end
            end
        end
    end
end

% assign element materials
idx = 1;
for b = 1:Nblocks
    nb = Nel_per_block(b);
    if nb>0
        elemmat(idx:idx+nb-1) = blockmatid(b);
        idx = idx + nb;
    end
end

% passing each element to the materials struct
numMaterials = length(materials);

for i = 1:numMaterials
    mat_id = materials(i).id;
    idx = (elemmat == mat_id);          % Find elements for this material
    materials(i).elements = elements(idx, :);  % Each row: [elemID node1 ... node8]
end


%% --------- Readin fibers from the model ----------
fiber_table = [];
[gauss_pts, ~] = get_gauss_points(3, 2); % 3D, order 2 ==> 8 Gauss points
num_gauss = size(gauss_pts, 1);

for i = 1:length(materials)
    if isfield(materials(i).parameters, 'fiber_type')
        fiber_type = materials(i).parameters.fiber_type;

        mat_elements = materials(i).elements;      % [#elems x 9]: elemID, n1,..n8
        elemIDs = mat_elements(:,1);
        elemNodes = mat_elements(:,2:end);         % node indices (#elems x 8)
        fiber_table = [];
        
        if strcmp(fiber_type, 'map')
            % --- MAP: copy the same fiber vector for all 8 gauss points for each element ---
            eldata_nodes = xDoc.getElementsByTagName('ElementData');
            for k = 0:eldata_nodes.getLength-1
                nd = eldata_nodes.item(k);
                if strcmp(char(nd.getAttribute('name')), 'fiber')
                    elem_set_name = char(nd.getAttribute('elem_set'));
                    tmp = read_febio_fibers_elementset(xDoc, elem_set_name);
                    % tmp: [elemID, fx, fy, fz]
                    
                    % Only consider elements in this material
                    [is_mat_elem, loc] = ismember(tmp(:,1), elemIDs);
                    for row = find(is_mat_elem)'
                        elemID = tmp(row,1);
                        fiber_dir = tmp(row,2:4);
                        for gp = 1:num_gauss
                            fiber_table = [fiber_table; elemID, gp, fiber_dir];
                        end
                    end
                end
            end
            materials(i).fibers = fiber_table;  % [elemID, gp, fiberX, fiberY, fiberZ]

        elseif strcmp(fiber_type, 'cylindrical')
            % --- CYLINDRICAL: compute fiber at each Gauss point ---
            m_center = materials(i).parameters.center(:)';
            m_axis = materials(i).parameters.axis(:)';
            m_vector = materials(i).parameters.vector(:)';

            for e = 1:size(elemIDs,1)
                elem_id = elemIDs(e);
                node_ids = elemNodes(e,:);
                % Get xyz for these 8 nodes.
                % nodes: [nodeID, x, y, z]. Let’s build node_xyz as 8x3.
                node_xyz = zeros(8,3);
                for k = 1:8
                    % Find the row in nodes corresponding to node_ids(k)
                    node_xyz(k,:) = nodes(nodes(:,1)==node_ids(k),2:4);
                end

                for gp = 1:num_gauss
                    xi_eta_zeta = gauss_pts(gp,:);
                    N = hex8_shape_functions(xi_eta_zeta); % [8x1]
                    N = N(:)'; % row
                    pos_gp = N * node_xyz; % interpolated xyz

                    % Cylindrical logic:
                    p = pos_gp - m_center;
                    a = m_axis / norm(m_axis);
                    b = p - dot(a,p)*a;
                    b = b / norm(b);

                    e1 = [1 0 0];
                    e1_rot = rodrigues_rot(e1, [0 0 1], a);
                    e1_rot = e1_rot / norm(e1_rot);
                    q = quat_between_vectors(e1_rot, b);
                    r = m_vector / norm(m_vector);
                    fiber_dir = rotate_by_quat(r, q);

                    fiber_table = [fiber_table; elem_id, gp, fiber_dir(:)'];
                end
            end
            materials(i).fibers = fiber_table;
        end
    end
end


%% --------- Readin mat-axis from the model ----------
mat_axis_table = [];
[gauss_pts, ~] = get_gauss_points(3, 2); % 3D, order 2 ==> 8 Gauss points
num_gauss = size(gauss_pts, 1);

for i = 1:length(materials)
    if isfield(materials(i).parameters, 'mat_axis_type')
        mat_axis_type = materials(i).parameters.mat_axis_type;

        mat_elements = materials(i).elements;      % [#elems x 9]: elemID, n1,..n8
        elemIDs = mat_elements(:,1);
        elemNodes = mat_elements(:,2:end);         % node indices (#elems x 8)
        mat_axis_table = [];
        
        if strcmp(mat_axis_type, 'cylindrical')
            % --- CYLINDRICAL: compute mat_axis at each Gauss point ---
            m_center = materials(i).parameters.center(:)';
            m_axis = materials(i).parameters.axis(:)';
            m_vector = materials(i).parameters.vector(:)';

            for e = 1:size(elemIDs,1)
                elem_id = elemIDs(e);
                node_ids = elemNodes(e,:);
                % Get xyz for these 8 nodes.
                % nodes: [nodeID, x, y, z]. Let’s build node_xyz as 8x3.
                node_xyz = zeros(8,3);
                for k = 1:8
                    % Find the row in nodes corresponding to node_ids(k)
                    node_xyz(k,:) = nodes(nodes(:,1)==node_ids(k),2:4);
                end

                for gp = 1:num_gauss
                    xi_eta_zeta = gauss_pts(gp,:);
                    N = hex8_shape_functions(xi_eta_zeta); % [8x1]
                    N = N(:)'; % row
                    pos_gp = N * node_xyz; % interpolated xyz

                    % Cylindrical logic:
                    p = pos_gp - m_center;
                    a = m_axis / norm(m_axis);
                    b = p - dot(a,p)*a;
                    b = b / norm(b);

                    e1 = [1 0 0];
                    e1_rot = rodrigues_rot(e1, [0 0 1], a);
                    e1_rot = e1_rot / norm(e1_rot);
                    q = quat_between_vectors(e1_rot, b);
                    r = m_vector / norm(m_vector);
                    fiber_dir = rotate_by_quat(r, q);

                    Q = fiber_to_axis(fiber_dir);
                    Q_flat = reshape(Q, 1, 9);
                    mat_axis_table = [mat_axis_table; elem_id, gp, Q_flat];

                end
            end
            materials(i).mat_axis = mat_axis_table;
        end
    end
end


%% --------- Tie-Contact (Primary-Secondary surfaces) from .feb XML ----------

% Find all contact definitions
contactList = xDoc.getElementsByTagName('contact');
tie_contact = [];
for i = 0:contactList.getLength-1
    contactNode = contactList.item(i);
    ctype = char(contactNode.getAttribute('type'));
    if strcmpi(ctype, 'tied-elastic')
        % -- Step 1: Get surface_pair name --
        surfPairName = char(contactNode.getAttribute('surface_pair'));

        % -- Step 2: Find <SurfacePair name="surfPairName"> --
        surfacePairList = xDoc.getElementsByTagName('SurfacePair');
        primaryName = '';
        secondaryName = '';
        for sp = 0:surfacePairList.getLength-1
            spNode = surfacePairList.item(sp);
            spName = char(spNode.getAttribute('name'));
            if strcmp(spName, surfPairName)
                primaryName   = strtrim(char(spNode.getElementsByTagName('primary').item(0).getTextContent));
                secondaryName = strtrim(char(spNode.getElementsByTagName('secondary').item(0).getTextContent));
                break;
            end
        end
        if isempty(primaryName) || isempty(secondaryName)
            warning('SurfacePair %s not found or incomplete!', surfPairName)
            continue
        end

        % -- Step 3: Find all <Surface name="..."> definitions --
        surfDefs = xDoc.getElementsByTagName('Surface');
        % --- Primary surface nodes ---
        prim_nodes = [];
        sec_nodes = [];
        sec_face_centers = [];
        for j = 0:surfDefs.getLength-1
            surfDef = surfDefs.item(j);
            surfDefName = strtrim(char(surfDef.getAttribute('name')));
            if strcmp(surfDefName, primaryName)
                quads = surfDef.getElementsByTagName('quad4');
                tris  = surfDef.getElementsByTagName('tri3');
                for q = 0:quads.getLength-1
                    quad = quads.item(q);
                    qid = str2double(quad.getAttribute('id'));
                    conn = sscanf(char(quad.getTextContent), '%d , %d , %d , %d');
                    prim_nodes = [prim_nodes;qid ,conn(:)'];
                end
                for t = 0:tris.getLength-1
                    tri = tris.item(t);
                    qid = str2double(quad.getAttribute('id'));
                    conn = sscanf(char(tri.getTextContent), '%d , %d , %d');
                    prim_nodes = [prim_nodes;qid,conn(:)'];
                end
            end
            if strcmp(surfDefName, secondaryName)
                quads = surfDef.getElementsByTagName('quad4');
                tris  = surfDef.getElementsByTagName('tri3');
                for q = 0:quads.getLength-1
                    quad = quads.item(q);
                    qid = str2double(quad.getAttribute('id'));
                    conn = sscanf(char(quad.getTextContent), '%d , %d , %d , %d');
                    sec_nodes = [sec_nodes; qid, conn(:)'];

                    % --- Compute center of the face for kd-tree ---
                    node_idx = conn(:)';  % node IDs
                    center = mean(nodes(arrayfun(@(nid) find(nodes(:,1)==nid,1), node_idx), 2:4), 1);
                    sec_face_centers = [sec_face_centers; qid, center];

                end
                for t = 0:tris.getLength-1
                    tri = tris.item(t);
                    qid = str2double(quad.getAttribute('id'));
                    conn = sscanf(char(tri.getTextContent), '%d , %d , %d');
                    sec_nodes = [sec_nodes; qid, conn(:)'];

                    % --- Compute center of the face for kd-tree ---
                    node_idx = conn(:)';  % node IDs
                    center = mean(nodes(arrayfun(@(nid) find(nodes(:,1)==nid,1), node_idx), 2:4), 1);
                    sec_face_centers = [sec_face_centers; center];
                end
            end
        end
        tie_contact.primary_name = primaryName;
        tie_contact.secondary_name = secondaryName;
        tie_contact.primary_nodes = find_surface_elements(prim_nodes, elements);
        tie_contact.secondary_elem = find_surface_elements(sec_nodes, elements);
        tie_contact.secondary_face_centers = sec_face_centers;

        break; % only one tied-elastic contact assumed for now
    end
end


%% --------- Output model ----------
model = struct();
model.nodes = nodes;
model.elements = elements;
model.elemmat = elemmat;
model.surfacesp = surfacesp;
model.blockname = blockname;
model.blockmat = blockmat;
model.blockmatid = blockmatid;
model.blockelem = elements_per_block;
model.matname = matname;
model.matmodel = matmodel;
model.matid = matid;
model.matprop = [matprop1;matprop2;matprop3;matprop4;matprop5;matprop6;matprop7]';
model.psurfid = psurfid;
model.tie_contact = tie_contact;
model.fibers = fiber_table;
model.mat_axis = mat_axis_table;
%model = findSharedSurfaceNode2Node(model, 2, 1); %outer surface to inner surface
model.material_model = materials;
model.all_surfaces = surfaces;


end



% --------- Helpers ---------
function val = readScalar(el, tagNames)
    val = [];
    for k = 1:numel(tagNames)
        tag = char(tagNames{k});
        nodes = el.getElementsByTagName(tag);
        if nodes.getLength > 0
            nd = nodes.item(0);
            txt = strtrim(char(nd.getTextContent));
            if ~isempty(txt)
                v = str2double(strrep(txt, ',', ' '));
                if ~isnan(v)
                    val = v;
                    return;
                end
            end
            vAttr = str2double(char(nd.getAttribute('value')));
            if ~isnan(vAttr)
                val = vAttr;
                return;
            end
        end
    end
end

function fiber = readFiber(el)
    fiber = struct;
    nodes = el.getElementsByTagName('fiber');
    if nodes.getLength == 0, return; end
    nd = nodes.item(0);
    typ = strtrim(lower(char(nd.getAttribute('type'))));
    if ~isempty(typ), fiber.type = string(typ); end
    txt = strtrim(char(nd.getTextContent));
    if ~isempty(txt)
        txt = strrep(txt, ',', ' ');
        ids = sscanf(txt, '%d');
        if ~isempty(ids)
            fiber.indices = ids(:).';
        end
    end
end

function mat_axis = readMatAxis(el)
    mat_axis = struct;
    nodes = el.getElementsByTagName('mat_axis');
    if nodes.getLength == 0, return; end
    nd = nodes.item(0);
    typ = strtrim(lower(char(nd.getAttribute('type'))));
    if ~isempty(typ), mat_axis.type = string(typ); end
    txt = strtrim(char(nd.getTextContent));
    if ~isempty(txt)
        txt = strrep(txt, ',', ' ');
        ids = sscanf(txt, '%d');
        if ~isempty(ids)
            mat_axis.indices = ids(:).';
        end
    end
end

function v = defaultZero(x)
    if isempty(x), v = 0; else, v = x; end 
end

function v_rot = rodrigues_rot(v, from, to)
    % Rotates vector v from direction "from" to "to" using Rodrigues' formula
    from = from(:)/norm(from);
    to = to(:)/norm(to);
    axis = cross(from, to);
    sin_ang = norm(axis);
    cos_ang = dot(from, to);
    if sin_ang < 1e-8
        v_rot = v; % no rotation needed
    else
        axis = axis / sin_ang;
        v_rot = v*cos_ang + cross(axis, v)*sin_ang + axis*dot(axis,v)*(1-cos_ang);
    end
end

function q = quat_between_vectors(u, v)
    % Returns a quaternion (as a vector [w x y z]) rotating u to v
    u = u/norm(u); v = v/norm(v);
    half = (u + v); half = half/norm(half);
    w = dot(u, half);
    xyz = cross(u, half);
    q = [w xyz];
end

function v_rot = rotate_by_quat(v, q)
    % v: 1x3 or 3x1 vector
    % q: quaternion as [w x y z]
    % Rotates vector v by quaternion q, returns rotated vector (1x3)
    
    % Turn v into a pure quaternion [0 v]
    vq = [0 v(:)'];
    q_conj = [q(1), -q(2), -q(3), -q(4)];
    vq_rot = quat_mult(quat_mult(q, vq), q_conj);
    v_rot = vq_rot(2:4);
end

function qout = quat_mult(q1,q2)
    % Quaternion multiplication q1*q2, both as [w x y z]
    w1=q1(1); x1=q1(2); y1=q1(3); z1=q1(4);
    w2=q2(1); x2=q2(2); y2=q2(3); z2=q2(4);
    w = w1*w2 - x1*x2 - y1*y2 - z1*z2;
    x = w1*x2 + x1*w2 + y1*z2 - z1*y2;
    y = w1*y2 - x1*z2 + y1*w2 + z1*x2;
    z = w1*z2 + x1*y2 - y1*x2 + z1*w2;
    qout = [w x y z];
end


function Q = fiber_to_axis(fiber)
%FIBER_TO_AXIS  Construct right-handed orthonormal basis Q from a single fiber vector.

fiber = fiber(:)/norm(fiber);

% Default reference axis: global Y
ref = [0;1;0];

% If 'fiber' parallels ref, fall back to global Z
if abs(dot(fiber, ref)) > 0.99
    ref = [0; 0; 1];
end

% Construct v2, the reference axis (orthogonal to fiber)
v2 = ref - dot(ref, fiber)*fiber;
v2 = v2 / norm(v2);

% Third axis, right-handed completion
v3 = cross(fiber, v2);

Q = [fiber v2 v3];
end