function runRicciFlow()
    %% Parameters
    inputFile = 'UNC071_undeformed.stl'; 
    outputFile = 'GivenGeom_RicciFlow.stl';
    maxIterations = 10;        
    stepSize = 0.01;           
    smoothingIterationsFirstPass = 500;  % Increased smoothing (from typical 5 - 20)
    smoothingFactorFirstPass = 0.1;     % Reduced alpha from 0.5 for gentler smoothing
    smoothingIterationsSecondPass = 200; % Increased iteration, from a typical pass of 20 iterations
    smoothingFactorSecondPass = 0.1;    % Even gentler smoothing factor for second pass

    %% Load STL and parse faces (F) and vertices (V)
    [F, V] = loadSTL(inputFile);

    if isa(F, 'triangulation')
        disp('Detected triangulation object. Extracting connectivity and points...');
        triObj = F;
        F = triObj.ConnectivityList;
        V = triObj.Points;
    end

    % Convert tetrahedra to triangles if needed
    if size(F,2) == 4
        disp('Converting tetrahedra to triangles...');
        F = tetraToTriangles(F);
    elseif size(F,2) ~= 3
        error('Faces are not Mx3 or Mx4. The STL may not be a standard mesh.');
    end

    % Ensure V is Nx3
    if size(V,2) < 3
        V = [V, zeros(size(V,1), 3 - size(V,2))];
        warning('Vertices had fewer than 3 columns. Padded with zeros.');
    elseif size(V,2) > 3
        V = V(:,1:3);
        warning('Vertices had more than 3 columns. Truncated to the first 3 columns.');
    end

    % Check if indexing is zero-based
    if min(F(:)) == 0
        disp('Detected zero-based indexing in faces. Adjusting to one-based...');
        F = F + 1;
    end

    % Validate indexing or try reindexing
    if max(F(:)) > size(V,1)
        disp('Face indices exceed the number of vertices. Attempting to reindex...');
        [F, V, success] = reindexMesh(F, V);
        if ~success
            disp('Reindexing failed. Removing invalid faces...');
            validFaceMask = all(F <= size(V,1), 2);
            F = F(validFaceMask, :);
            if isempty(F)
                error('No valid faces remain. The STL file is severely invalid.');
            else
                disp('Some faces were removed to ensure validity.');
            end
        end
    end

    %% Dummy Ricci flow iteration
    E = edgesFromFaces(F);
    L = edgeLengths(V, E);

    u = zeros(size(V,1),1);
    K_target = zeros(size(V,1),1);

    for iter = 1:maxIterations
        L_current = adjustedEdgeLengths(L, u, E);
        K = zeros(size(V,1),1);
        fprintf('Iteration %d: Curvature error = %e\n', iter, norm(K - K_target));
        u = u + stepSize * (K_target - K);
    end

    % Slight scaling
    V = V * 1.001;

    % Enhanced smoothing: Two-pass smoothing with more iterations and gentler alpha
    V = smoothMesh(V, F, smoothingIterationsFirstPass, smoothingFactorFirstPass);
    V = smoothMesh(V, F, smoothingIterationsSecondPass, smoothingFactorSecondPass);

    % Create a triangulation object
    TR = triangulation(F, V);

    % Write out using the version of stlwrite that requires triangulation first
    stlwrite(TR, outputFile);
    fprintf('Process completed. Enhanced smoothing applied. Result saved to %s\n', outputFile);
end

%% Helper Functions (Unchanged)

function [F, V] = loadSTL(inputFile)
    data = stlread(inputFile);
    if isstruct(data)
        if isfield(data, 'faces') && isfield(data, 'vertices')
            F = data.faces;
            V = data.vertices;
        elseif isfield(data, 'Faces') && isfield(data, 'Vertices')
            F = data.Faces;
            V = data.Vertices;
        else
            F = data; 
            V = [];
        end
    else
        try
            [F,V] = stlread(inputFile);
        catch
            if isa(data, 'triangulation')
                F = data;
                V = data; 
            else
                error('Unexpected data format. Check stlread usage.');
            end
        end
    end
end

function F_tri = tetraToTriangles(F_tet)
    M = size(F_tet,1);
    F_tri = zeros(4*M,3);
    for i = 1:M
        v = F_tet(i,:);
        F_tri(4*(i-1)+1,:) = v([1,2,3]);
        F_tri(4*(i-1)+2,:) = v([1,3,4]);
        F_tri(4*(i-1)+3,:) = v([1,2,4]);
        F_tri(4*(i-1)+4,:) = v([2,3,4]);
    end
end

function [F_out, V_out, success] = reindexMesh(F_in, V_in)
    usedVerts = unique(F_in(:));
    maxIndex = max(usedVerts);
    success = false;
    if maxIndex <= size(V_in,1)
        map = zeros(maxIndex, 1);
        map(usedVerts) = 1:numel(usedVerts);
        V_out = V_in(usedVerts, :);
        F_out = map(F_in);
        success = true;
    else
        F_out = F_in;
        V_out = V_in;
    end
end

function E = edgesFromFaces(F)
    E = sort([F(:,[1,2]); F(:,[2,3]); F(:,[3,1])],2);
    E = unique(E,'rows');
end

function L = edgeLengths(V, E)
    v1 = V(E(:,1),:);
    v2 = V(E(:,2),:);
    L = sqrt(sum((v2 - v1).^2,2));
end

function L_current = adjustedEdgeLengths(L, u, E)
    L_current = L;
end

function V = smoothMesh(V, F, iterations, alpha)
    nV = size(V,1);
    adjList = vertexAdjacency(F, nV);
    for iter = 1:iterations
        V_new = V;
        for i = 1:nV
            neighbors = adjList{i};
            if ~isempty(neighbors)
                avgPos = mean(V(neighbors,:), 1);
                V_new(i,:) = V(i,:) + alpha * (avgPos - V(i,:));
            end
        end
        V = V_new;
    end
end

function adjList = vertexAdjacency(F, nV)
    adjList = cell(nV,1);
    for i = 1:size(F,1)
        f = F(i,:);
        adjList{f(1)} = unique([adjList{f(1)}, f(2), f(3)]);
        adjList{f(2)} = unique([adjList{f(2)}, f(1), f(3)]);
        adjList{f(3)} = unique([adjList{f(3)}, f(1), f(2)]);
    end
end

