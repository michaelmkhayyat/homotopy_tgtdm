function output = callgurobi(interfacedata,map,psdims,dimoutOrig)

% Gurobi 8 does not like seeing options for gurobi 9
if ~isempty(interfacedata.options.gurobi)
    if isequal(interfacedata.options.gurobi,interfacedata.options.default.gurobi)
        interfacedata.options.gurobi = [];
    end
end

options = interfacedata.options;
nOriginal = nnz(interfacedata.variabletype == 0);
model = yalmip2gurobi(interfacedata);

if length(dimoutOrig)==5
    map_solver = map;
    mapc_solver = zeros(size(map_solver));  % Initialize complement map
    
    % Calculate complement map
    for idx = 1:length(map_solver)
        mapc_solver(map_solver(idx)) = idx;
    end
    
    tempVec = zeros(size(interfacedata.hashCache));
    % Note that allu{1} is the joint state vector
    %           allu{2} is the joint control vector
    %           allu{3} are the homotopy variables
    %           allu{4} are the sigma variables
    %           allu{5} are the parameters
    psdims_N = psdims(1);
    psdims_n = psdims(2);
    
    %%%%%%%%%%%%
    
    allu_vec = [zeros(1,max(dimoutOrig{1}))';zeros(1,max(dimoutOrig{2}))';100*ones(1,max(dimoutOrig{3}))';zeros(1,max(dimoutOrig{4}))';zeros(1,max(dimoutOrig{5}))'];
    
    tempVec = allu_vec(mapc_solver);
    tempVec(interfacedata.removethese) = [];
    model.branchingpriority = tempVec;
end
if interfacedata.options.savedebug
    save gurobidebug model
end


if options.showprogress;showprogress('Calling GUROBI',options.showprogress);end
%solvertime = tic;
result = gurobi(model,model.params);
%solvertime = toc(solvertime);
solvertime = result.runtime;%toc(solvertime);

% Gurobi assumes semi-continuous variables only can take negative values so
% we negate semi-continuous violating this
if isfield(result,'x')
    x = result.x;
    if length(x) == length(model.obj)
        if ~isempty(model.NegativeSemiVar)
            x(model.NegativeSemiVar) = -x(model.NegativeSemiVar);
        end
    end
    x = x(1:nOriginal);
else
    x = zeros(nOriginal,1);
end

% On nonconvex models, monomials are included in the list of variables
% simply set those terms to zero, not used anyway
xtemp = zeros(length(interfacedata.c),1);
xtemp(find(interfacedata.variabletype == 0)) = x;
x = xtemp;

problem = 0;
qcDual = [];
if isfield(result,'pi')
    % Gurobi has reversed sign-convention
    D_struc = -result.pi;
    if sum(interfacedata.K.q) > 0
        % YALMIP has appended equalities to model normalized SOCP variables
        % These duals should be moved to the assumed SOCP dual position
        % order: socp, equality, linear element
        % order: equality, linear element, socp
        m = sum(interfacedata.K.q);
        socpDuals = D_struc(1:m);
        D_struc(1:m)=[];
        D_struc = [D_struc(1:interfacedata.K.f+interfacedata.K.l);socpDuals];        
    end
else
    D_struc = [];
end
if isfield(result,'qcpi')
    qcDual = result.qcpi;
end

switch result.status
    case 'OPTIMAL'
        problem =0;
    case 'INFEASIBLE'
        problem = 1;
    case 'UNBOUNDED'
        problem = 2;
    case {'ITERATION_LIMIT','TIME_LIMIT','NODE_LIMIT'}
        problem = 3;
    case {'NUMERIC','SUBOPTIMAL'}
        problem = 4;
    case 'INF_OR_UNBD'
        problem = 12;
    case 'INTERRUPTED';
        problem = 16;
    otherwise
        problem = -1;
end

% Save all data sent to solver?
if interfacedata.options.savesolverinput
	solverinput.model = model;
else
	solverinput = [];
end

% Save all data from the solver?
if interfacedata.options.savesolveroutput    
	solveroutput.result = result; 
else
	solveroutput = [];
end

% Standard interface 
output = createOutputStructure(x,D_struc,[],problem,interfacedata.solver.tag,solverinput,solveroutput,solvertime);
output.qcDual      = qcDual;



