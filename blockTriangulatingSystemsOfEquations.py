from collections import deque, defaultdict


def tarjanStronglyConnectedComponents(edge_list):

    """
    Finds all strongly connected components (SCCs) in a directed graph using Tarjan's algorithm.

    Parameters:
        edge_list (List[List[int]]): List of directed edges [u, v].

    Returns:
        List[List[int]]: A list of SCCs, where each SCC is a list of vertex indices.

    Time Complexity:
        O(V + E), where:
            - V = number of unique vertices in the graph
            - E = number of edges
    """

    # Step 1: Collect all unique vertices
    vertices = set()
    for u, v in edge_list:
        vertices.add(u)
        vertices.add(v)
    
    # Step 2: Build adjacency list
    graph = defaultdict(list)
    for u, v in edge_list:
        graph[u].append(v)
    
    # Step 3: Tarjan's algorithm
    index = 0
    indices = {}
    lowlink = {}
    stack = []
    on_stack = set()
    sccs = []

    def strongconnect(v):
        nonlocal index
        indices[v] = lowlink[v] = index
        index += 1
        stack.append(v)
        on_stack.add(v)

        for w in graph[v]:
            if w not in indices:
                strongconnect(w)
                lowlink[v] = min(lowlink[v], lowlink[w])
            elif w in on_stack:
                lowlink[v] = min(lowlink[v], indices[w])

        if lowlink[v] == indices[v]:
            scc = []
            while True:
                w = stack.pop()
                on_stack.remove(w)
                scc.append(w)
                if w == v:
                    break
            sccs.append(scc)

    for v in vertices:
        if v not in indices:
            strongconnect(v)

    return sccs


def hopcraftKartMaximumMatching(U, V, edges):
    """
    Computes the maximum matching in a bipartite graph using the Hopcroft-Karp algorithm.

    The graph is assumed to be bipartite with partitions U and V.
    Each edge is directed from a node in U to a node in V.

    Parameters:
        U (set): Set of vertices in the left partition.
        V (set): Set of vertices in the right partition.
        edges (List[List[int]]): List of edges [u, v], where u ∈ U and v ∈ V.

    Returns:
        List[List[int]]: A list of matched pairs [[u1, v1], [u2, v2], ...].

    Time Complexity:
        O(sqrt(N) * E), where:
            - N = |U| + |V| (total number of vertices)
            - E = number of edges
    """
    
    # Build adjacency list
    graph = defaultdict(list)
    for u, v in edges:
        graph[u].append(v)

    pair_u = {u: None for u in U}
    pair_v = {v: None for v in V}
    dist = {}

    def bfs():
        queue = deque()
        for u in U:
            if pair_u[u] is None:
                dist[u] = 0
                queue.append(u)
            else:
                dist[u] = float('inf')
        dist[None] = float('inf')

        while queue:
            u = queue.popleft()
            if dist[u] < dist[None]:
                for v in graph[u]:
                    if dist[pair_v[v]] == float('inf'):
                        dist[pair_v[v]] = dist[u] + 1
                        queue.append(pair_v[v])
        return dist[None] != float('inf')

    def dfs(u):
        if u is not None:
            for v in graph[u]:
                if dist[pair_v[v]] == dist[u] + 1:
                    if dfs(pair_v[v]):
                        pair_u[u] = v
                        pair_v[v] = u
                        return True
            dist[u] = float('inf')
            return False
        return True

    matching = 0
    while bfs():
        for u in U:
            if pair_u[u] is None:
                if dfs(u):
                    matching += 1

    return [[u, v] for u, v in pair_u.items() if v is not None]


def blockTriangulateJacobian(n, JacobianEntries, showGraph):

    """
    Performs block triangularization of a Jacobian matrix based on its sparsity pattern.

    The method constructs a bipartite graph linking equations and variables, finds a 
    maximum matching, reverses matched edges, identifies strongly connected components 
    (SCCs), and groups them into blocks representing subsystems of equations and variables 
    that can be solved together.

    This is useful for decomposing large systems of nonlinear equations into smaller 
    coupled blocks, improving the efficiency and stability of numerical solvers.

    Parameters:
        n (int): Number of variables (and equations) in the system.
        JacobianEntries (List[Tuple[int, int, float]]): List of nonzero entries in the 
            Jacobian, each as a tuple (equation_index, variable_index, value).

    Returns:
        List[List[List[int]]]: List of blocks, where each block is a pair:
            [equations_list, variables_list], each sorted in ascending order.

    Method Overview:
        1. Constructs a bipartite graph with variable nodes [0..n-1] and equation nodes [n..2n-1].
        2. Computes a maximum matching using Hopcroft-Karp algorithm on the bipartite graph.
        3. Reverses edges in the bipartite graph that belong to the matching.
        4. Finds strongly connected components (SCCs) on the directed graph formed.
        5. Groups variables and equations from each SCC into blocks.

    Time Complexity:
        - Building the bipartite graph: O(m), where m is the number of nonzero Jacobian entries.
        - Maximum matching (Hopcroft-Karp): O(sqrt(n) * m).
        - Reversing edges: O(m).
        - Tarjan's SCC algorithm: O(n + m).
        Overall complexity is dominated by Hopcroft-Karp: O(sqrt(n) * m).

    Notes:
        - Assumes the system has n equations and n variables.
        - JacobianEntries must accurately reflect the sparsity pattern.
        - The function depends on externally defined `hopcraftKartMaximumMatching` and 
          `tarjanStronglyConnectedComponents` functions.
    """

    ## create a bipartite graph which links the equaitons to the variables
    ## variables will be denoted with their index
    ## equations will be denoted with their index plus n
    ## this is done in order to differentiate line
    BipartiteGraph = [[variable, equation+n] for (equation, variable, value) in JacobianEntries]


    ## Create a maximum matching of the bipartite graph
    U = list(range(n))
    V = list(range(n, 2*n))
    M = hopcraftKartMaximumMatching(U, V, BipartiteGraph)


    ## Flip the edges in the Bipartite graph if the edge is in the maximum matching
    for edge in M:
        index = BipartiteGraph.index(edge)
        BipartiteGraph[index] = [edge[1], edge[0]]


    ## find the strongly connected components of the masked bipartite graph
    stronglyConnectedComponents = tarjanStronglyConnectedComponents(BipartiteGraph)


    ## format the strongly connected components in to blocks
    blocks, currentVariables, currentEquations = [], [], []
    for stronglyConnectedComponent in reversed(stronglyConnectedComponents):

        for elem in stronglyConnectedComponent:
            if elem < n:
                ## this is a variable
                currentVariables.append(elem)
            else:
                ## this is an equation
                currentEquations.append(elem - n)
        
        ## check if a block has been made
        if (len(stronglyConnectedComponent) != 1) or ((len(currentEquations) == len(currentVariables)) and (len(currentEquations) != 0)):
            ## either:
            ## There are more than one element
            ## therefore this must be a proper strongly connected components
            ## this represents n equations with n unknwons
            ## bank the current equations and variables and create new lists
            ## or:
            ## this is not a strongly connected component
            ## but the length of the variables and the len of the equations are equal
            ## therefore these n equaitons and variable can be solved together

            ## sort the lists of equations and variables
            currentEquations.sort()
            currentVariables.sort()
            
            ## store the lists in the list of blocks
            blocks.append([currentEquations, currentVariables])
            
            ## reset the lists
            currentVariables = []
            currentEquations = []
        
    if showGraph:
        import networkx as nx
        from pyvis.network import Network

        ## create the graph
        G = nx.DiGraph()
        for i in range(2*n):
            shape = 'circle' if i < n else 'square'
            if i < n:
                name = f"x{i + 1}"
            else:
                name = f"f{i-n + 1}"
            G.add_node(i, shape = shape, label = str(name))


        for edge in BipartiteGraph:
            color = 'blue'
            if [edge[1], edge[0]] in M: color = 'red'
            G.add_edge(edge[0], edge[1], color = color)

        nx.draw(G, with_labels = True)
        nt = Network(directed=True)
        nt.from_nx(G)
        nt.toggle_physics(False)
        nt.options.edges.smooth.enabled = False
        nt.show('nx.html', notebook=False)

    return blocks




# http://www.mmrc.iss.ac.cn/~xgao/paper/00-order.pdf

n = 8
def JacobianEntries(x):
    ## the jacobian is in general dependent on x
    ## however, in this example it is not

    jacobianEntries = [
        [0, 2, +1],
        [0, 7, +1],
        [1, 6, +1],
        [2, 4, +1],
        [2, 5, -1],
        [2, 6, -1],
        [3, 0, +1],
        [3, 3, +1],
        [3, 5, -1],
        [4, 1, +1],
        [4, 7, +1],
        [5, 2, +1],
        [5, 4, -1],
        [5, 6, +1],
        [6, 3, +1],
        [7, 0, +1],
        [7, 5, +1],
        [7, 6, +1],
    ]
    return jacobianEntries

def Equations(x):
    equations = []
    equations.append((x[2] + x[7]) - 11)
    equations.append((x[6]) - 7)
    equations.append((x[4] - x[5] - x[6]) + 8)
    equations.append((x[0] + x[3] - x[5]) + 1)
    equations.append((x[1] + x[7]) - 10)
    equations.append((x[2] - x[4] + x[6]) - 5)
    equations.append((x[3]) - 4)
    equations.append((x[0] + x[5] + x[6]) - 14)
    
    
    return equations

x = [0] * n

blocks = blockTriangulateJacobian(n,JacobianEntries(x), False)


for equationsIndex, variableIndex in blocks:   

    ## pull out the equations of the block
    equations = Equations(x)
    equations = [equations[elem] for elem in equationsIndex]


    
    ## pull out the jacobian entries of the block
    jacobianEntries = JacobianEntries(x)
    _jacobianEntries = []
    for entry in jacobianEntries:

        equation, variable, value = entry

        if not equation in equationsIndex: continue

        if not variable in variableIndex: continue

        _jacobianEntries.append(entry)
    jacobianEntries = _jacobianEntries

    ## solve the system
    import numpy as np
    n = len(equations)
    J = np.zeros([n,n])
    for entry in jacobianEntries:
        
        equation, variable, value = entry
        
        i = equationsIndex.index(equation)
        j = variableIndex.index(variable)
        J[i,j] = value
    
    ## TODO why is it -J in stead of J?
    var = np.linalg.solve(-J, equations)

    ## insert the solved variables in to x
    for v, index in zip(var, variableIndex):
        x[index] = float(v)
    
print('x', x)
