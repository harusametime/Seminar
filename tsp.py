import numpy as np
import scipy.spatial.distance
from pip._vendor.pyparsing import delimitedList
from ortools.linear_solver import pywraplp
from operator import pos
import itertools
from dask.array.creation import arange
from blaze.server.tests.test_server import cities
import copy
import math


'''
TSP (Traveling Salesman Problem)
'''

class problem:

    def __init__(self, n_cities):

        self.n_cities = n_cities

        self.solver = pywraplp.Solver('SolveIntegerProblem',
                           pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)

        # x, y position of n cities
        pos = np.random.rand(n_cities, 2)

        # distance matrix
        distance = scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(pos, "euclidean"))


        '''
        decision variable
        x_ij indicates whether the path from city i to j is used or not.
        '''
        self.x = []
        for i in range(n_cities):
            for j in range(n_cities):
                x_name = 'x' + str(i) +',' +str(j)
                self.x.append(self.solver.IntVar(0, 1, x_name))

        '''
        objective function to minimize cost
        (objective function to maximize profit is treated as constraint)
        '''
        objective = self.solver.Objective()
        for i in range(n_cities):
            for j in range(n_cities):
                objective.SetCoefficient(self.x[i * n_cities + j], distance[i,j])

        objective.SetMinimization()


        '''
        Other constraints
        const1) only one path comes to each city.
        const2) only one path departs from each city.
        const3) subtour elimination
        '''

        #only one path comes to each city.
        const1 =[]
        for i in range(n_cities):
            const = self.solver.Constraint(1, 1)
            for j in range(n_cities):
                const.SetCoefficient(self.x[i*n_cities + j], 1)
            const1.append(const)

        #only one path departs from each city.
        const2 =[]
        for j in range(n_cities):
            const = self.solver.Constraint(1, 1)
            for i in range(n_cities):
                const.SetCoefficient(self.x[i*n_cities + j], 1)
            const2.append(const)

        # subtour elimination
        const3= []
        for i in range(n_cities):
            # subset is not null and a set of all vertexes
            if i == 0 or i == n_cities-1:
                continue
            else:
                # S is a sub set of all vertexes V
                S1=itertools.combinations(range(n_cities), i)
                for s1 in S1:
                    s2 = copy.deepcopy(s1)
                    const = self.solver.Constraint(-self.solver.infinity(), i-1)
                    for ss1 in list(s1):
                        for ss2 in list(s2):
                            const.SetCoefficient(self.x[ss1*n_cities + ss2], 1)
                    const3.append(const)


    '''
    Solve TSPP problem where one objective function regarding cost is considered
    by epsilon-constraint method,

    verbose = 0 : nothing is shown, 1(default): solution is shown, 2: variable, constraint, solution are shown.
    '''

    def solve(self, verbose = 1):
        result_status = self.solver.Solve()

        # The problem has an optimal solution.
        assert result_status == pywraplp.Solver.OPTIMAL

        # The solution looks legit (when using self.solvers other than
        # GLOP_LINEAR_PROGRAMMING, verifying the solution is highly recommended!).
        assert self.solver.VerifySolution(1e-7, True)

        if verbose == 2:
            print'Number of variables =', self.solver.NumVariables()
            print'Number of constraints =', self.solver.NumConstraints()

            
        if verbose >=1 :
            # The objective value of the solution.
            print 'Optimal objective value = %f' % self.solver.Objective().Value()

            for variable in self.x:
                if variable.solution_value() ==1:
                    print variable.name(),
                    print "=",
                    print variable.solution_value()
            

if __name__ == '__main__':

    n_cities = 15
    p = problem(n_cities)
    p.solve()
