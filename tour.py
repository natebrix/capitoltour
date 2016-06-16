# Copyright 2016, Gurobi Optimization, Inc.
# Modified by Nathan Brixius, June 2016. Thanks Gurobi!

# Solve a traveling salesman problem. The base MIP model only includes
# 'degree-2' constraints, requiring each node to have exactly
# two incident edges.  Solutions to this model may contain subtours -
# tours that don't visit every city. We add constraints to eliminate
# subtours until they're all gone. In Gurobi's version of this code
# they use 'lazy constraints', which is more efficient. This version
# of the code runs on the cloud, where lazy constraints are not
# supported. Even this more inefficient version is able to find
# probably optimal solutions much more quickly than a genetic algorithm
# can find decent suboptimal ones.
#
# I have also added an option to find the best tour involving K of the N
# state capitols.

import sys
import math
#from gurobipy import *
import csv

tsp_path = r'C:\Users\Nathan\OneDrive\Public\TSP\\'
tsp_outpath = tsp_path + "output\\"
tsp_file = tsp_path + r'my-waypoints-dist-dur.tsv'

class CapitolTour:
    def __init__(self, tsp_file=r'C:\Users\Nathan\OneDrive\Public\TSP\my-waypoints-dist-dur.tsv'):
        self.n = 0 # todo
        self.k = 0
        self.find_k = False
        self.tsp_file = tsp_file

    # Adds a constraint to eliminate subtours (if necessary).
    #
    # It is preferable to do this using lazy constraints (see the Gurobi example tsp.py)
    # but I am running on the cloud and Gurobi Cloud doesn't support that.
    def subtour_elimination(self, model):
        n = self.n
        selected = []
        # make a list of edges selected in the solution
        for i in range(n):
            sol = model.getAttr('x', [model._x[i, j] for j in range(n)])
            selected += [(i,j) for j in range(n) if sol[j] > 0.5]

        # find the shortest cycle in the selected edge list
        tour = self.subtour(selected, self.get_subtour_visited(model))
        if len(tour) < self.k:
            # add a subtour elimination constraint
            expr = 0
            for i in range(len(tour)):
                for j in range(i+1, len(tour)):
                    expr += model._x[tour[i], tour[j]]
            # In real life you should do this instead and call m.optimize() using this
            # funciton as a callback (see the included Gurobi example)
            # model.cbLazy(expr <= len(tour) - 1)
            print("***** Adding constraint to eliminate tour of length %d" % len(tour))
            model.addConstr(expr <= len(tour) - 1)
            return True
        print("***** No subtour found")
        return False

    # Initializes which cities should be examined for subtours. We need this because
    # sometimes we are looking for an optimal tour of size K among N cities.
    def get_subtour_visited(self, model):
        visited = [False] * self.n
        if self.find_k:
            y = model.getAttr('x', [model._y[i] for i in range(self.n)])
            visited = [y_i < 0.5 for y_i in y] # disregard inactive cities
        return visited

    # Given a list of edges, finds the shortest subtour
    def subtour(self, edges, visited):
        cycles = []
        lengths = []
        selected = [[] for i in range(self.n)]
        for x, y in edges:
            selected[x].append(y)
        while True:
            current = visited.index(False)
            thiscycle = [current]
            while True:
                visited[current] = True
                neighbors = [x for x in selected[current] if not visited[x]]
                if len(neighbors) == 0:
                    break
                current = neighbors[0]
                thiscycle.append(current)
            cycles.append(thiscycle)
            lengths.append(len(thiscycle))
            if sum(lengths) == self.k:
                break
        return cycles[lengths.index(min(lengths))]

    def get_dist(self, tsp):
        with open(tsp, 'rb') as tspfile:
            r = csv.reader(tspfile, delimiter='\t')
            d = [row for row in r]

        d = d[1:] # skip header row
        locs = set([r[0] for r in d]) | set([r[1] for r in d])
        loc_map = {l:i for i, l in enumerate(locs)}
        idx_map = {i:l for i, l in enumerate(locs)}
        dist = [(loc_map[r[0]], loc_map[r[1]], r[2]) for r in d]
        return dist, idx_map

    def get_tour(self, tsp):
        with open(tsp, 'rb') as tspfile:
            r = csv.reader(tspfile, delimiter='\t')
            d = [row[0] for row in r]
        return d

    def get_latlong(self, latlong):
        with open(latlong, 'rb') as f:
            r = csv.reader(f, delimiter=',')
            next(r, None)
            return [(float(row[1]), float(row[2])) for row in r]


    def dist_from_coords(self, dist, n):
        points = []
        for i in range(n):
            points.append([0] * n)
        for i, j, v in dist:
            points[i][j] = points[j][i] = float(v)
        return points

    def leg_distances(self, tour, dist):
        legs = []
        i = tour[-1]
        for j in tour:
            legs.append(dist[i][j])
            i = j
        return legs

    # Find an optimal tour among the cities specified in tsp_file.
    # k is an optional parameter. If specified, an optimal tour of length k
    # will be returned, otherwise a tour with all cities will be returned.
    def find_tour(self, k=-1):
        coords, idx_map = self.get_dist(self.tsp_file)
        self.n = n = len(idx_map)
        dist = self.dist_from_coords(coords, n)
        self.find_k = k > 0
        self.k = k if k > 0 else self.n

        m = Model()
        m.params.OutputFlag = 0

        # Create variables.
        # x[i, j] is 1 if the edge i->j is on the optimal tour, and 0 otherwise.
        x = {}
        for i in range(n):
            for j in range(i+1):
                x[i,j] = m.addVar(obj=dist[i][j], vtype=GRB.BINARY,
                                     name='x'+str(i)+'_'+str(j))
                x[j,i] = x[i,j]

        # Nathan added this variable. y[i] is 1 if state capitol i is on the
        # optimal tour. We only need this variable if we are looking for
        # an optimal tour that involves K of the N state capitols.
        y = {}
        if self.find_k:
            for i in range(n):
                y[i] = m.addVar(vtype=GRB.BINARY, name='y'+str(i))
        m.update()

        if not self.find_k:
            # Add degree-2 constraint, and forbid loops
            for i in range(n):
                m.addConstr(quicksum(x[i, j] for j in range(n)) == 2)
                x[i,i].ub = 0
        else:
            # choose exactly k capitols.
            m.addConstr(quicksum(y[i] for i in range(n)) == k)
            # if a state capitol is on a tour, the degree-2 constraint applies.
            # a state capitol is on a tour only when y[i]==1.
            for i in range(n):
                m.addConstr(quicksum(x[i, j] for j in range(n)) == 2 * y[i])
                x[i, i].ub = 0

        m.update()

        m._x = x
        m._y = y
        m.optimize()

        for i in range(1000): # <-- bad; don't care.
            m.update()
            m.optimize()
            go = self.subtour_elimination(m)
            if not go:
                break

        solution = m.getAttr('x', x)
        selected = [(i,j) for i in range(n) for j in range(n) if solution[i,j] > 0.5]
        tour = self.subtour(selected, self.get_subtour_visited(m))
        assert len(tour) == self.k

        print('')
        print('Optimal tour: %s' % str(tour))
        print('Optimal cost: %g' % m.objVal)
        print('')
        for loc in [idx_map[i] for i in tour]:
            print(loc)

        # let's check our math.
        leg_dist = self.leg_distances(tour, dist)
        tour_dist = sum(leg_dist)
        print(tour_dist)
        assert abs(tour_dist - m.objVal) <= 1e-4

        return tour, leg_dist, [idx_map[i] for i in tour]

    # Finds 'efficient frontier': minimal length tour involving K cities for K = 3..48.
    def find_efficient_frontier(self):
        r = []
        for k in range(3, 48 + 1):
            print("Finding optimal size %d tour" % k)
            r.append(self.find_tour(k))
        return zip(*r)

    def write_results(self, r, path=tsp_path):
        for t, dist, names in zip(*r):
            n = len(t)
            n_s = str(n).zfill(2)
            with open(tsp_path + "tour_id_%s.txt" % n_s, 'wb') as f:
                f.write('\n'.join(map(str, t)))
            with open(tsp_path + "dist_%s.txt" % n_s, 'wb') as f:
                f.write('\n'.join(map(str, dist)))
            with open(tsp_path + "tour_names_%s.txt" % n_s, 'wb') as f:
                f.write('\n'.join(names))

