# material selection model

import pyomo.environ as pyo

# data
materials = ['steel', 'alum', 'carbon', 'cheese']

density =   {   'steel' : 1.2,
                'alum'  : 0.8,
                'carbon': 1.8,
                'cheese': 0.7}

conductivity = {'steel' : 6.4,
                'alum'  : 3.1,
                'carbon': 4.4,
                'cheese': 0.3}

price =     {   'steel' : 2.3,
                'alum'  : 3.5,
                'carbon': 5.8,
                'cheese': 6.0}

lengths = [1,2,3,5]

m = pyo.ConcreteModel('material selector')

# SETS (used to index the decision variable and the parameters)
m.matl = pyo.Set(initialize=materials)

# VARIABLES
m.x = pyo.Var(m.matl, domain=pyo.Binary)   # a binary decision variable representing the selection of matl

# PARAMETERS
m.density = pyo.Param(m.matl, initialize=density)
m.conductivity = pyo.Param(m.matl, initialize=conductivity)
m.price = pyo.Param(m.matl, initialize=price)


# OBJ (minimize price)
m.obj = pyo.Objective(expr=sum(m.x[i] * m.price[i] for i in m.matl))

# Constraints
m.c1 = pyo.Constraint(expr=(sum(m.x[i] * m.density[i] for i in m.matl) >= 1.0))     # min density
m.c2 = pyo.Constraint(expr=(sum(m.x[i] * m.conductivity[i] for i in m.matl) <= 5.0)) # max cond.

# solve it
solver = pyo.SolverFactory('glpk')
result = solver.solve(m)
m.display()