import numpy as np
from scipy.optimize import linprog

# --------------------------
# Input Parameters
# --------------------------
# cost = [a, b, c] for a + bP + cP²
cost = np.array([
    [561, 7.92, 0.001562],
    [310, 7.85, 0.00194],
    [78,  7.97, 0.00482]
])

Glim = np.array([
    [150, 600],
    [100, 400],
    [50,  200]
])

numSegment = int(input("Input number of segments: "))
Pload = 850
num_of_gen = 3

# --------------------------
# Build LP constraints
# --------------------------
Aeq = np.ones((1, numSegment * cost.shape[0]))
beq = np.array([Pload - np.sum(Glim[:, 0])])

Pedge = np.zeros((cost.shape[0], numSegment + 1))
slope = np.zeros((cost.shape[0], numSegment))
delPgen = np.zeros((cost.shape[0], numSegment))

# --------------------------
# Segment slopes & lengths
# --------------------------
for i in range(cost.shape[0]):
    Pedge[i, :] = np.linspace(Glim[i, 0], Glim[i, 1], numSegment + 1)
    for j in range(numSegment):
        limCost = cost[i, :] @ np.array([
            [1, 1],
            [Pedge[i, j+1], Pedge[i, j]],
            [Pedge[i, j+1] ** 2, Pedge[i, j] ** 2]
        ])
        delCost = limCost[0] - limCost[1]
        delP = Pedge[i, j+1] - Pedge[i, j]

        delPgen[i, j] = delP
        slope[i, j] = delCost / delP

# --------------------------
# Flatten slope & delPgen
# --------------------------
fmin = slope.flatten()
b = np.concatenate([delPgen.flatten(), np.zeros_like(delPgen.flatten())])

# --------------------------
# Inequality Constraints A x ≤ b
# --------------------------
A = np.zeros((2 * delPgen.size, delPgen.size))

# Upper bounds
for i in range(delPgen.size):
    A[i, i] = 1

# Lower bounds (-x ≤ 0)
for i in range(delPgen.size):
    A[i + delPgen.size, i] = -1

# --------------------------
# Solve LP using SciPy
# --------------------------
print("Solving using SciPy linprog...")
res = linprog(
    c=fmin,
    A_ub=A,
    b_ub=b,
    A_eq=Aeq,
    b_eq=beq,
    bounds=[(0, None)] * fmin.size,
    method="highs"
)

if not res.success:
    print("Optimization failed:", res.message)
    exit()

P = res.x

# --------------------------
# Compute Generator Outputs
# --------------------------
Pgen = np.zeros(num_of_gen)
idx = 0
for i in range(num_of_gen):
    Pgen[i] = np.sum(P[idx:idx + numSegment]) + Glim[i, 0]
    idx += numSegment

print("\nGenerator Outputs (MW):")
print("P1 =", Pgen[0])
print("P2 =", Pgen[1])
print("P3 =", Pgen[2])

# --------------------------
# Compute Cost
# --------------------------
nHours = 1
sub_Hourly_Cost = np.zeros((nHours, num_of_gen))

for h in range(nHours):
    for g in range(num_of_gen):
        Pg = Pgen[g]
        a, b_coef, c = cost[g, :]
        sub_Hourly_Cost[h, g] = c * Pg**2 + b_coef * Pg + a

Hourly_Cost = np.sum(sub_Hourly_Cost, axis=1)
totalcost = np.sum(Hourly_Cost)

print("\nTotal Cost of Operation =", totalcost)
