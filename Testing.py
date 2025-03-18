# Part 2
import math
import matplotlib.pyplot as plot

# f(x) = e^x
def f(x):
    return math.exp(x)

# Compute the divided differences table by the lower triangle
def divided_diff(x, y):
    n = len(x)
    table = [[0] * n for _ in range(n)]
    
    # First column is just y-values, f(x_i)
    for i in range(n):
        table[i][0] = y[i]
    
    # Fill the lower triangle of the matrix
    for j in range(1, n):
        for i in range(j, n):
            table[i][j] = (table[i][j-1] - table[i-1][j-1]) / (x[i] - x[i-j])
    
    # The diagonal of the lower triangle contains the Newton coefficients
    return [table[i][i] for i in range(n)]

# Evaluate the Newton polynomial at point x_val
def newton_polynomial(x, coeff, x_val):
    n = len(coeff)
    result = coeff[0]
    product = 1.0
    
    # Defining each term in the Newton polynomial
    for i in range(1, n):
        product *= (x_val - x[i - 1])
        
        # Computing the y_val given x_val
        result += coeff[i] * product
    return result

# Compute the maximum error for a set of n + 1 points
def max_error(n):
    # Equally spaced points in the interval [-1, 1]
    xList = [-1 + 2 * i / n for i in range(n + 1)] # For n degree + 1
    yList = [f(x) for x in xList]
    
    # Get divided differences coefficients
    coeff = divided_diff(xList, yList)
    
    # Sample 501 points from [-1, 1]
    t_values = [-1 + 2 * i / 500 for i in range(501)]
    max_error = 0
    
    # Find max error
    for t in t_values:
        p = newton_polynomial(xList, coeff, t)
        error = abs(f(t) - p)
        if error > max_error:
            max_error = error
    return max_error

# Compute errors for different n values and plot
n_values = [2, 4, 8, 16, 32]
errors = []
for n in n_values:
    errors.append(max_error(n))
    print("E(" + str(n)+ ") = " + str(max_error(n)))

# Plot maximum error vs n
plot.plot(n_values, errors, marker='X')
plot.xlabel('n nodes')
plot.ylabel('Max Error, E(n)')
plot.title('Max Error vs n for Newton Interpolation of e^x')
plot.show()