from skopt import gp_minimize


def objective(params):
    # Set the parameters for HypoDD
    # Run HypoDD
    # Compute the score based on your criteria (e.g., accuracy of locations)
    return score

res = gp_minimize(objective,                  # the function to minimize
                  [(value1, value2),          # the bounds on the first parameter
                   (value3, value4)],         # the bounds on the second parameter
                  acq_func="EI",              # the acquisition function
                  n_calls=100)                # the number of evaluations of f

print("Best parameters: {}".format(res.x))