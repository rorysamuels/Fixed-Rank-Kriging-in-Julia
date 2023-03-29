# Run from Julia_FRK directory

# This is a short walkthrough of the implemented FRK functions using simiulated data over CONUS. Tuning parameters have been chosen arbitrarily so no seed was set.

# Load required packages and include the FRK modules.
using NetCDF
using Distributions
include("../FRK.jl")



# Generate data
n = 1500
lat = rand(n) * 50
lon = rand(n) * -60
U = hcat(lat,lon)
Z = rand(Normal(270,7),n)
X = lat


# We can first fit the model using FRK_fit().
# We need to pass the data locations, response, and covariates (in that order)
# followed by a vector of arguments to control the basis function construction.
# This vector of arguments will be passed to the general function Create_basis().
# We can pass any vector that would result in a valid Create_basis() method, for
# example:

my_basis_args = [120,3]

# will generate call Create_basis(U, lowres_centers::Integer, n_res::Integer).
# Although we could  have defined a boundary for the basis functions to be
# generated within,  here I am specifying that I want 115 basis functions
# (at the coarsest resolution) to be generated globally, with 2 subsequent
# resolutions. This sounds like a lot of basis functions but Create_basis()
# will automaically remove basis functions  which evaluate near zero
# for every observation.

# We pass these arguments to FRK_fit which returns a list of model information.
model = FRK_fit(U,Z,X, basis_args = my_basis_args);


# We can pass the fitted model to FRK_predict() to make predicitons and
# calculate uncertainity estiamtes. In addition to the list outputted from FRK_fit(),
# we also need to specify arguements to control the construction of the BAUs.
# We can pass any vector that would result in a valid BAU_grid() method.
# For example, here I will specify a resolution, L_P, and a prediction region.

L_P = 5; lat_range = [33,50]; lon_range = [-118,-80]

my_BAU_args = [L_P, lat_range, lon_range]

# Below, I have also specified a few additional arguemnts which are optional.
# I have asked for the mspe to be returned and have set parallel = true,
# which will distribute the computation onto the number of prcessors according
# to Threads.nthreads() (this can be set with the Julia environment variable
# JULIA_NUM_THREADS or by running julia with the option julia --threads).


results = FRK_predict(model, BAU_args = my_BAU_args, parallel = true, mspe = true);

results

# Alternatively, we can use the function FRK() to accomplish everything in one shot.
# The down side is we don't have a 'model object' if we want to make more predictions
# at new locations or at a different resolution for example. In this example,
# I specify that I want the results saved in the directory 'exampeles/results'.

FRK(U,Z,X, basis_args = my_basis_args, BAU_args = my_BAU_args, parallel=true,
    mspe = true, save = true, path = "examples/results", dataname = "AIRS_NST",
    print = false)

