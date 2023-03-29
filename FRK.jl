# Run from Julia_FRK directory

# Load necessary packages
using StatsBase
using Statistics
using LinearAlgebra
using DataFrames
using SparseArrays
using BlockArrays
using H3.API
using CSV
using Pipe



# Source scripts containing functions for FRK
modules = readdir("functions")

for script in modules
	include("functions/$script")
end
