import numpy as np
import copy
import os
import warnings
cwd = os.getcwd()

# Retrieve simulation parameters

simulation = np.genfromtxt('../temp-data/simulationParameters.csv', delimiter=',')

simulationParameters = simulation[:, 1]

NSITES = int(simulationParameters[0])
dt = simulationParameters[1]
beta = simulationParameters[2]
L = int(simulationParameters[3])
t = simulationParameters[4]
U = simulationParameters[5]
mu = simulationParameters[6]
totalMCSweeps = int(simulationParameters[7])
freq = int(simulationParameters[8])
intsize = int(simulationParameters[9])
geom = int(simulationParameters[10])
ny = int(simulationParameters[11])

# Load configuration weights, etc.

weights = np.loadtxt('../temp-data/Log-weights.csv', skiprows = 1)
nEl = np.loadtxt('../temp-data/Electron-density.csv', skiprows = 1)
nUp_nDw = np.loadtxt('../temp-data/Double-occupancy.csv', skiprows = 1)
magCorr = np.loadtxt('../temp-data/Sz10.csv', skiprows = 1)
uneqMagCorr = np.loadtxt('../temp-data/UneqSz10.csv', skiprows = 1)
signs = np.loadtxt('../temp-data/Signs.csv', skiprows = 1)

directory1 = ('data/' + str(NSITES) + \
              'sites_L=' + str(L) + \
              '_beta=' + str(beta) + \
              '_dt_' + str(dt) + '_t_' + \
              str(t) + '_U_'+ str(U) + '_mu_' + str(mu))

directory2 = (directory1 + '/data-to-reproduce/' + \
              'totalMCSweeps_' + str(totalMCSweeps) + \
              '_freq_' + str(freq) + '_intsize_' + str(intsize) + \
              '_geom_' + str(geom) + '_ny_' + str(ny) )

directory3 = (directory1 + '/plots/' + \
              'totalMCSweeps_' + str(totalMCSweeps) + \
              '_freq_' + str(freq) + '_intsize_' + str(intsize) + \
              '_geom_' + str(geom) + '_ny_' + str(ny) )


if not os.path.exists(directory1):
    os.makedirs(directory1)

if not os.path.exists(directory2):
    os.makedirs(directory2)

if not os.path.exists(directory3):
    os.makedirs(directory3)

np.savetxt(directory2 + '/Log-weights.csv', (weights), header = "Configuration log weight")
np.savetxt(directory2 + '/Electron-density.csv', (nEl), header = "Electron density <n>")
np.savetxt(directory2 + '/Double-occupancy.csv', (nUp_nDw), header = "Double occupancy <n+ n->,")
np.savetxt(directory2 + '/Sz10.csv', (magCorr), header = "<Sz_0 Sz_1 >")
np.savetxt(directory2 + '/UneqSz10.csv', (uneqMagCorr), header = "int_0^beta dt <Sz_0 (0) Sz_1 (1) >")
np.savetxt(directory2 + '/Signs.csv', (signs), header = "Average sign")

paramNames  = np.array(['Number of sites,', 'dtau,', 'beta,', 'L,',\
't,', 'U,', 'mu,', 'Number of MC Sweeps,', 'Frequency of recomputing G,',\
'Number of multiplied Bs after stabilization,', 'Geometry,', 'Ny,'])
namedParams = np.zeros(paramNames.size, dtype=[('pnames', 'U50'), ('params', float)])
namedParams['pnames'] = paramNames
namedParams['params'] = simulationParameters

np.savetxt(directory2 + '/simulationParameters.csv', namedParams, fmt="%50s% 10.7f" )
