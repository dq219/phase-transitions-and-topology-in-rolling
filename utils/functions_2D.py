import warnings
warnings.filterwarnings('ignore')

import numpy as np
import sympy as sym
import pandas as pd
from sympy.utilities.lambdify import lambdify
from sympy.parsing.mathematica import parse_mathematica

# get mathematica expressions for the various functions

# r0 - r(θ); r1 - r'(θ); r2 - r''(θ); r3 - r'''(θ)
r_inputs = [sym.symbols('r0'), sym.symbols('r1'), sym.symbols('r2'), sym.symbols('r3')]

# curvature kappa
kappa_expression = parse_mathematica('(r0^2 + 2 r1^2 - r0 r2)/(r0^2 + r1^2)^(3/2)')

# differential of kappa (with respect to θ)
DkaDt_expression = parse_mathematica('(3 r0^2 r1 r2 - 3 r1^3 r2 - r0^3 (r1 + r3) - r0 r1 (4 r1^2 - 3 r2^2 + r1 r3))/(r0^2 + r1^2)^(5/2)')

# arc length
arcLn_expression = parse_mathematica('Sqrt[r0^2 + r1^2]')

# square of the arc length
aLnSq_expression = parse_mathematica('r0^2 + r1^2')

kappa_expression = lambdify(r_inputs, kappa_expression)
DkaDt_expression = lambdify(r_inputs, DkaDt_expression)
arcLn_expression = lambdify(r_inputs, arcLn_expression)
aLnSq_expression = lambdify(r_inputs, aLnSq_expression)

# d^2 θ/dt^2 expression, used for inertia rolling
DtDt2_expression = parse_mathematica('(-theta1 (k1 m (1 + 2 r0^2) (r0^2 + r1^2) theta1 + k0 (4 m r0^3 r1 theta1 + m r0 r1 (1 + 2 r1^2) theta1 + r1 (2 r1 + m r2 theta1) + 2 r0^2 (1 + m r1 r2 theta1))) - 2 r0 r1 Cos[a] + 2 r0^2 Sin[a])/(k0 m (1 + 2 r0^2) (r0^2 + r1^2))')

DtDt2_inputs = [sym.symbols('r0'), sym.symbols('r1'), sym.symbols('r2'), sym.symbols('r3'), sym.symbols('k0'), sym.symbols('k1'), sym.symbols('theta1'), sym.symbols('m'), sym.symbols('a')]
DtDt2_expression = lambdify(DtDt2_inputs, DtDt2_expression)

# generate radius function using Fourier modes
def get_r(seed, num_modes = 5, std = 0.10):

	np.random.seed(seed) # set seed
	An = np.random.normal(0, std, num_modes) # random amplitudes
	dθ = np.random.random(num_modes) * np.pi * 2 # random phase
	
	# function r that takes in the polar angle θ as argument
	# 'order' is the order of differential to take
	# order = 0: r(θ)
	# order = 1: r'(θ), etc.
	def r(order, θ):

	    if order == 0:
	        _ = θ * 0 + 1 # baseline of 1
	    else:
	        _ = θ * 0

	    for i in range(len(An)):
	        l = (i + 2) # 
	        if (order%2) == 0:
	            _ += (-1) ** (order // 2) * An[i] * np.sin(l * (θ - dθ[i])) * l ** order / l ** 2
	        else: 
	            _ += (-1) ** ((order - 1) // 2) * An[i] * np.cos(l * (θ - dθ[i])) * l ** order / l ** 2
	    return _

	return r, An

# test whether a realisation is convex
def test_r(r, θ = np.linspace(0, 2 * np.pi, 1000), tag = None):

	R0 = r(0, θ)
	R1 = r(1, θ)
	R2 = r(2, θ)
	R3 = r(3, θ)
	
	Rs = [R0, R1, R2, R3]
	
	# calculate curvature and elementary arc length
	kappa = kappa_expression(*Rs)
	
	if min(kappa) < 0.1:
		if tag:
			print('The realisation {0} is nearly concave! Bad!'.format(tag))
		else:
			print('The realisation is nearly concave! Bad!')
		return False

	return True

# inverse contact point angular velocity in viscous rolling
# integration over this leads to the period
def dtdθ(θ, r, alpha):
	
	R0 = r(0, θ)
	R1 = r(1, θ)
	R2 = r(2, θ)
	R3 = r(3, θ)
	
	Rs = [R0, R1, R2, R3]
	
	kappa = kappa_expression(*Rs)
	
	top = (R0 ** 2 + R1 ** 2) * kappa
	bot = R0 ** 2 * np.sin(alpha) - R0 * R1 * np.cos(alpha)
	
	return top / bot

# equation of motion for inertia rolling
# t0 - θ(t); t1 - θ'(t); t2 - θ''(t)
# we have to define this after r(θ) is generated, so we can use solve_ivp of scipy
def propagate(time, state, a, m, r):
	[θ0, θ1] = state

	R0 = r(0, θ0)
	R1 = r(1, θ0)
	R2 = r(2, θ0)
	R3 = r(3, θ0)
	Rs = [R0, R1, R2, R3]

	k0 = kappa_expression(*Rs)
	k1 = DkaDt_expression(*Rs)
	θ2 = DtDt2_expression(*Rs, k0, k1, θ1, m, a)

	return np.array([θ1, θ2])

# fitting function for inertia rolling
def fit_fn(time, crit_time, slope, offset):
    mask = time > crit_time
    yy = np.ones(len(time)) * offset
    yy[mask] += slope * (time[mask] - crit_time)
    return yy
