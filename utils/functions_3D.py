import numpy as np
import scipy as sp
import sympy as sym
from sympy import Ynm
from sympy.utilities.lambdify import lambdify
from sympy.parsing.mathematica import parse_mathematica
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import quaternion



# spherical angles θ and φ
sym_t = sym.Symbol("t", real=True)
sym_p = sym.Symbol("p", real=True)

# curvature kappa
kappa_expression = parse_mathematica('(-(2 rp rt Sin[t] + r0 (rp Cos[t] - rtp Sin[t]))^2 + (r0^2 + 2 rt^2 -  r0 rtt) Sin[t]^2 (2 rp^2 + r0^2 Sin[t]^2 -  r0 (rpp + rt Cos[t] Sin[t])))/(r0^2 (rp^2 + r0^2 Sin[t]^2 + rt^2 Sin[t]^2)^2)')

# fundamental forms
FundE_expression = parse_mathematica('r0^2 + rt^2')
FundF_expression = parse_mathematica('rp rt')
FundG_expression = parse_mathematica('rp^2 + r0^2 Sin[t]^2')

FundL_expression = parse_mathematica('-((r0 (r0^2 + 2 rt^2 - r0 rtt) Sin[t])/Sqrt[r0^2 (rp^2 + r0^2 Sin[t]^2 + rt^2 Sin[t]^2)])')
FundM_expression = parse_mathematica('(r0 (-2 rp rt Sin[t] + r0 (-rp Cos[t] + rtp Sin[t])))/Sqrt[r0^2 (rp^2 + r0^2 Sin[t]^2 + rt^2 Sin[t]^2)]')
FundN_expression = parse_mathematica('-(r0 Sin[t] (2 rp^2 + r0^2 Sin[t]^2 - r0 (rpp + rt Cos[t] Sin[t])))/Sqrt[r0^2 (rp^2 + r0^2 Sin[t]^2 + rt^2 Sin[t]^2)]')

# Q factor that determines when halting states disappear
QFact_expression = parse_mathematica('Sqrt[rt^2 + rp^2 Csc[t]^2]/r0')

# body frame contact vector r components
x0_expression = parse_mathematica('r0 Cos[p] Sin[t]')
y0_expression = parse_mathematica('r0 Sin[p] Sin[t]')
z0_expression = parse_mathematica('r0 Cos[t]')

# body frame dr/dθ components
xt_expression = parse_mathematica('r0 Cos[p] Cos[t] + rt Cos[p] Sin[t]')
yt_expression = parse_mathematica('r0 Cos[t] Sin[p] + rt Sin[p] Sin[t]')
zt_expression = parse_mathematica('rt Cos[t] - r0 Sin[t]')

# body frame dr/dφ components
xp_expression = parse_mathematica('rp Cos[p] Sin[t] - r0 Sin[p] Sin[t]')
yp_expression = parse_mathematica('r0 Cos[p] Sin[t] + rp Sin[p] Sin[t]')
zp_expression = parse_mathematica('rp Cos[t]')

INPUTS = [sym_t,sym_p, sym.Symbol("r0"), sym.Symbol("rt"), sym.Symbol("rtt"), sym.Symbol("rp"), sym.Symbol("rpp"), sym.Symbol("rtp")]
		
kappa_expression = lambdify(INPUTS, kappa_expression)

FundE_expression = lambdify(INPUTS, FundE_expression)
FundF_expression = lambdify(INPUTS, FundF_expression)
FundG_expression = lambdify(INPUTS, FundG_expression)

FundL_expression = lambdify(INPUTS, FundL_expression)
FundM_expression = lambdify(INPUTS, FundM_expression)
FundN_expression = lambdify(INPUTS, FundN_expression)

QFact_expression = lambdify(INPUTS, QFact_expression)

x0_expression = lambdify(INPUTS, x0_expression)
y0_expression = lambdify(INPUTS, y0_expression)
z0_expression = lambdify(INPUTS, z0_expression)

xt_expression = lambdify(INPUTS, xt_expression)
yt_expression = lambdify(INPUTS, yt_expression)
zt_expression = lambdify(INPUTS, zt_expression)

xp_expression = lambdify(INPUTS, xp_expression)
yp_expression = lambdify(INPUTS, yp_expression)
zp_expression = lambdify(INPUTS, zp_expression)

# EOM for viscous rolling
# omega vectors are functions of: q0 - q3, x0, y0, z0, F, and G
# F is force parallel to ramp, G is force perpendicular to ramp
omgx_expression = parse_mathematica('G (2 q1 q2 x0 + 2 q0 q3 x0 + q0^2 y0 - q1^2 y0 + q2^2 y0 - q3^2 y0 - 2 q0 q1 z0 + 2 q2 q3 z0)')
omgy_expression = parse_mathematica('-F (-2 q0 q2 x0 + 2 q1 q3 x0 + 2 q0 q1 y0 +    2 q2 q3 y0 + q0^2 z0 - q1^2 z0 - q2^2 z0 + q3^2 z0) - G (q0^2 x0 + q1^2 x0 - (q2^2 + q3^2) x0 + q0 (-2 q3 y0 + 2 q2 z0) +    2 q1 (q2 y0 + q3 z0))')
omgz_expression = parse_mathematica('F (2 q1 q2 x0 + 2 q0 q3 x0 + q0^2 y0 - q1^2 y0 + q2^2 y0 - q3^2 y0 - 2 q0 q1 z0 + 2 q2 q3 z0)')
omg_inputs = [sym.Symbol("q0"),
              sym.Symbol("q1"),
              sym.Symbol("q2"),
              sym.Symbol("q3"),
              sym.Symbol("x0"),
              sym.Symbol("y0"),
              sym.Symbol("z0"),
              sym.Symbol("F"),
              sym.Symbol("G")]

omgx_expression = lambdify(omg_inputs, omgx_expression)
omgy_expression = lambdify(omg_inputs, omgy_expression)
omgz_expression = lambdify(omg_inputs, omgz_expression)

# propagating q with omega vectors as inputs
Q0Dt_expression = parse_mathematica('1/2 (-omgx q1 - omgy q2 - omgz q3)')
Q1Dt_expression = parse_mathematica('1/2 (omgx q0 - omgz q2 + omgy q3)')
Q2Dt_expression = parse_mathematica('1/2 (omgy q0 + omgz q1 - omgx q3)')
Q3Dt_expression = parse_mathematica('1/2 (omgz q0 - omgy q1 + omgx q2)')
QDt_inputs = [sym.Symbol("q0"),
              sym.Symbol("q1"),
              sym.Symbol("q2"),
              sym.Symbol("q3"),
              sym.Symbol("omgx"),
              sym.Symbol("omgy"),
              sym.Symbol("omgz")]

Q0Dt_expression = lambdify(QDt_inputs, Q0Dt_expression)
Q1Dt_expression = lambdify(QDt_inputs, Q1Dt_expression)
Q2Dt_expression = lambdify(QDt_inputs, Q2Dt_expression)
Q3Dt_expression = lambdify(QDt_inputs, Q3Dt_expression)

# unit vector for convenience
xVec = np.array([1, 0, 0])
yVec = np.array([0, 1, 0])
zVec = np.array([0, 0, 1])

# general class for storing the 3D ball information
# also handles simulation
# state vector = [θ, φ, q0, q1, q2, q3, x, y, z] for viscous rolling
# state vector = [θ, φ, q0, q1, q2, q3, x, y, z, ωx, ωy, ωz] for inertia rolling
class ball:

	# at initialisation, generate radius function in spherical coordinate
	def __init__(self, maxL, As, Ps):
		
		sym_r = 1
		ct = 0
		for i in range(2, maxL + 1):
		    for j in np.arange(-i, i + 1):
		        sym_r += As[ct] * Ynm(i, j, sym_t, sym_p + Ps[ct]) / (i * (i + 1))
		        ct += 1
		        
		sym_r = sym.expand_trig(sym.re(sym_r).expand(func=True))
		sym_rt  = sym.diff(sym_r, sym_t)
		sym_rtt = sym.diff(sym_r, sym_t, sym_t)
		sym_rp  = sym.diff(sym_r, sym_p)
		sym_rpp = sym.diff(sym_r, sym_p, sym_p)
		sym_rtp = sym.diff(sym_r, sym_t, sym_p)	

		self.r0_expression  = lambdify([sym_t,sym_p], sym_r)
		self.rt_expression  = lambdify([sym_t,sym_p], sym_rt )
		self.rtt_expression = lambdify([sym_t,sym_p], sym_rtt)
		self.rp_expression  = lambdify([sym_t,sym_p], sym_rp )
		self.rpp_expression = lambdify([sym_t,sym_p], sym_rpp)
		self.rtp_expression = lambdify([sym_t,sym_p], sym_rtp)

	# useful function to get INPUTS (see above) from θ and φ
	def get_inputs(self, _t, _p):
		return [_t,_p,
				self.r0_expression(_t, _p),
				self.rt_expression(_t, _p),
				self.rtt_expression(_t, _p),
				self.rp_expression(_t, _p),
				self.rpp_expression(_t, _p),
				self.rtp_expression(_t, _p)]

	# render the ball
	# if keep == True, these data are used to visualise the ball
	# see draw_ball, draw_kappa, and draw_QFact below
	def make_ball(self, grid_point = 500, keep = False):

		theta, phi = np.linspace(1e-4, np.pi - 1e-4, grid_point), np.linspace(0, 2 * np.pi, grid_point)
		THETA, PHI = np.meshgrid(theta, phi)
		
		_inputs = self.get_inputs(THETA, PHI)

		X = x0_expression(_inputs[0], _inputs[1], _inputs[2], _inputs[3], _inputs[4], _inputs[5], _inputs[6], _inputs[7])
		Y = y0_expression(_inputs[0], _inputs[1], _inputs[2], _inputs[3], _inputs[4], _inputs[5], _inputs[6], _inputs[7])
		Z = z0_expression(_inputs[0], _inputs[1], _inputs[2], _inputs[3], _inputs[4], _inputs[5], _inputs[6], _inputs[7])
		
		kappa = kappa_expression(_inputs[0], _inputs[1], _inputs[2], _inputs[3], _inputs[4], _inputs[5], _inputs[6], _inputs[7])
		QFact = QFact_expression(_inputs[0], _inputs[1], _inputs[2], _inputs[3], _inputs[4], _inputs[5], _inputs[6], _inputs[7])

		if keep:
			self.THETA = THETA
			self.PHI = PHI 
			self.X = X 
			self.Y = Y 
			self.Z = Z 
			self.kappa = kappa 
			self.QFact = QFact

		return X, Y, Z, THETA, PHI, kappa, QFact

	# check the curvature across the ball
	def check_kappa(self, grid_point = 100):

		X, Y, Z, THETA, PHI, kappa, QFact = self.make_ball(grid_point = grid_point)
		self.kmin = np.min(kappa)
		self.kmax = np.max(kappa)
		
		print('Maximum and minimum curvatures are: {0:.2f} and {1:.2f}. These will be used to set upper and lower bounds on cmaps'.format(np.max(kappa), np.min(kappa)))
	
	# manually setting the colour map vmin/vmax
	def set_kmimmax(self, kmin, kmax):
		self.kmin = kmin
		self.kmax = kmax

	# draw the ball in 3D
	def draw_ball(self, ax, cmap = plt.cm.RdBu_r):

		fcolors = cmap((self.kappa - self.kmin) / (self.kmax - self.kmin))
		ax.plot_surface(self.X, self.Y, self.Z, facecolors = fcolors, rstride = 1, cstride = 1, antialiased=False)
		ax.set_box_aspect((1,1,1))
		ax.set_title('3D ball')
		
		ax.set_xlabel('x')
		ax.set_ylabel('y')
		ax.set_zlabel('z')

	# plot curvature
	def draw_kappa(self, ax, cmap = 'RdBu_r'):

		img = ax.imshow(self.kappa.T, cmap = cmap, origin = 'upper', vmin=self.kmin, vmax=self.kmax, extent = [0, 2 * np.pi, 1 * np.pi, 0])
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)
		plt.colorbar(img,cax=cax)
		ax.set_title('map of curvature')

	# plot Q factor
	def draw_QFact(self, ax, cmap = 'RdBu_r'):

		img = ax.imshow(self.QFact.T, cmap = cmap, origin = 'upper', extent = [0, 2 * np.pi, 1 * np.pi, 0])
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)
		plt.colorbar(img,cax=cax)
		ax.set_title('map of Q factor')

	# compute the quaternion that rotates the ball
	# such that the chosen (θ, φ) point is in contact with the ramp
	# and the ball is tangent to ramp
	def compute_direct_Q(self, _inputs):

		# first generate quaternion for the normal vector
		R0 = np.array([x0_expression(*_inputs),
					   y0_expression(*_inputs),
					   z0_expression(*_inputs)])

		Rt = np.array([xt_expression(*_inputs),
		               yt_expression(*_inputs),
		               zt_expression(*_inputs)])
		
		Rp = np.array([xp_expression(*_inputs),
		               yp_expression(*_inputs),
		               zp_expression(*_inputs)])
		
		EE = FundE_expression(*_inputs)
		FF = FundF_expression(*_inputs)
		GG = FundG_expression(*_inputs)
		
		NN = np.cross(Rt, Rp) / np.sqrt(EE * GG - FF ** 2)
		NN = np.quaternion(0, *NN)
		
		# target quaternion: pointing into the ramp
		target = - np.quaternion(0,0,0,1)
		
		# get the Q quaternion by combining the angle bisector and z axis
		q1 = (NN + target)
		q1 = q1 / np.sqrt(q1.norm())
		q2 = target

		Q = q2 * q1

		return(Q)
	
	# initialise the ball to sit on a given (θ, φ) point
	# the spin parameter further spins the ball on the spot
	# if spin == 'AUTO', the code attempts to spin the ball such that the initial Q has zero real component
	# i.e. the vector part lies on the surface of the 3-ball, in contrast to the interior
	def initialisation(self, THETA_CT = np.pi / 2, PHI_CT = 0, spin = 0, cmap = plt.cm.RdBu_r, DomegaDt = [0, 0, 0], inertia = False):

		_inputs = self.get_inputs(THETA_CT, PHI_CT)

		R0 = np.array([x0_expression(*_inputs),
					   y0_expression(*_inputs),
					   z0_expression(*_inputs)])

		Q = self.compute_direct_Q(_inputs)

		if spin == 'AUTO':
			_Q = quaternion.as_float_array(Q)
			if _Q[0] == 0:
				spin = 0 
			elif _Q[3] == 0:
				spin = np.pi 
			else:
				spin = 2 * np.arctan(_Q[0] / _Q[3])
		else:
			pass
		q3 = np.quaternion(np.cos(spin / 2), *(np.sin(spin / 2)) * np.array([0,0,1]))
		Q = q3 * Q

		CM_Pos = - Q * np.quaternion(0, *R0) * Q.conjugate()

		COOR = np.array([np.zeros(self.X.shape), self.X, self.Y, self.Z])
		COOR = COOR.reshape(4, -1).T
		COOR = quaternion.as_quat_array(COOR)
		COOR_T = Q * COOR * Q.conjugate()
		COOR_T = quaternion.as_float_array(COOR_T).T

		# starting position 3D ball
		self.X_T = COOR_T[1].reshape(self.THETA.shape)
		self.Y_T = COOR_T[2].reshape(self.THETA.shape)
		self.Z_T = COOR_T[3].reshape(self.THETA.shape)

		# draw a black spot on the ball that marks the initial contact point
		self.fcolors_T = cmap((self.kappa - self.kmin) / (self.kmax - self.kmin))
		CT_LOC = ((self.THETA - THETA_CT) ** 2 + np.sin(self.THETA) ** 2 * ((self.PHI - PHI_CT) % (np.pi * 2)) ** 2) < 0.01
		self.fcolors_T[CT_LOC] = [0, 0, 0, 1]
		if inertia:
			self.s0 = [THETA_CT, PHI_CT, *quaternion.as_float_array(Q), *quaternion.as_float_array(CM_Pos)[1:], *DomegaDt]
		else:
			self.s0 = [THETA_CT, PHI_CT, *quaternion.as_float_array(Q), *quaternion.as_float_array(CM_Pos)[1:]]

	# initialise the ball without computing some of the other things
	# the 'other things' are for checking the initialisation makes sense using draw_initial
	def initialisation_simple(self, THETA_CT = np.pi / 2, PHI_CT = 0, spin = 0, DomegaDt = [0, 0, 0], inertia = False):

		_inputs = self.get_inputs(THETA_CT, PHI_CT)

		R0 = np.array([x0_expression(*_inputs),
					   y0_expression(*_inputs),
					   z0_expression(*_inputs)])

		Q = self.compute_direct_Q(_inputs)

		if spin == 'AUTO':
			_Q = quaternion.as_float_array(Q)
			if _Q[0] == 0:
				spin = 0 
			elif _Q[3] == 0:
				spin = np.pi 
			else:
				spin = 2 * np.arctan(_Q[0] / _Q[3])
		else:
			pass
		q3 = np.quaternion(np.cos(spin / 2), *(np.sin(spin / 2)) * np.array([0,0,1]))
		Q = q3 * Q

		CM_Pos = - Q * np.quaternion(0, *R0) * Q.conjugate()

		COOR = np.array([np.zeros(self.X.shape), self.X, self.Y, self.Z])
		COOR = COOR.reshape(4, -1).T
		COOR = quaternion.as_quat_array(COOR)
		COOR_T = Q * COOR * Q.conjugate()
		COOR_T = quaternion.as_float_array(COOR_T).T

		if inertia:
			self.s0 = [THETA_CT, PHI_CT, *quaternion.as_float_array(Q), *quaternion.as_float_array(CM_Pos)[1:], *DomegaDt]
		else:
			self.s0 = [THETA_CT, PHI_CT, *quaternion.as_float_array(Q), *quaternion.as_float_array(CM_Pos)[1:]]

	def draw_initial(self, ax1, ax2):

		if ax1:

			ax1.plot_surface(self.X, self.Y, self.Z, facecolors = self.fcolors_T, rstride = 1, cstride = 1, antialiased=False)
			ax1.set_xlabel('x')
			ax1.set_ylabel('y')
			ax1.set_zlabel('z')
			ax1.set_title('body frame')
			ax1.set_box_aspect((1,1,1))

		if ax2:
			ax2.plot_surface(self.X_T, self.Y_T, self.Z_T, facecolors = self.fcolors_T, rstride = 1, cstride = 1, antialiased=False)
			ax2.set_xlabel('x')
			ax2.set_ylabel('y')
			ax2.set_zlabel('z')
			ax2.set_title('lab frame')
			ax2.set_box_aspect((1,1,1))

	# function that implements the viscous EOM
	# f is force parallel to ramp, g is perpendicular to ramp
	def propagate_viscous(self, time, state, f, g):
	
		[_t, _p, _q0, _q1, _q2, _q3, _x, _y, _z] = state
		
		_inputs = self.get_inputs(_t, _p)

		_x0 = x0_expression(*_inputs)
		_y0 = y0_expression(*_inputs)
		_z0 = z0_expression(*_inputs)
		
		_xt = xt_expression(*_inputs)
		_yt = yt_expression(*_inputs)
		_zt = zt_expression(*_inputs)
		
		_xp = xp_expression(*_inputs)
		_yp = yp_expression(*_inputs)
		_zp = zp_expression(*_inputs)
		
		_r0_body = np.array([_x0, _y0, _z0])
		_rt_body = np.array([_xt, _yt, _zt])
		_rp_body = np.array([_xp, _yp, _zp])
		
		_q = np.quaternion(_q0, _q1, _q2, _q3)
		
		_r0_Labf = _q * np.quaternion(0, *_r0_body) * _q.conjugate()
		_rt_Labf = _q * np.quaternion(0, *_rt_body) * _q.conjugate()
		_rp_Labf = _q * np.quaternion(0, *_rp_body) * _q.conjugate()
		
		_r0_Labf = quaternion.as_float_array(_r0_Labf)[1:]
		_rt_Labf = quaternion.as_float_array(_rt_Labf)[1:]
		_rp_Labf = quaternion.as_float_array(_rp_Labf)[1:]
		
		_omgx = omgx_expression(_q0, _q1, _q2, _q3, _x0, _y0, _z0, f, g)
		_omgy = omgy_expression(_q0, _q1, _q2, _q3, _x0, _y0, _z0, f, g)
		_omgz = omgz_expression(_q0, _q1, _q2, _q3, _x0, _y0, _z0, f, g)
		
		_omegave = [_omgx, _omgy, _omgz]
		
		_Q0Dt = Q0Dt_expression(_q0, _q1, _q2, _q3, _omgx, _omgy, _omgz)
		_Q1Dt = Q1Dt_expression(_q0, _q1, _q2, _q3, _omgx, _omgy, _omgz)
		_Q2Dt = Q2Dt_expression(_q0, _q1, _q2, _q3, _omgx, _omgy, _omgz)
		_Q3Dt = Q3Dt_expression(_q0, _q1, _q2, _q3, _omgx, _omgy, _omgz)
		
		_L = FundL_expression(*_inputs)
		_M = FundM_expression(*_inputs)
		_N = FundN_expression(*_inputs)
		
		_mat1 = np.array([[ _N, -_M],
		                  [-_M, _L]]) / (_L * _N - _M ** 2)
		
		_mat2 = np.array([np.dot(_rt_Labf, np.cross(_omegave, -zVec)),
		                  np.dot(_rp_Labf, np.cross(_omegave, -zVec))])
		
		_DDtp = np.matmul(_mat1,_mat2)
		
		_v = - np.cross(_omegave, _r0_Labf)
		
		return np.array([_DDtp[0], _DDtp[1],
		                 _Q0Dt, _Q1Dt, _Q2Dt, _Q3Dt,
		                 _v[0], _v[1], _v[2]])

	def propagate_inertia(self, time, state, f, g, m):

		[_t, _p, _q0, _q1, _q2, _q3, _x, _y, _z, _omgx, _omgy, _omgz] = state
		
		_inputs = np.array(self.get_inputs(_t, _p)).astype(float)

		_r0 = _inputs[2]

		_x0 = x0_expression(*_inputs)
		_y0 = y0_expression(*_inputs)
		_z0 = z0_expression(*_inputs)
		
		_xt = xt_expression(*_inputs)
		_yt = yt_expression(*_inputs)
		_zt = zt_expression(*_inputs)
		
		_xp = xp_expression(*_inputs)
		_yp = yp_expression(*_inputs)
		_zp = zp_expression(*_inputs)
		
		_r0_body = np.array([_x0, _y0, _z0])
		_rt_body = np.array([_xt, _yt, _zt])
		_rp_body = np.array([_xp, _yp, _zp])
		
		_q = np.quaternion(_q0, _q1, _q2, _q3)
		
		_r0_Labf = _q * np.quaternion(0, *_r0_body) * _q.conjugate()
		_rt_Labf = _q * np.quaternion(0, *_rt_body) * _q.conjugate()
		_rp_Labf = _q * np.quaternion(0, *_rp_body) * _q.conjugate()
		
		_r0_Labf = quaternion.as_float_array(_r0_Labf)[1:]
		_rt_Labf = quaternion.as_float_array(_rt_Labf)[1:]
		_rp_Labf = quaternion.as_float_array(_rp_Labf)[1:]
		
		_omegave = np.array([_omgx, _omgy, _omgz])
		
		_Q0Dt = Q0Dt_expression(_q0, _q1, _q2, _q3, _omgx, _omgy, _omgz)
		_Q1Dt = Q1Dt_expression(_q0, _q1, _q2, _q3, _omgx, _omgy, _omgz)
		_Q2Dt = Q2Dt_expression(_q0, _q1, _q2, _q3, _omgx, _omgy, _omgz)
		_Q3Dt = Q3Dt_expression(_q0, _q1, _q2, _q3, _omgx, _omgy, _omgz)
		
		_L = FundL_expression(*_inputs)
		_M = FundM_expression(*_inputs)
		_N = FundN_expression(*_inputs)
		
		_mat1 = np.array([[ _N, -_M],
		                  [-_M, _L]]) / (_L * _N - _M ** 2)
		
		_mat2 = np.array([np.dot(_rt_Labf, np.cross(_omegave, -zVec)),
		                  np.dot(_rp_Labf, np.cross(_omegave, -zVec))])
		
		_DDtp = np.matmul(_mat1,_mat2)
		
		_v = - np.cross(_omegave, _r0_Labf)
		
		# calculate rate of change of angular velocity
		_DrDt = np.cross(_omegave, _r0_Labf) + _rt_Labf * _DDtp[0] + _rp_Labf * _DDtp[1]
		_term1 = m * np.cross(_r0_Labf, np.cross(_DrDt, _omegave))
		_term2 = - 5 /2 * _r0_Labf * np.dot(_r0_Labf, _omegave)
		_term3 = - np.cross(_r0_Labf, [f, 0, - g])
		_term4 = - _omegave
		_DomDt = (_term1 + _term2 + _term3 + _term4) / (2 * m / 5 + m * _r0 ** 2)
		
		return np.array([_DDtp[0], _DDtp[1],
		                 _Q0Dt, _Q1Dt, _Q2Dt, _Q3Dt,
		                 _v[0], _v[1], _v[2], _DomDt[0], _DomDt[1], _DomDt[2]])

	# simulation
	def simulate(self, f = 0, g = 0, ti = 0, tf = 3, outputN = 1000, m = 0):
		if m == 0:
			res = sp.integrate.solve_ivp(self.propagate_viscous, (ti, tf), self.s0, t_eval = np.linspace(ti, tf, outputN), args = (f, g), method = 'DOP853')
			self.s = res.y.T
			self.t = res.t
		else:
			res = sp.integrate.solve_ivp(self.propagate_inertia, (ti, tf), self.s0, t_eval = np.linspace(ti, tf, outputN), args = (f, g, m), method = 'DOP853')
			self.s = res.y.T
			self.t = res.t
			
	# check that the error did not accummulate
	def results_quality_control(self, axes):
		
		self.CM_Z = [self.s0[8]]
		
		self.CT_X = [0]
		self.CT_Y = [0]
		
		self.nn = [[0, 0, 0, -1]]
		self.normQ = [1]
		
		for i in range(len(self.s) - 1):
		    i = i + 1
		    
		    _t = self.s[i][0]
		    _p = self.s[i][1]
		    _q = np.quaternion(*self.s[i][2:6])
		    _x = self.s[i][6]
		    _y = self.s[i][7]

		    _inputs = self.get_inputs(_t, _p)
		    
		    _R0 = np.array([x0_expression(*_inputs),
		                    y0_expression(*_inputs),
		                    z0_expression(*_inputs)])
		    
		    _R0 = np.quaternion(0, *_R0)
		    _R0 = quaternion.as_float_array(_q * _R0 * _q.conjugate())[1:]
		    
		    self.CM_Z.append(- _R0[2])
		    
		    self.CT_X.append(_x + _R0[0])
		    self.CT_Y.append(_y + _R0[1])
		    
		    _Rt = np.array([xt_expression(*_inputs),
		                    yt_expression(*_inputs),
		                    zt_expression(*_inputs)])
		    
		    _Rp = np.array([xp_expression(*_inputs),
		                    yp_expression(*_inputs),
		                    zp_expression(*_inputs)])
		    
		    EE = FundE_expression(*_inputs)
		    FF = FundF_expression(*_inputs)
		    GG = FundG_expression(*_inputs)
		    
		    NN = np.cross(_Rt, _Rp) / np.sqrt(EE * GG - FF ** 2)
		    NN = np.quaternion(0, *NN)
		    NN = _q * NN * _q.conjugate()
		    
		    self.nn.append(quaternion.as_float_array(NN))
		    self.normQ.append(_q.norm())
		
		self.CM_Z = np.array(self.CM_Z)
		
		self.CT_X = np.array(self.CT_X)
		self.CT_Y = np.array(self.CT_Y)

		self.normQ = np.array(self.normQ)
		self.nn = np.array(self.nn)
		
		for i in range(3):
			axes[i].plot(self.t, self.nn[:,1 + i])
			axes[i].set_title('qnq* component {0}'.format(i + 1))

		axes[3].plot(self.t, self.normQ)
		axes[3].set_title('Q norm')

		axes[4].set_title('z position')
		axes[4].plot(self.t, self.CM_Z, label = 'holonomic')
		axes[4].plot(self.t, self.s[:,8], label = 'non-slip')
		axes[4].legend()

		axes[5].plot(self.CT_X, self.CT_Y, label = 'contact')
		axes[5].plot(self.s[:, 6], self.s[:, 7], label = 'c. o. mass')
		axes[5].set_title('xy trajectory')
		axes[5].legend()
					
					