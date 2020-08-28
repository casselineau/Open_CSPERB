import numpy as N

class HC():
	def __init__(self, T_min, T_max, force_T_limits=False):
		self.Tmin = T_min
		self.Tmax = T_max
		self.force_T_limits = force_T_limits

	@staticmethod
	def check_valid(prop):
		def wrapper_check_valid(self, T):
			if hasattr(T,'__len__'):
				Tlow = (T<self.Tmin).any()
				Thigh = (T>self.Tmax).any()
			else:
				Tlow = (T<self.Tmin)
				Thigh = (T>self.Tmax)
			if self.force_T_limits == False:
				if Tlow or Thigh:
					print "Temperature outside of correlations range"
					return N.nan
			return prop(self, T)
		return wrapper_check_valid


class Na(HC):
	'''	
	Fink JK, Leibowitz L. Thermodynamic and Transport Properties of Sodium Liquid and Vapour. Reactor Engineering Division, Argonne National Laboratory; 1995.
	'''
	def __init__(self, force_T_limits=False, force_Re_limits=False, force_Pe_limits=False, htc_correlation='Skupinski'):
		print 'Na HC'
		Tmin = 371.
		Tmax = 1255.
		HC.__init__(self, Tmin, Tmax, force_T_limits)
		self.force_Re_limits = force_Re_limits
		self.force_Pe_limits = force_Pe_limits
		self.htc_correlation = htc_correlation

	@HC.check_valid
	def rho(self, T):
		return 219.+275.32*(1.-T/2503.7)+511.58*N.sqrt(1.-T/2503.7)

	@HC.check_valid
	def mu(self, T):
		return N.exp(-6.4406-0.3958*N.log(T)+556.835/T)

	@HC.check_valid
	def k(self, T):
		return 124.67-0.11382*T+5.5226e-5*T**2.-1.1842e-8*T**3.

	@HC.check_valid
	def Cp(self, T):
		return (1.6582-8.4790e-4*T+4.4541e-7*T**2.-2992.6*T**(-2.))*1e3
	@HC.check_valid
	def h(self, T):

		return (-365.77+1.6582*T-4.2395e-4*T**2.+1.4847e-7*T**3.+2992.6/T)*1e3

	@HC.check_valid
	def s(self,T):
		return (1.6582*N.log(T)-8.4790e-4*T+4.4541e-7/2.*T**2.+2992.6/2.*T**(-2.))*1e3

	def V_tube(self, mf_HC_tube, T, D_tube_in):
		rho = self.rho(T)
		# Velocity in the tubes:
		V = mf_HC_tube/(N.pi*(D_tube_in/2.)**2.*rho)
		return V

	def h_conv_tube(self, mf_HC_tube, T, T_w, D_tube_in):
		rho = self.rho(T)
		mu = self.mu(T)
		Cp = self.Cp(T)
		k = self.k(T)
		# Velocity in the tubes:
		V = mf_HC_tube/(rho*N.pi*(D_tube_in/2.)**2.)
		# Reynolds:
		Re = D_tube_in*rho*V/mu
		# Prandtl:
		Pr = Cp*mu/k

		if self.htc_correlation == 'Lyon-Martinelli':
			if self.force_Re_limits == False:
				if ~(N.logical_and((Re>=1e4), (Re<=1e6)).all()):
					print 'Re:', Re
					print 'V:',V, '[m.s-1]'
					print 'Liquid sodium reynolds number outside of Lyon-Martinelli correlation range (1e4 <= Re <= 1e6)'
					return N.nan

			if (Pr<=0.1).all():
				# Nusselt:
				Nu = 7.+0.025*(Re*Pr)**0.8	
			else:
				print 'Pr:',Pr
				print 'Liquid sodium Prandtl number outside of Lyon-Martinelli correlation range(Pr <= 0.1)'
				return N.nan

		if self.htc_correlation == 'Skupinski':
			Pe = Re*Pr
			if self.force_Pe_limits == False:
				if ~(N.logical_and((Pe>=5.8e1), (Pe<=1.31e4)).all()):
					print 'Pe:', Pe
					print 'Liquid sodium Peclet number outside of Skupinski correlation range(5.8e1 <= Pe <= 1.31e4)'
					return N.nan
			if self.force_Re_limits == False:
				if ~(N.logical_and((Re>=3.6e3), (Re<=90.5e5)).all()):
					print 'Re:', Re
					print 'V:',V, '[m.s-1]'
					print 'Liquid sodium reynolds number outside of Skupinski correlation range (3.6e3 <= Re <= 90.5e5)'
					return N.nan

			Nu = 4.82+0.0185*Pe**0.827
		# Heat transfer coefficient:
		u = Nu*k/(D_tube_in)
		return u

	def h_conv_tube_V(self, V, T, T_w, D_tube_in):
		rho = self.rho(T)
		mu = self.mu(T)
		Cp = self.Cp(T)
		k = self.k(T)
		# Reynolds:
		Re = D_tube_in*rho*V/mu
		# Prandtl:
		Pr = Cp*mu/k
		if (Pr<=0.1).all():
			if N.logical_and((Re>=1e4), (Re<=1e6)).all():
				# Nusselt:
				Nu = 7.+0.025*(Re*Pr)**0.8
			else:
				print 'Re:',Re
				print 'V:',V, '[m.s-1]'
				print 'Liquid sodium reynolds number outside of Lyon-Martinelli correlation range (1e4 <= Re <= 1e6)'
				return N.nan
				
		else:
			print Pr
			print 'Liquid sodium Prandtl number outside of Lyon-Martinelli correlation range(Pr <= 0.1)'
			return N.nan
		# Heat transfer coefficient:
		u = Nu*k/(D_tube_in)
		return u

	def f_D(self, mf_HC_tube, T, D_tube_in, tube_roughness=45e-6): # Darcy Weissbach
		rho = self.rho(T)
		mu = self.mu(T)
		Cp = self.Cp(T)
		k = self.k(T)

		Pr = Cp*mu/k

		V = mf_HC_tube/(N.pi*(D_tube_in/2.)**2.*rho)
		Re = D_tube_in*rho*V/mu
		#f_D = (1.8*N.log10(Re)-1.5)**-2.
		S = N.log(Re/(1.816*N.log(1.1*Re/(N.log(1.+1.1*Re)))))
		f_D = (-2.*N.log10(tube_roughness/(3.71*D_tube_in)+2.18*S/Re))**(-2.) # Brkic using Lambert W-function approximation to solve Colebrook's implicit equation.
		return f_D

	def p_drop(self, mf_HC_tube, T, D_tube_in, tube_lengths, tube_roughness=45e-6):
		f_D = self.f_D(mf_HC_tube, T, D_tube_in, tube_roughness)
		rho = self.rho(T)
		V = mf_HC_tube/(N.pi*(D_tube_in/2.)**2.*rho)
		Dp = f_D*tube_lengths/D_tube_in*rho*V**2./2.
		return Dp

class Solar_salt(HC):
	'''	
	H. Benoit et al. Renewable and Sustainable Energy Reviews 55 (2016) 298-315
	'''
	def __init__(self, force_T_limits=False):
		print 'Solar salt HC'
		Tmin = 533.
		Tmax = 873.
		HC.__init__(self, Tmin, Tmax, force_T_limits)

	@HC.check_valid
	def rho(self, T):
		TC = T-273.15
		return 2090. - 0.636*TC

	@HC.check_valid
	def mu(self, T):
		TC = T-273.15
		return 2.2714e-2-1.2e-4*TC+2.281e-7*TC**2.-1.474e-10*TC**3.

	@HC.check_valid		
	def k(self, T):
		TC = T-273.15
		return 0.00019*TC + 0.443

	@HC.check_valid
	def Cp(self, T):
		TC = T-273.15
		return 1443. + 0.172*TC

	@HC.check_valid
	def h(self, T):
		TC = T-273.15
		return 1443.*TC+0.172/2.*TC**2.		

	@HC.check_valid
	def s(self, T):
		return (1443.-0.172*273.15)*N.log(T) + 0.172*T

	def V_tube(self, mf_HC_tube, T, D_tube_in):
		rho = self.rho(T)
		# Velocity in the tubes:
		V = mf_HC_tube/(N.pi*(D_tube_in/2.)**2.*rho)
		return V
		
	def h_conv_tube(self, mf_HC_tube, T, T_w, D_tube_in):
		rho = self.rho(T)
		mu = self.mu(T)
		Cp = self.Cp(T)
		k = self.k(T)
		# Velocity in the tubes:
		V = mf_HC_tube/(N.pi*(D_tube_in/2.)**2.*rho)
		# Reynolds:
		Re = D_tube_in*rho*V/mu
		# Prandtl:
		Pr = Cp*mu/k
		Nu = N.zeros(len(T))

		Wu_Pr = N.logical_and(Pr>=1.6,Pr<=23.9)
		Wu_trans_Re = N.logical_and(Re>=2300.,Re<=1e4)
		Wu_trans = N.logical_and(Wu_trans_Re,Wu_Pr)

		Wu_turb_Re = N.logical_and(Re>=1e4,Re<=4.3e4)
		Wu_turb = N.logical_and(Wu_turb_Re,Wu_Pr)

		Gnielinski_Pr = N.logical_and(Pr>=0.5,Pr<=2000.)
		Gnielinski_Re = N.logical_and(Re>=4e3,Re<=5e6)
		Gnielinski = N.logical_and(Gnielinski_Re,Gnielinski_Pr)

		if self.force_Re_limits==False:
			valid = N.logical_or(N.logical_or(Wu_trans,Wu_turb),Gnielinski)
			if ~valid.all():
				print 'Solar salt Prandtl or Reynolds number outside of Wu et al. correlation range'
				return N.nan

		Nu[Wu_trans] = 0.00154*Re[Wu_trans]**1.1*Pr[Wu_trans]**(1./3.) # Wu et al.

		Nu[Wu_turb] = 0.02948*Re[Wu_turb]**0.787*Pr[Wu_turb]**(1./3.) # Wu et al.

		f_D = self.f_D(mf_HC_tube, T, D_tube_in)
		Pr_w = self.Cp(T_w)*self.mu(T_w)/self.k(T_w)
		K = (Pr/Pr_w)**0.11
		Nu_3 = ((f_D/8.)*(Re-1000.)*Pr/(1.+12.7*N.sqrt(f_D/8.)*(Pr**(2./3.)-1.)))*K # Gnielinski
		Nu[Gnielinski] = Nu_3[Gnielinski]

		# Heat transfer coefficient:
		u = Nu*k/(D_tube_in)
		return u

	def f_D(self, mf_HC_tube, T, D_tube_in, tube_roughness=45e-6):
		rho = self.rho(T)
		mu = self.mu(T)
		Cp = self.Cp(T)
		k = self.k(T)

		Pr = Cp*mu/k

		V = mf_HC_tube/(N.pi*(D_tube_in/2.)**2.*rho)
		Re = D_tube_in*rho*V/mu
		#f_D = (1.8*N.log10(Re)-1.5)**-2.
		
		S = N.log(Re/(1.816*N.log(1.1*Re/(N.log(1.+1.1*Re)))))
		f_D = (-2.*N.log10(tube_roughness/(3.71*D_tube_in)+2.18*S/Re))**(-2.) # Brkic using Lambert W-function approximation to solve Colebrook's implicit equation.
		return f_D

	def p_drop(self, mf_HC_tube, T, D_tube_in, tube_lengths):
		f_D = self.f_D(mf_HC_tube, T, D_tube_in)
		rho = self.rho(T)
		V = mf_HC_tube/(N.pi*(D_tube_in/2.)**2.*rho)
		Dp = f_D*tube_lengths/D_tube_in*rho*V**2./2.
		return Dp
	

