def Funcion_Qr_1(x1):	
	import numpy as np
	import scipy.io
	xof=scipy.io.loadmat('Fcn_Qr.mat')['xxoffset']
	xgain=scipy.io.loadmat('Fcn_Qr.mat')['xgain']
	xmin=scipy.io.loadmat('Fcn_Qr.mat')['xymin']
	W1=scipy.io.loadmat('Fcn_Qr.mat')['IW1_1']
	W2=scipy.io.loadmat('Fcn_Qr.mat')['LW2_1']
	W3=scipy.io.loadmat('Fcn_Qr.mat')['LW3_2']
	W4=scipy.io.loadmat('Fcn_Qr.mat')['LW4_3']
	b1=scipy.io.loadmat('Fcn_Qr.mat')['b1']
	b2=scipy.io.loadmat('Fcn_Qr.mat')['b2']
	b3=scipy.io.loadmat('Fcn_Qr.mat')['b3']
	b4=scipy.io.loadmat('Fcn_Qr.mat')['b4']
	yof=scipy.io.loadmat('Fcn_Qr.mat')['yxoffset']
	ygain=scipy.io.loadmat('Fcn_Qr.mat')['ygain']
	ymin=scipy.io.loadmat('Fcn_Qr.mat')['ymin']
		
	def mapminmaxapply(x):
		y = x-xof
		y = y*xgain
		y = y+xmin
		return y
	
	# Sigmoid Symmetric Transfer Function
	def tansigapply(n):
		a = 2 / (1 + np.exp(-2*n)) - 1;
		return a
	
	# Map Minimum and Maximum Output Reverse-Processing Function
	def mapminmaxreverse(y):
		x = y-ymin;
		x = x/ygain;
		x = x+yof;
		return x
	  # Input 1
	xp1 = mapminmaxapply(x1);
	
	# Layer 1
	a1 = tansigapply(b1 + W1@xp1);
	
	# Layer 2
	a2 = tansigapply(b2 + W2@a1);
	
	# Layer 3
	a3 = tansigapply(b3 + W3@a2);
	
	# Layer 4
	a4 = b4 + W4@a3;
	
	# Output 1
	y1 = mapminmaxreverse(a4);
	return y1	
