def Funcion_Qp_1(x1):  
  import numpy as np
  import scipy.io
  xof=scipy.io.loadmat('Fcn_Qp.mat')['xxoffset']
  xgain=scipy.io.loadmat('Fcn_Qp.mat')['xgain']
  xmin=scipy.io.loadmat('Fcn_Qp.mat')['xymin']
  W1=scipy.io.loadmat('Fcn_Qp.mat')['IW1_1']
  W2=scipy.io.loadmat('Fcn_Qp.mat')['LW2_1']
  W3=scipy.io.loadmat('Fcn_Qp.mat')['LW3_2']
  b1=scipy.io.loadmat('Fcn_Qp.mat')['b1']
  b2=scipy.io.loadmat('Fcn_Qp.mat')['b2']
  b3=scipy.io.loadmat('Fcn_Qp.mat')['b3']
  yof=scipy.io.loadmat('Fcn_Qp.mat')['yxoffset']
  ygain=scipy.io.loadmat('Fcn_Qp.mat')['ygain']
  ymin=scipy.io.loadmat('Fcn_Qp.mat')['ymin']
  
  
  
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
  a3=b3 + W3@a2;
  
  # Output 1
  y1 = mapminmaxreverse(a3);
  return y1
