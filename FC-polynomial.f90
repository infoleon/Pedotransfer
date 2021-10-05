real*8 function FnFC (VGa, VGl, VGn, Ks, qfc, zfc)
    
    !Returns THETA-FC according to https://doi.org/10.1016/j.geoderma.2021.115308 (Inforsato & De Jong van Lier, 2021)
    !VGa (m-1) is Van Genuchten - Mualem alpha parameter
    !VGl is Van Genuchten - Mualem l parameter
    !VGn is Van Genuchten - Mualem n parameter
    !Ks (m/d) is saturated hydraulic conductivity
    !qfc (mm/d) is FC flux criterion
    !zfc (m) is depth for FC evaluation
    
    !returns FnFC (m3/m3) = field capacity effective saturation 
    
    !Quirijn de Jong van Lier, October 2021
	
    implicit none
	real*8 VGa, VGl, Vgn, Ks, qfc, zfc
	real*8 a, n, l, k, z, q, a2, n2, l2, k2, z2, q2
    real*8 Prm(0:55)

    a = log10(VGa/100.)
    n = 1./VGn
    l = VGl
    k = log10(Ks*100.)
    z = log10(zfc*100.)
    q = log10(qfc)
    a2 = a**2
    n2 = n**2
    l2 = l**2
    k2 = k**2
    z2 = z**2
    q2 = q**2
    
    
    Data Prm /0.86547, -5.4347e-2, 2.0148e-2, 0, 0.26228, 1.2148e-2, 3.3358e-2, -9.4132e-2, -2.1040e-2, 6.1327e-2, -0.13323, -2.9878e-2, 2.2671e-2, 8.0371e-2, -4.6663e-2, 8.7628e-3, -0.24900, -0.19466, -0.10400, 8.4741e-3, -4.2722e-2, -0.15654, 0.21908, -0.14345, -3.0160e-2, 2.7437e-2, 0.12831, -3.8070e-2, 1.6280, 0.70434, 9.0383e-2, -0.54178, 0.59028, -1.4531, -0.29144, 3.0435e-2, -7.5374e-2, 2.8861e-2, 2.1007e-2, 3.8902e-3, -6.2775e-2, 5.2356e-2, -1.8087e-2, 2.6784e-2, 9.6550e-2, -3.1786e-2, 3.7102e-2, -7.4344e-3, 0.35107, 0, -0.62512, 0.19079/
    
    FnFC = Prm(0) + a*Prm(1) + a*l*Prm(2) + a*z*Prm(3) + a*n2*Prm(4) + a*k2*Prm(5) + a*z2*Prm(6) + a*n*k*Prm(7) + a*n*l*Prm(8) + a*n*q*Prm(9) + a*n*z*Prm(10) + a*k*q*Prm(11) + a*k*z*Prm(12) + a2*Prm(13) + a2*n*Prm(14) + a*a2*Prm(15)
    FnFC = FnFC   + n*k*Prm(16) + n*l*Prm(17) + n*k2*Prm(18) + n*l2*Prm(19) + n*q2*Prm(20) + n*z2*Prm(21) + n*k*q*Prm(22) + n*k*z*Prm(23) + n*k*l*Prm(24) + n*q*l*Prm(25) + n*q*z*Prm(26) + n*z*l*Prm(27) + n2*Prm(28) + n2*k*Prm(29) + n2*l*Prm(30) + n2*q*Prm(31) + n2*z*Prm(32) + n*n2*Prm(33)
    FnFC = FnFC   + k*Prm(34) + k*l*Prm(35) + k*q*Prm(36) + k*q2*Prm(37) + k*z2*Prm(38) + k*q*l*Prm(39) + k*q*z*Prm(40) + k2*Prm(41) + k2*q*Prm(42) + k2*z*Prm(43)
    FnFC = FnFC   + l*Prm(44) + l*q*Prm(45) + l*z*Prm(46) + l2*Prm(47) + q*Prm(48) + q2*z*Prm(49) + z*Prm(50) + z2*Prm(51)
    
    
    return
	end function