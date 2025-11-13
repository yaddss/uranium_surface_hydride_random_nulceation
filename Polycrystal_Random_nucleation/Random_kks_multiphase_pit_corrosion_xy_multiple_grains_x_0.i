#  length_scale = 5 nm, time_scale = 1 s, energy_scale = 1.25e-16 J;
#  1 m = 1e+3 * 1 mm = 2e+8 * 5e-9 m = 2e+8 * 5 nm; 1 s = 1 * 1 s; 1 J = 8e+15 * 1.25e-16 J
#  RT= 8.3145*(240+273.15) J/mol
#  KKS multiphase simple example, constant mobility, periodic BCs
#  Notes: mobility can't be customized by setting Laa,Lab,Lac,.....
#  T = 240 +273 K

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 600
  ny = 600
  nz = 0
  xmin = -300
  xmax = 300
  ymin = -300
  ymax = 300
  zmin = 0
  zmax = 0
  elem_type = QUAD4  #3D HEX8, 2D QUAD4
[]

[GlobalParams]
  # displacement
  displacements = 'disp_x disp_y'
  # phase-field potential barrier
  wi = 48 # Height of potential barrier
[]


[AuxVariables]
  [./phase]
    order = FIRST
    family = MONOMIAL
  [../]
  # Chemical free Energy
  [./f_ch]
    order = FIRST
    family = MONOMIAL
  [../]
  # Grain boundary
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  # driving force
  [./chdf_12]
    order = FIRST
    family = MONOMIAL
  [../]
  [./eldf_eta2]
    order = FIRST
    family = MONOMIAL
  [../]
  [./tot_delta_eta2]
    order = FIRST
    family = MONOMIAL
  [../]
  # stress
  [./stress_11]
    order = FIRST
    family = MONOMIAL
  [../]
  [./stress_12]
    order = FIRST
    family = MONOMIAL
  [../]
  [./stress_22]
    order = FIRST
    family = MONOMIAL
  [../]
  [./vonmises]
    order = FIRST
    family = MONOMIAL
  [../]
[]
# []

[AuxKernels]
  # AuxKernel that calculates the elastic strain energy
  [./phase]
    type = MaterialRealAux   #The computed value will be the volume-averaged quantity over the element.
    variable = phase
    property = phase_sum
    execute_on = 'initial timestep_end'
  [../]
  # AuxKernel that calculates the chemical free Energy
  [./f_ch]
    type = MaterialRealAux   #The computed value will be the volume-averaged quantity over the element.
    variable = f_ch
    property =   f_ch_mat
    execute_on = 'initial timestep_end'
  [../]
  # AuxKernel that calculates the driving force
  [./chemical_driving_force_12]
    type = MaterialRealAux   #The computed value will be the volume-averaged quantity over the element.
    variable = chdf_12
    property =   deltaG_ch_12
    execute_on = 'initial timestep_end'
  [../]
  [./elastic_driving_force_eta2]
    type = MaterialRealAux   #The computed value will be the volume-averaged quantity over the element.
    variable = eldf_eta2
    property =   deltaG_el_eta2
    execute_on = 'initial timestep_end'
  [../]
  [./tot_driving_force_eta2]
    type = MaterialRealAux   #The computed value will be the volume-averaged quantity over the element.
    variable = tot_delta_eta2
    property =   deltaG_eta2
    execute_on = 'initial timestep_end'
  [../]
  # AuxKernel that calculates the GB term
  [bnds]
    type = BndsCalcAux       
    variable = bnds
    execute_on = 'initial'
	v = 'gr0 gr1 gr2 gr3 gr4 gr5'
  [../]
  [./sigma_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_11
    index_i = 0
    index_j = 0
	execute_on = 'initial timestep_end'
  [../]
  [./sigma_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_12
    index_i = 0
    index_j = 1
	execute_on = 'initial timestep_end'
  [../]
  [./sigma_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_22
    index_i = 1
    index_j = 1
	execute_on = 'initial timestep_end'
  [../]
  [./vonmises]
    type = RankTwoScalarAux
    variable = vonmises
    rank_two_tensor = stress
    scalar_type = VonMisesStress
    execute_on = 'initial timestep_end'
  [../]
[]

[Variables]
  # grains 
  [./gr0]
    order = FIRST
    family = LAGRANGE
  [../]
  [./gr1]
    order = FIRST
    family = LAGRANGE
  [../]
  [./gr2]
    order = FIRST
    family = LAGRANGE
  [../]
  [./gr3]
    order = FIRST
    family = LAGRANGE
  [../]
  [./gr4]
    order = FIRST
    family = LAGRANGE
  [../]
  [./gr5]
    order = FIRST
    family = LAGRANGE
  [../]
  
  # initial stress
  [./in_stress_11]
    order = FIRST
    family = LAGRANGE
  [../]
  [./in_stress_12]
    order = FIRST
    family = LAGRANGE
  [../]
  [./in_stress_22]
    order = FIRST
    family = LAGRANGE
  [../]
  
  
  # order parameter of metal matrix
  [./eta1]
    order = FIRST
    family = LAGRANGE
  [../]
  # order parameter of hydride precipitate
  [./eta2]
    order = FIRST
    family = LAGRANGE
  [../]
  # order parameter of carbide inclusion
  [./eta3]
    order = FIRST
    family = LAGRANGE
  [../]
  
  # solute hydrogen concentration
  [./c]
    order = FIRST
    family = LAGRANGE
  [../]

  # matrix phase solute concentration
  [./c1]
    order = FIRST
    family = LAGRANGE
  [../]
  # precipitate phase solute concentration
  [./c2]
    order = FIRST
    family = LAGRANGE
  [../]
  # carbide phase solute concentration
  [./c3]
    order = FIRST
    family = LAGRANGE
  [../]
  # displacement 
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
   
  # Lagrange multiplier
  [./lambda]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  [../]
[]

[UserObjects]
  [./inserter]
    type = DiscreteNucleationInserter
    hold_time = 1
    probability =P  #0.00005
    radius = 7
  [../]
  [./map]
    type = DiscreteNucleationMap
    int_width = 6
    periodic = eta2
    inserter = inserter
  [../]
  [exo_variables]
    type = SolutionUserObject
    mesh = 'poly_grain_growth_2D_eldrforce_x_-0.1_out.e'
    system_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 in_stress_11 in_stress_12 in_stress_22'
    timestep = LATEST
  []
[]

[Functions]
  # # grains gr0 gr1 gr2
  # [./ic_fun_gr0]
    # type = ParsedFunction
    # value = '0.5*(1.0+tanh(-2.0*alpha/gb_width*(x-x1)))+0.5*(1.0-tanh(-2.0*alpha/gb_width*(x-x2)))'
    # vars = 'alpha gb_width  x1 x2'
    # vals = '2.0    24.0  -80.0  80.0'
  # [../]
  # [./ic_fun_gr1_gr2]
    # type = ParsedFunction
    # value = '1-(0.5*(1.0+tanh(-2.0*alpha/gb_width*(x-x1)))+0.5*(1.0-tanh(-2.0*alpha/gb_width*(x-x2))))'
    # vars = 'alpha gb_width  x1 x2'
    # vals = '2.0    24.0  -80.0  80.0'
  # [../]
  # [./ic_fun_gr1]
    # type = ParsedFunction
    # value = '0.5*(1.0-tanh(-2.0*alpha/gb_width*(x-x1)))-0.5*(1.0-tanh(-2.0*alpha/gb_width*(x-x3)))'
    # vars = 'alpha gb_width  x1 x3'
    # vals = '2.0    24.0  -80.0  0'
  # [../]
  # [./ic_fun_gr2]
    # type = ParsedFunction
    # value = '(1-(0.5*(1.0+tanh(-2.0*alpha/gb_width*(x-x1)))+0.5*(1.0-tanh(-2.0*alpha/gb_width*(x-x2)))))-(0.5*(1.0-tanh(-2.0*alpha/gb_width*(x-x1)))-0.5*(1.0-tanh(-2.0*alpha/gb_width*(x-x3))))'
    # vars = 'alpha gb_width  x1 x2 x3'
    # vals = '2.0    24.0  -80.0 80.0 0'
  # [../]
  
  # # test single nucleation
 [./ic_fun_gr0]
    type = SolutionFunction
    solution = exo_variables
    from_variable = gr0
  [../]
  [./ic_fun_gr1]
    type = SolutionFunction
    solution = exo_variables
    from_variable = gr1
  [../]
  [./ic_fun_gr2]
    type = SolutionFunction
    solution = exo_variables
    from_variable = gr2
  [../]
  [./ic_fun_gr3]
    type = SolutionFunction
    solution = exo_variables
    from_variable = gr3
  [../]
  [./ic_fun_gr4]
    type = SolutionFunction
    solution = exo_variables
    from_variable = gr4
  [../]
  [./ic_fun_gr5]
    type = SolutionFunction
    solution = exo_variables
    from_variable = gr5
  [../]
  [./ic_fun_in_stress_11]
    type = SolutionFunction
    solution = exo_variables
    from_variable = in_stress_11
  [../]
  [./ic_fun_in_stress_12]
    type = SolutionFunction
    solution = exo_variables
    from_variable = in_stress_12
  [../]
  [./ic_fun_in_stress_22]
    type = SolutionFunction
    solution = exo_variables
    from_variable = in_stress_22
  [../]
[]

[ICs]
  # [./eta1]
    # variable = eta1
    # type = FunctionIC
    # function = ic_fun_eta1
  # [../]
  # [./eta2]
    # variable = eta2
    # type = FunctionIC
    # function = ic_fun_eta2
  # [../]
 [./eta1]
    type = ConstantIC
    variable = eta1
    value = 1
  [../]
  [./eta2]
    type = ConstantIC
    variable = eta2
    value = 0
  [../]
  [./eta3]
    type = ConstantIC
    variable = eta3
    value = 0
  [../]
  # [./c]
    # variable = c
    # type = FunctionIC
    # function = ic_fun_c
  # [../]
  [./c]
    type = ConstantIC
    variable = c
	value = 0.05
  [../]
  # grains gr0 gr1 gr2 gr3 gr4 gr5
  [./gr0]
    variable = gr0
    type = FunctionIC
    function = ic_fun_gr0
  [../]
  [./gr1]
    variable = gr1
    type = FunctionIC
    function = ic_fun_gr1
  [../]
  [./gr2]
    variable = gr2
    type = FunctionIC
    function = ic_fun_gr2
  [../]
  [./gr3]
    variable = gr3
    type = FunctionIC
    function = ic_fun_gr3
  [../]
  [./gr4]
    variable = gr4
    type = FunctionIC
    function = ic_fun_gr4
  [../]
  [./gr5]
    variable = gr5
    type = FunctionIC
    function = ic_fun_gr5
  [../]
  
  [./in_stress_11]
    variable = in_stress_11
    type = FunctionIC
    function = ic_fun_in_stress_11
  [../]
  [./in_stress_12]
    variable = in_stress_12
    type = FunctionIC
    function = ic_fun_in_stress_12
  [../]
  [./in_stress_22]
    variable = in_stress_22
    type = FunctionIC
    function = ic_fun_in_stress_22
  [../]
[]

[BCs]  #default noflux
  # # periodic BC
  # [./Periodic]
    # [./all_bcs]
      # auto_direction = 'x'
    # [../]
  # [../]
[]
 

[Materials]
  [./probability]
    # This is a made up toy nucleation rate it should be replaced by
    # classical nucleation theory in a real simulation.
    type = ParsedMaterial
    f_name = P
    args = 'c'
    function = '0.5*exp(-16*3.14*(0.5)^3/3/(deltaG_eta2*1e9)^2/kbT)'
	material_property_names = 'deltaG_eta2 kbT'
    outputs = exodus
  [../]
  # nucleation feta
  [./nucleationfpn]
    type = DiscreteNucleation
	op_names  = eta2
    op_values = 1
    map = map
	penalty = 200
	penalty_mode = MATCH
    f_name = feta2
    outputs = exodus
  [../]	
  # Test the sum of eta is equal to one
  [./etasummat]
    type = ParsedMaterial
    f_name = etasum
    args = 'eta1 eta2 eta3'
    material_property_names = 'h1 h2 h3'
    function = 'h1+h2+h3-1'
    outputs = exodus
	# use_displaced_mesh = true
  [../]
  # This parsed material creates a single property for visualization purposes.
  # It will be 0 for phase 1, -1 for phase 2, and 1 for phase 3
  [./phasemap]
    type = ParsedMaterial
    f_name = phase_sum
    args = 'eta2 eta3 bnds'
    function = 'if(eta3>0.5,1,0)-if(eta2>0.5,1,0)'
  [../]
 # Chemical free energy of the matrix
  [./f1]
    type = DerivativeParsedMaterial
    f_name = F1
    args = 'c1'
    function = '(235.0916*c1^2-14.6697*c1+0.2288)' #(233.9625*c1^2-14.1313*c1+0.2134)
    outputs = exodus
  [../]
  # Chemical Free energy of the precipitate phase
  [./f2]
    type = DerivativeParsedMaterial
    f_name = F2
    args = 'c2'
    function = '(368.7712*c2^2-556.8445*c2+207.4808)' #(185.3169*c2^2-281.6816*c2+104.3110)
    outputs = exodus
  [../]
  # Chemical Free energy of the carbide phase
  [./f3]
    type = DerivativeParsedMaterial
    f_name = F3
    args = 'c3'
    function = '(338.5319*c3^2-34.6657*c3+1.0987)'
    outputs = exodus
  [../]
  # solute-grain boundary interaction energy
  [./f_sg_mat]
    type = DerivativeParsedMaterial
    f_name = f_sg_mat
    args = 'gr0 gr1 gr2 gr3 gr4 gr5 c'
    function = '-p0*c*((-0.5*gr0^2+0.25*gr0^4)+(-0.5*gr1^2+0.25*gr1^4)+(-0.5*gr2^2+0.25*gr2^4)+(-0.5*gr3^2+0.25*gr3^4)+(-0.5*gr4^2+0.25*gr4^4)+(-0.5*gr5^2+0.25*gr5^4)+gamma*(gr0^2*gr1^2+gr0^2*gr2^2+gr0^2*gr3^2+gr0^2*gr4^2+gr0^2*gr5^2+gr1^2*gr2^2+gr1^2*gr3^2+gr1^2*gr4^2+gr1^2*gr5^2+gr2^2*gr3^2+gr2^2*gr4^2+gr2^2*gr5^2+gr3^2*gr4^2+gr3^2*gr5^2+gr4^2*gr5^2)+0.25)'
    material_property_names = 'gamma p0'
  [../]
  [./gamma]
    type = DerivativeParsedMaterial
    f_name = gamma
    args = 'gr0 gr1 gr2 gr3 gr4 gr5'
    function = '4*zeta/theta_ref*(1-plog(zeta/theta_ref,0.005))'
    material_property_names = 'zeta theta_ref'
  [../]
  [./zeta]
    type = DerivativeParsedMaterial
    f_name = zeta
    args = 'gr0 gr1 gr2 gr3 gr4 gr5' 
    function = '(gr0^2*gr1^2*theta01+gr0^2*gr2^2*theta02+gr0^2*gr3^2*theta03+gr0^2*gr4^2*theta04+gr0^2*gr5^2*theta05+gr1^2*gr2^2*theta12+gr1^2*gr3^2*theta13+gr1^2*gr4^2*theta14+gr1^2*gr5^2*theta15+gr2^2*gr3^2*theta23+gr2^2*gr4^2*theta24+gr2^2*gr5^2*theta25+gr3^2*gr4^2*theta34+gr3^2*gr5^2*theta35+gr4^2*gr5^2*theta45)/(gr0^2*gr1^2+gr0^2*gr2^2+gr0^2*gr3^2+gr0^2*gr4^2+gr0^2*gr5^2+gr1^2*gr2^2+gr1^2*gr3^2+gr1^2*gr4^2+gr1^2*gr5^2+gr2^2*gr3^2+gr2^2*gr4^2+gr2^2*gr5^2+gr3^2*gr4^2+gr3^2*gr5^2+gr4^2*gr5^2)'
    material_property_names = 'theta01 theta02 theta03 theta04 theta05 theta12 theta13 theta14 theta15 theta23 theta24 theta25 theta34 theta35 theta45'
  [../]
  # constant properties
  [./constants_gr]
    type = GenericFunctionMaterial  
    prop_names  = 'theta01 theta02 theta03 theta04 theta05 theta12 theta13 theta14 theta15 theta23 theta24 theta25 theta34 theta35 theta45 theta_ref p0'   # p0 = 1*omega
    prop_values = '3       6       9       12      15      3       6       9       12      3       6       9       3       6       3
				   15
				   1*24*60/3e-8*8e+15/(2e+8)^3*1'  
  [../]
  
  
  
  # Switching functions for each phase
  # It's the fifth order polynomial,
  # hi = etai^5+5etai^4(1-etai)+10etai^3(1-etai)^2+15/4[etai^2(1-etai)^3-etai^2(1-etai)(etaj-etak)^2]
  # sum<i=1,N=3>hi = 1
  # h1(eta1, eta2, eta3)
  [./h1]
    type = SwitchingFunction3PhaseMaterial
    eta_i = eta1
    eta_j = eta2
    eta_k = eta3
    f_name = h1
  [../]
  # h2(eta1, eta2, eta3)
  # defualt function
  [./h2]
    type = SwitchingFunction3PhaseMaterial
    eta_i = eta2
    eta_j = eta3
    eta_k = eta1
    f_name = h2
  [../]
  # h3(eta1, eta2, eta3)
  [./h3]
    type = SwitchingFunction3PhaseMaterial
    eta_i = eta3
    eta_j = eta1
    eta_k = eta2
    f_name = h3
  [../]

  # Potential barrier functions for each phase
  # defualt function wi*etai^2(1-etai)^2
  # sum_g = (w1+w2)*eta1^2*eta2^2+(w1+w3)*eta1^2*eta3^2+(w2+w3)*eta2^2*eta3^2
  #       +  2*(w1*eta1+w2*eta2+w3*eta3)*eta1*eta2*eta3
  # if w1=w2=w3=w, then
  # sum_g = 2w(eta1^2*eta2^2+eta1^2*eta3^2+eta2^2*eta3^2+eta1*eta2*eta3*(eta1+eta2+eta3))
  #       = 2w*(sum<i\neq j,i=1,N=3>(etai*etaj)+eta1*eta2*eta3)
  [./g1]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta1
    function_name = g1
  [../]
  [./g2]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta2
    function_name = geta2
  [../]
  [./g3]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta3
    function_name = g3
  [../]
  
  
  # g(eta)
  [./g_eta_nucleation]
  # add the chemical and nucleation free energy contributions together
    type = DerivativeSumMaterial
    derivative_order = 2
    args = 'eta2'
    sum_materials = 'geta2 feta2'
	f_name = g2
  [../] 
  
  
  # constant properties
  [./constants1]
    #diffusion coefficient, evolution coefficient, gradient energy coefficient
	# D = (m^2/s) * length_scale^2/time_scale
	# M = D*Vm/R/T(m^5/s/J) * length_scale^5/time_scale/energy_scale
	# L = (m^3/s/J) * length_scale^3/time_scale/energy_scale
	# kappa = 3/4*sigma*lambda(J/m) * energy_scale/length_scale
	# omega = 24*sigma/lambda(J/m^3) * energy_scale/length_scale^3 
	# rho = kg/m^3 = J*s^2/m^5 *energy_scale*time_scale/length_scale^5 	
	# D1=D2=0.019exp(-11.6[kcal/mol]/T/R)
	# D3=0.037exp(-14.3[kcal/mol]/R/T)
	# rho_metal = 19.07 g/cm^3
	# rho_precipitate = 10.92 g/cm^3
    type = GenericFunctionMaterial  
    prop_names  = 'Ref_D1 Ref_D2 Ref_D3 Ref_D_GB  L  kappa  omega  mu  kbT'   # mu = 3*omega
    prop_values = '2.1727e-11*(2e+08)^2/1 
	               2.1727e-11*(2e+08)^2/1 
				   2.1727e-11*(2e+08)^2/1 
				   2.1727e-11*100*(2e+08)^2/1 
	               2.1727e-10*(2e+08)^3/1/8e+15    
				   3/4*60*3e-8*8e+15/2e+08
				   24*60/3e-8*8e+15/(2e+8)^3
				   3*24*60/3e-8*8e+15/(2e+8)^3*1
				   1.380649e-23*(240+273)'  
  [../]

  # constant properties
  [./constants2]
    type = GenericFunctionMaterial  
    prop_names  = 'gr0_eih11  gr0_eih22  gr0_eih12   gr1_eih11 gr1_eih22 gr1_eih12   gr2_eih11 gr2_eih22 gr2_eih12   gr3_eih11 gr3_eih22 gr3_eih12   gr4_eih11 gr4_eih22 gr4_eih12   gr5_eih11 gr5_eih22 gr5_eih12'   
    prop_values = '0.1259     0.1259     0           0.1259    0.1259    0           0.1259    0.1259    0           0.1259    0.1259    0           0.1259    0.1259    0           0.1259    0.1259    0'  
  [../]
  

  # d^2F/dc^2
  [./second_Derivative_chemical_free_energy]
    type = ParsedMaterial
    f_name = fcc
    args = 'c1 c2 c3 eta1 eta2 eta3'
    function = 'd2F1*d2F2*d2F3/(h1*d2F2*d2F3+h2*d2F1*d2F3+h3*d2F1*d2F2)'
    material_property_names = 'd2F1:=D[F1,c1,c1] d2F2:=D[F2,c2,c2] d2F3:=D[F3,c3,c3] h1 h2 h3'
  [../]
  # Mobility M
  [./Mobility_M]
    type = ParsedMaterial
    f_name = M
    function = 'D/fcc'
    material_property_names = 'D fcc'
    outputs = exodus
  [../]
  # fsg,c,gr0*D/fcc
  [./Mobility_Msgr0]
    type = ParsedMaterial
    f_name = M_sgr0
    function = 'd2f_sg_matdcdgr0*M'
    args = 'c c1 c2 c3 eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    material_property_names = 'd2f_sg_matdcdgr0:=D[f_sg_mat,c,gr0] M'
    outputs = exodus
  [../]
  # fsg,c,gr1*D/fcc
  [./Mobility_Msgr1]
    type = ParsedMaterial
    f_name = M_sgr1
    function = 'd2f_sg_matdcdgr1*M'
    args = 'c c1 c2 c3 eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    material_property_names = 'd2f_sg_matdcdgr1:=D[f_sg_mat,c,gr1] M'
    outputs = exodus
  [../]
  # fsg,c,gr2*D/fcc
  [./Mobility_Msgr2]
    type = ParsedMaterial
    f_name = M_sgr2
    function = 'd2f_sg_matdcdgr2*M'
    args = 'c c1 c2 c3 eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    material_property_names = 'd2f_sg_matdcdgr2:=D[f_sg_mat,c,gr2] M'
    outputs = exodus
  [../]
   # fsg,c,gr3*D/fcc
  [./Mobility_Msgr3]
    type = ParsedMaterial
    f_name = M_sgr3
    function = 'd2f_sg_matdcdgr3*M'
    args = 'c c1 c2 c3 eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    material_property_names = 'd2f_sg_matdcdgr3:=D[f_sg_mat,c,gr3] M'
    outputs = exodus
  [../]
  # fsg,c,gr1*D/fcc
  [./Mobility_Msgr4]
    type = ParsedMaterial
    f_name = M_sgr4
    function = 'd2f_sg_matdcdgr4*M'
    args = 'c c1 c2 c3 eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    material_property_names = 'd2f_sg_matdcdgr4:=D[f_sg_mat,c,gr4] M'
    outputs = exodus
  [../]
  # fsg,c,gr2*D/fcc
  [./Mobility_Msgr5]
    type = ParsedMaterial
    f_name = M_sgr5
    function = 'd2f_sg_matdcdgr5*M'
    args = 'c c1 c2 c3 eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    material_property_names = 'd2f_sg_matdcdgr5:=D[f_sg_mat,c,gr5] M'
    outputs = exodus
  [../]
 
 
  # !!!! Note h1,h2, and h3 cannot be used to replace the interpolation function, because only the function value can be passed, and the other derivative values cannot be passed.
  # !!!! that is dh1_gr0deta1, dh1_gr0deta2, dh1_gr0deta3 do not exist.
  # h_gr0
  [./h1_gr0]
    type = DerivativeParsedMaterial
    f_name = h1_gr0
    args = 'eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    function = 'eta1^2/4*(15*(1-eta1)*(1+eta1-(eta3-eta2)^2)+eta1*(9*eta1^2-5))*gr0' 
  [../]
  [./h2_gr0]
    type = DerivativeParsedMaterial
    f_name = h2_gr0
    args = 'eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    function = 'eta2^2/4*(15*(1-eta2)*(1+eta2-(eta1-eta3)^2)+eta2*(9*eta2^2-5))*gr0' 
  [../]
  [./h3_gr0]
    type = DerivativeParsedMaterial
    f_name = h3_gr0
    args = 'eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    function = 'eta3^2/4*(15*(1-eta3)*(1+eta3-(eta2-eta1)^2)+eta3*(9*eta3^2-5))*gr0' 
  [../]
  
  # h_gr1
  [./h1_gr1]
    type = DerivativeParsedMaterial
    f_name = h1_gr1
    args = 'eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    function = 'eta1^2/4*(15*(1-eta1)*(1+eta1-(eta3-eta2)^2)+eta1*(9*eta1^2-5))*gr1' 
  [../]
  [./h2_gr1]
    type = DerivativeParsedMaterial
    f_name = h2_gr1
    args = 'eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    function = 'eta2^2/4*(15*(1-eta2)*(1+eta2-(eta1-eta3)^2)+eta2*(9*eta2^2-5))*gr1' 
  [../]
  [./h3_gr1]
    type = DerivativeParsedMaterial
    f_name = h3_gr1
    args = 'eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    function = 'eta3^2/4*(15*(1-eta3)*(1+eta3-(eta2-eta1)^2)+eta3*(9*eta3^2-5))*gr1' 
  [../]
 
   # h_gr2
  [./h1_gr2]
    type = DerivativeParsedMaterial
    f_name = h1_gr2
    args = 'eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    function = 'eta1^2/4*(15*(1-eta1)*(1+eta1-(eta3-eta2)^2)+eta1*(9*eta1^2-5))*gr2' 
  [../]
  [./h2_gr2]
    type = DerivativeParsedMaterial
    f_name = h2_gr2
    args = 'eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    function = 'eta2^2/4*(15*(1-eta2)*(1+eta2-(eta1-eta3)^2)+eta2*(9*eta2^2-5))*gr2' 
  [../]
  [./h3_gr2]
    type = DerivativeParsedMaterial
    f_name = h3_gr2
    args = 'eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    function = 'eta3^2/4*(15*(1-eta3)*(1+eta3-(eta2-eta1)^2)+eta3*(9*eta3^2-5))*gr2' 
  [../]
 
  # h_gr3
  [./h1_gr3]
    type = DerivativeParsedMaterial
    f_name = h1_gr3
    args = 'eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    function = 'eta1^2/4*(15*(1-eta1)*(1+eta1-(eta3-eta2)^2)+eta1*(9*eta1^2-5))*gr3' 
  [../]
  [./h2_gr3]
    type = DerivativeParsedMaterial
    f_name = h2_gr3
    args = 'eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    function = 'eta2^2/4*(15*(1-eta2)*(1+eta2-(eta1-eta3)^2)+eta2*(9*eta2^2-5))*gr3' 
  [../]
  [./h3_gr3]
    type = DerivativeParsedMaterial
    f_name = h3_gr3
    args = 'eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    function = 'eta3^2/4*(15*(1-eta3)*(1+eta3-(eta2-eta1)^2)+eta3*(9*eta3^2-5))*gr3' 
  [../]
 
  # h_gr4
  [./h1_gr4]
    type = DerivativeParsedMaterial
    f_name = h1_gr4
    args = 'eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    function = 'eta1^2/4*(15*(1-eta1)*(1+eta1-(eta3-eta2)^2)+eta1*(9*eta1^2-5))*gr4' 
  [../]
  [./h2_gr4]
    type = DerivativeParsedMaterial
    f_name = h2_gr4
    args = 'eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    function = 'eta2^2/4*(15*(1-eta2)*(1+eta2-(eta1-eta3)^2)+eta2*(9*eta2^2-5))*gr4' 
  [../]
  [./h3_gr4]
    type = DerivativeParsedMaterial
    f_name = h3_gr4
    args = 'eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    function = 'eta3^2/4*(15*(1-eta3)*(1+eta3-(eta2-eta1)^2)+eta3*(9*eta3^2-5))*gr4' 
  [../]
 
  # h_gr5
  [./h1_gr5]
    type = DerivativeParsedMaterial
    f_name = h1_gr5
    args = 'eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    function = 'eta1^2/4*(15*(1-eta1)*(1+eta1-(eta3-eta2)^2)+eta1*(9*eta1^2-5))*gr5' 
  [../]
  [./h2_gr5]
    type = DerivativeParsedMaterial
    f_name = h2_gr5
    args = 'eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    function = 'eta2^2/4*(15*(1-eta2)*(1+eta2-(eta1-eta3)^2)+eta2*(9*eta2^2-5))*gr5' 
  [../]
  [./h3_gr5]
    type = DerivativeParsedMaterial
    f_name = h3_gr5
    args = 'eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    function = 'eta3^2/4*(15*(1-eta3)*(1+eta3-(eta2-eta1)^2)+eta3*(9*eta3^2-5))*gr5' 
  [../]
 
  [./D_GB]
    type = DerivativeParsedMaterial
    f_name = D_GB
    args = 'gr0 gr1 gr2 gr3 gr4 gr5'
    function = 'Ref_D_GB/2*(1+tanh(7.5*(zeta/theta_ref-2/3)))'
    material_property_names = 'Ref_D_GB zeta theta_ref'
  [../]
 
 
  # Interpolated Coefficients for Diffusion 
  [./D]
    type = ParsedMaterial
    material_property_names = 'Ref_D1 Ref_D2 Ref_D3 h1 h2 h3 D_GB'
    function = Ref_D1*h1+Ref_D2*h2+Ref_D3*h3 #+16*(gr0^2*gr1^2+gr0^2*gr2^2+gr0^2*gr3^2+gr0^2*gr4^2+gr0^2*gr5^2+gr1^2*gr2^2+gr1^2*gr3^2+gr1^2*gr4^2+gr1^2*gr5^2+gr2^2*gr3^2+gr2^2*gr4^2+gr2^2*gr5^2+gr3^2*gr4^2+gr3^2*gr5^2+gr4^2*gr5^2)*D_GB
    f_name = D
	args = 'gr0 gr1 gr2 gr3 gr4 gr5 eta1 eta2 eta3'
	outputs = exodus
  [../]
  
  # Interpolation for diffusion equation
  [./Dh1]
    type = DerivativeParsedMaterial
    material_property_names = 'D h1'
    function = D*h1
    f_name = Dh1
  [../]
  [./Dh2]
    type = DerivativeParsedMaterial
    material_property_names = 'D h2'
    function = D*h2
    f_name = Dh2
  [../]
  [./Dh3]
    type = DerivativeParsedMaterial
    material_property_names = 'D h3'
    function = D*h3
    f_name = Dh3
  [../]
  
  # Chemical Free energy of the system
  [./chemical_free_energy_p]
    type = ParsedMaterial
    args = 'eta1 eta2 eta3 c1 c2 c3'
    f_name = f_ch_mat
    function = h1*F1+h2*F2+h3*F3+omega*(g1+g2+g3)
    material_property_names = 'h1 h2 h3 F1 F2 F3 g1 g2 g3 omega'
  [../]
  
  # Elastic energy of the system
  [./elastic_free_energy_p]
    type = ElasticEnergyMaterial
    args = 'eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    f_name = f_el_mat
    derivative_order = 2
  [../]
  
  #Driving force
  [./ch_driving_force_12]
    type = DerivativeParsedMaterial
    args = 'eta1 eta2 eta3 c1 c2 c3'
    f_name = deltaG_ch_12
    function = 'F2-F1-chmu*(c2-c1)'
    material_property_names = 'chmu:=D[F1,c1] F1 F2'
  [../]
  [./el_driving_force_eta2]
    type = DerivativeParsedMaterial
    args = 'eta1 eta2 eta3 c1 c2 c3 stress_11 stress_12 stress_22 gr0 gr1 gr2 gr3 gr4 gr5'
    f_name = deltaG_el_eta2
    function = '-(stress_11)*(gr0*gr0_eih11+gr1*gr1_eih11+gr2*gr2_eih11+gr3*gr3_eih11+gr4*gr4_eih11+gr5*gr5_eih11)-(stress_22)*(gr0*gr0_eih22+gr1*gr1_eih22+gr2*gr2_eih22+gr3*gr3_eih22+gr4*gr4_eih22+gr5*gr5_eih22)-2*(stress_12)*(gr0*gr0_eih12+gr1*gr1_eih12+gr2*gr2_eih12+gr3*gr3_eih12+gr4*gr4_eih12+gr5*gr5_eih12)'
    material_property_names = 'gr0_eih11 gr0_eih22 gr0_eih12 gr1_eih11 gr1_eih22 gr1_eih12 gr2_eih11 gr2_eih22 gr2_eih12 gr3_eih11 gr3_eih22 gr3_eih12 gr4_eih11 gr4_eih22 gr4_eih12 gr5_eih11 gr5_eih22 gr5_eih12'
  [../]
  [./tot_driving_force_eta2]
    type = DerivativeParsedMaterial
    f_name = deltaG_eta2
    function = '-deltaG_ch_12-deltaG_el_eta2'
    material_property_names = 'deltaG_ch_12 deltaG_el_eta2'
  [../] 

  # Mechanical properties gr0
  [./Stiffness_matrix_gr0]
    type = ComputeElasticityTensor
    fill_method = symmetric9     #select symmetry
    base_name = C_matrix_gr0   #C11 C12 C13 C22 C23 C33 C44 C55 C66
    C_ijkl = '207.6 50.6 24.5 185.9	102.7 242.7 108.4 55.7 63.2' #T=250C
    # Rotation Matrix
    euler_angle_1 = 0. 
    euler_angle_2 = 0.
    euler_angle_3 = 0.
  [../]
  [./Stiffness_precipitate_gr0]
    type = ComputeElasticityTensor
    fill_method = symmetric9     #select symmetry
    base_name = C_precipitate_gr0  #C11 C12 C13 C22 C23 C33 C44 C55 C66
    C_ijkl = '222 70 70 222 70 222 58 58 58'
	# Rotation Matrix
    euler_angle_1 = 0. 
    euler_angle_2 = 0.
    euler_angle_3 = 0.
  [../]
  [./Stiffness_carbide_gr0]
    type = ComputeElasticityTensor
    fill_method = symmetric9     #select symmetry
    base_name = C_carbide_gr0
    C_ijkl = '395 121 121 395 121 395 64.1 64.1 64.1'
    # Rotation Matrix
    euler_angle_1 = 0. 
    euler_angle_2 = 0.
    euler_angle_3 = 0.
  [../] 

  # Mechanical properties gr1
  [./Stiffness_matrix_gr1]
    type = ComputeElasticityTensor
    fill_method = symmetric9     #select symmetry
    base_name = C_matrix_gr1   #C11 C12 C13 C22 C23 C33 C44 C55 C66
    C_ijkl = '207.6 50.6 24.5 185.9	102.7 242.7 108.4 55.7 63.2' #T=250C
    # Rotation Matrix
    euler_angle_1 = 3. 
    euler_angle_2 = 0.
    euler_angle_3 = 0.
  [../]
  [./Stiffness_precipitate_gr1]
    type = ComputeElasticityTensor
    fill_method = symmetric9     #select symmetry
    base_name = C_precipitate_gr1  #C11 C12 C13 C22 C23 C33 C44 C55 C66
    C_ijkl = '222 70 70 222 70 222 58 58 58'
	# Rotation Matrix
    euler_angle_1 = 3. 
    euler_angle_2 = 0.
    euler_angle_3 = 0.
  [../]
  [./Stiffness_carbide_gr1]
    type = ComputeElasticityTensor
    fill_method = symmetric9     #select symmetry
    base_name = C_carbide_gr1
    C_ijkl = '395 121 121 395 121 395 64.1 64.1 64.1'
    # Rotation Matrix
    euler_angle_1 = 3. 
    euler_angle_2 = 0.
    euler_angle_3 = 0.
  [../]    
  
  # Mechanical properties gr2
  [./Stiffness_matrix_gr2]
    type = ComputeElasticityTensor
    fill_method = symmetric9     #select symmetry
    base_name = C_matrix_gr2   #C11 C12 C13 C22 C23 C33 C44 C55 C66
    C_ijkl = '207.6 50.6 24.5 185.9	102.7 242.7 108.4 55.7 63.2' #T=250C
    # Rotation Matrix
    euler_angle_1 = 6. 
    euler_angle_2 = 0.
    euler_angle_3 = 0.
  [../]
  [./Stiffness_precipitate_gr2]
    type = ComputeElasticityTensor
    fill_method = symmetric9     #select symmetry
    base_name = C_precipitate_gr2  #C11 C12 C13 C22 C23 C33 C44 C55 C66
    C_ijkl = '222 70 70 222 70 222 58 58 58'
	# Rotation Matrix
    euler_angle_1 = 6. 
    euler_angle_2 = 0.
    euler_angle_3 = 0.
  [../]
  [./Stiffness_carbide_gr2]
    type = ComputeElasticityTensor
    fill_method = symmetric9     #select symmetry
    base_name = C_carbide_gr2
    C_ijkl = '395 121 121 395 121 395 64.1 64.1 64.1'
    # Rotation Matrix
    euler_angle_1 = 6. 
    euler_angle_2 = 0.
    euler_angle_3 = 0.
  [../]    
  
    # Mechanical properties gr3
  [./Stiffness_matrix_gr3]
    type = ComputeElasticityTensor
    fill_method = symmetric9     #select symmetry
    base_name = C_matrix_gr3   #C11 C12 C13 C22 C23 C33 C44 C55 C66
    C_ijkl = '207.6 50.6 24.5 185.9	102.7 242.7 108.4 55.7 63.2' #T=250C
    # Rotation Matrix
    euler_angle_1 = 9. 
    euler_angle_2 = 0.
    euler_angle_3 = 0.
  [../]
  [./Stiffness_precipitate_gr3]
    type = ComputeElasticityTensor
    fill_method = symmetric9     #select symmetry
    base_name = C_precipitate_gr3  #C11 C12 C13 C22 C23 C33 C44 C55 C66
    C_ijkl = '222 70 70 222 70 222 58 58 58'
	# Rotation Matrix
    euler_angle_1 = 9. 
    euler_angle_2 = 0.
    euler_angle_3 = 0.
  [../]
  [./Stiffness_carbide_gr3]
    type = ComputeElasticityTensor
    fill_method = symmetric9     #select symmetry
    base_name = C_carbide_gr3
    C_ijkl = '395 121 121 395 121 395 64.1 64.1 64.1'
    # Rotation Matrix
    euler_angle_1 = 9. 
    euler_angle_2 = 0.
    euler_angle_3 = 0.
  [../]   
  
    # Mechanical properties gr4
  [./Stiffness_matrix_gr4]
    type = ComputeElasticityTensor
    fill_method = symmetric9     #select symmetry
    base_name = C_matrix_gr4   #C11 C12 C13 C22 C23 C33 C44 C55 C66
    C_ijkl = '207.6 50.6 24.5 185.9	102.7 242.7 108.4 55.7 63.2' #T=250C
    # Rotation Matrix
    euler_angle_1 = 12. 
    euler_angle_2 = 0.
    euler_angle_3 = 0.
  [../]
  [./Stiffness_precipitate_gr4]
    type = ComputeElasticityTensor
    fill_method = symmetric9     #select symmetry
    base_name = C_precipitate_gr4  #C11 C12 C13 C22 C23 C33 C44 C55 C66
    C_ijkl = '222 70 70 222 70 222 58 58 58'
	# Rotation Matrix
    euler_angle_1 = 12. 
    euler_angle_2 = 0.
    euler_angle_3 = 0.
  [../]
  [./Stiffness_carbide_gr4]
    type = ComputeElasticityTensor
    fill_method = symmetric9     #select symmetry
    base_name = C_carbide_gr4
    C_ijkl = '395 121 121 395 121 395 64.1 64.1 64.1'
    # Rotation Matrix
    euler_angle_1 = 12. 
    euler_angle_2 = 0.
    euler_angle_3 = 0.
  [../]   
  
    # Mechanical properties gr5
  [./Stiffness_matrix_gr5]
    type = ComputeElasticityTensor
    fill_method = symmetric9     #select symmetry
    base_name = C_matrix_gr5   #C11 C12 C13 C22 C23 C33 C44 C55 C66
    C_ijkl = '207.6 50.6 24.5 185.9	102.7 242.7 108.4 55.7 63.2' #T=250C
    # Rotation Matrix
    euler_angle_1 = 15. 
    euler_angle_2 = 0.
    euler_angle_3 = 0.
  [../]
  [./Stiffness_precipitate_gr5]
    type = ComputeElasticityTensor
    fill_method = symmetric9     #select symmetry
    base_name = C_precipitate_gr5  #C11 C12 C13 C22 C23 C33 C44 C55 C66
    C_ijkl = '222 70 70 222 70 222 58 58 58'
	# Rotation Matrix
    euler_angle_1 = 15. 
    euler_angle_2 = 0.
    euler_angle_3 = 0.
  [../]
  [./Stiffness_carbide_gr5]
    type = ComputeElasticityTensor
    fill_method = symmetric9     #select symmetry
    base_name = C_carbide_gr5
    C_ijkl = '395 121 121 395 121 395 64.1 64.1 64.1'
    # Rotation Matrix
    euler_angle_1 = 15. 
    euler_angle_2 = 0.
    euler_angle_3 = 0.
  [../]   
  
  
  [./C]
    type = CompositeElasticityTensor
    args = 'eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    tensors = 'C_matrix_gr0   C_precipitate_gr0  C_carbide_gr0   C_matrix_gr1   C_precipitate_gr1  C_carbide_gr1 C_matrix_gr2   C_precipitate_gr2  C_carbide_gr2 C_matrix_gr3   C_precipitate_gr3  C_carbide_gr3 C_matrix_gr4   C_precipitate_gr4  C_carbide_gr4 C_matrix_gr5   C_precipitate_gr5  C_carbide_gr5'
    weights = 'h1_gr0  h2_gr0  h3_gr0   h1_gr1   h2_gr1   h3_gr1  h1_gr2   h2_gr2   h3_gr2   h1_gr3   h2_gr3   h3_gr3  h1_gr4   h2_gr4   h3_gr4  h1_gr5   h2_gr5   h3_gr5'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
  [./strain]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    eigenstrain_names = 'composite_eigenstrain'
  [../]
  
  [./eigen_precipitate_gr0]
    type = GenericConstantRankTwoTensor
    tensor_values = '0.1259  0.1259  0 0 0 0' #ei11,ei22,ei33,ei23,ei13,ei12
    tensor_name = eigen_precipitate_gr0 
  [../]   
  [./eigen_precipitate_gr1]
    type = GenericConstantRankTwoTensor
    tensor_values = '0.1259  0.1259  0 0 0 0' #ei11,ei22,ei33,ei23,ei13,ei12 
    tensor_name = eigen_precipitate_gr1 
  [../]   
  [./eigen_precipitate_gr2]
    type = GenericConstantRankTwoTensor
    tensor_values = '0.1259  0.1259  0 0 0 0' #ei11,ei22,ei33,ei23,ei13,ei12 
    tensor_name = eigen_precipitate_gr2 
  [../] 
  [./eigen_precipitate_gr3]
    type = GenericConstantRankTwoTensor
    tensor_values = '0.1259  0.1259  0 0 0 0' #ei11,ei22,ei33,ei23,ei13,ei12
    tensor_name = eigen_precipitate_gr3
  [../]   
  [./eigen_precipitate_gr4]
    type = GenericConstantRankTwoTensor
    tensor_values = '0.1259  0.1259  0 0 0 0' #ei11,ei22,ei33,ei23,ei13,ei12 
    tensor_name = eigen_precipitate_gr4 
  [../]   
  [./eigen_precipitate_gr5]
    type = GenericConstantRankTwoTensor
    tensor_values = '0.1259  0.1259  0 0 0 0' #ei11,ei22,ei33,ei23,ei13,ei12 
    tensor_name = eigen_precipitate_gr5 
  [../] 
  
  
  [./eigen_carbide_gr0]
    type = GenericConstantRankTwoTensor
    tensor_values = '0.0145 0.0145 0 0 0 0'   #ei11,ei22,ei33,ei23,ei13,ei12   
    tensor_name = eigen_carbide_gr0
  [../]
  [./eigen_carbide_gr1]
    type = GenericConstantRankTwoTensor
    tensor_values = '0.0145 0.0145 0 0 0 0'   #ei11,ei22,ei33,ei23,ei13,ei12  
    tensor_name = eigen_carbide_gr1
  [../]
  [./eigen_carbide_gr2]
    type = GenericConstantRankTwoTensor
    tensor_values = '0.0145 0.0145 0 0 0 0'   #ei11,ei22,ei33,ei23,ei13,ei12  
    tensor_name = eigen_carbide_gr2
  [../]
  [./eigen_carbide_gr3]
    type = GenericConstantRankTwoTensor
    tensor_values = '0.0145 0.0145 0 0 0 0'   #ei11,ei22,ei33,ei23,ei13,ei12   
    tensor_name = eigen_carbide_gr3
  [../]
  [./eigen_carbide_gr4]
    type = GenericConstantRankTwoTensor
    tensor_values = '0.0145 0.0145 0 0 0 0'   #ei11,ei22,ei33,ei23,ei13,ei12  
    tensor_name = eigen_carbide_gr4
  [../]
  [./eigen_carbide_gr5]
    type = GenericConstantRankTwoTensor
    tensor_values = '0.0145 0.0145 0 0 0 0'   #ei11,ei22,ei33,ei23,ei13,ei12  
    tensor_name = eigen_carbide_gr5
  [../]
  # [./eigen_bar]
    # type = GenericConstantRankTwoTensor
    # tensor_values = '-0.05 0 0 0 0 0'   #ei11,ei22,ei33,ei23,ei13,ei12  
    # tensor_name = eigen_bar
  # [../]
  [./composite_eigenstrain]
    type = CompositeEigenstrain
    tensors ='eigen_precipitate_gr0 eigen_precipitate_gr1 eigen_precipitate_gr2 eigen_precipitate_gr3 eigen_precipitate_gr4 eigen_precipitate_gr5 eigen_carbide_gr0 eigen_carbide_gr1 eigen_carbide_gr2 eigen_carbide_gr3 eigen_carbide_gr4 eigen_carbide_gr5'
    weights = 'h2_gr0 h2_gr1 h2_gr2 h2_gr3 h2_gr4 h2_gr5 h3_gr0 h3_gr1 h3_gr2 h3_gr3 h3_gr4 h3_gr5'
    args = 'eta1 eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    eigenstrain_name = 'composite_eigenstrain'
  [../]
  # [./bar]
    # type = DerivativeParsedMaterial
    # f_name = bar
    # function = '-1' 
  # [../]
[]

[Kernels]
  [./TensorMechanics]
  [../]

  #Kernels for diffusion equation
  [./diff_time]
    type = TimeDerivative
    variable = c
  [../]
  [./diff_c1]
    type = MatDiffusion
    variable = c
    diffusivity = Dh1
    v = c1
  [../]
  [./diff_c2]
    type = MatDiffusion
    variable = c
    diffusivity = Dh2
    v = c2
  [../]
  [./diff_c3]
    type = MatDiffusion
    variable = c
    diffusivity = Dh3
    v = c3
  [../]
  [./diff_gr0]
    type = MatDiffusion
    variable = c
    diffusivity = M_sgr0
    v = gr0
  [../]
  [./diff_gr1]
    type = MatDiffusion
    variable = c
    diffusivity = M_sgr1
    v = gr1
  [../]
  [./diff_gr2]
    type = MatDiffusion
    variable = c
    diffusivity = M_sgr2
    v = gr2
  [../]
  [./diff_gr3]
    type = MatDiffusion
    variable = c
    diffusivity = M_sgr3
    v = gr3
  [../]
  [./diff_gr4]
    type = MatDiffusion
    variable = c
    diffusivity = M_sgr4
    v = gr4
  [../]
  [./diff_gr5]
    type = MatDiffusion
    variable = c
    diffusivity = M_sgr5
    v = gr5
  [../] 
  
  # Kernels for Allen-Cahn equation of eta1
  [./deta1dt]
    type = TimeDerivative
    variable = eta1
  [../]
  [./ACElaeta1]
    type = AllenCahn
    variable = eta1
    f_name = f_el_mat
    args = 'eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    mob_name = L   
  [../]
  [./ACBulkFeta1]
    type = KKSMultiACBulkF
    variable  = eta1
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    gi_name   = g1
    eta_i     = eta1
    #wi        = 24       #double well height
    args      = 'c1 c2 c3 eta2 eta3'
    mob_name  = L
  [../]
  [./ACBulkCeta1]
    type = KKSMultiACBulkC
    variable  = eta1
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    cj_names  = 'c1 c2 c3'
    eta_i     = eta1
    args      = 'eta2 eta3'
    mob_name  = L
  [../]
  [./ACInterfaceeta1]
    type = ACInterface
    variable = eta1
    kappa_name = kappa
    mob_name = L
  [../]
  [./multipler1]
    type = MatReaction
    variable = eta1
    v = lambda
    mob_name = L
  [../]
 
  # Kernels for Allen-Cahn equation of eta2
  [./deta2dt]
    type = TimeDerivative
    variable = eta2
  [../]
  [./ACElaeta2]
    type = AllenCahn
    variable = eta2
    f_name = f_el_mat
    args = 'eta1 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    mob_name = L
  [../]
  [./ACBulkFeta2]
    type = KKSMultiACBulkF
    variable  = eta2
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    gi_name   = g2
    eta_i     = eta2
    #wi        = 24       #double well height
    args      = 'c1 c2 c3 eta1 eta3'
    mob_name  = L
  [../]
  [./ACBulkCeta2]
    type = KKSMultiACBulkC
    variable  = eta2
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    cj_names  = 'c1 c2 c3'
    eta_i     = eta2
    args      = 'eta1 eta3'
    mob_name  = L
  [../]
  [./ACInterfaceeta2]
    type = ACInterface
    variable = eta2
    kappa_name = kappa
    mob_name = L
  [../]
  [./multipler2]
    type = MatReaction
    variable = eta2
    v = lambda
    mob_name = L
  [../]
  
  
  # Kernels for Allen-Cahn equation of eta3
  [./deta3dt]
    type = TimeDerivative
    variable = eta3
  [../]
  [./ACElaeta3]
    type = AllenCahn
    variable = eta3
    f_name = f_el_mat
    args = 'eta1 eta2 gr0 gr1 gr2 gr3 gr4 gr5'
    mob_name = 0
  [../]
  [./ACBulkFeta3]
    type = KKSMultiACBulkF
    variable  = eta3
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    gi_name   = g3
    eta_i     = eta3
    #wi        = 24             #double well height
    args      = 'c1 c2 c3 eta1 eta2'
    mob_name  = 0
  [../]
  [./ACBulkCeta3]
    type = KKSMultiACBulkC
    variable  = eta3
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    cj_names  = 'c1 c2 c3'
    eta_i     = eta3
    args      = 'eta1 eta2'
    mob_name  = 0
  [../]
  [./ACInterfaceeta3]
    type = ACInterface
    variable = eta3
    kappa_name = kappa
    mob_name = 0
  [../]
  [./multipler3]
    type = MatReaction
    variable = eta3
    v = lambda
    mob_name = 0
  [../]
  
  # Kernels for the Lagrange multiplier equation
  [./mult_lambda]
    type = MatReaction
    variable = lambda
    mob_name = 2
  [../]
  [./mult_ACEla_1]
    type = CoupledAllenCahn
    variable = lambda
    v = eta1
    f_name = f_el_mat
    args = 'eta2 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    mob_name = 1
  [../]
  [./mult_ACBulkF_1]
    type = KKSMultiACBulkF
    variable  = lambda
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    gi_name   = g1
    eta_i     = eta1
    #wi        = 24                #double well height
    mob_name  = 1
    args      = 'c1 c2 c3 eta2 eta3'
  [../]
  [./mult_ACBulkC_1]
    type = KKSMultiACBulkC
    variable  = lambda
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    cj_names  = 'c1 c2 c3'
    eta_i     = eta1
    args      = 'eta2 eta3'
    mob_name  = 1
  [../]
  [./mult_CoupledACint_1]
    type = SimpleCoupledACInterface
    variable = lambda
    v = eta1
    kappa_name = kappa
    mob_name = 1
  [../]
  
  
  [./mult_ACEla_2]
    type = CoupledAllenCahn
    variable = lambda
    v = eta2
    f_name = f_el_mat
    args = 'eta1 eta3 gr0 gr1 gr2 gr3 gr4 gr5'
    mob_name = 1
  [../]
  [./mult_ACBulkF_2]
    type = KKSMultiACBulkF
    variable  = lambda
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    gi_name   = g2
    eta_i     = eta2
    #wi        = 24                #double well height
    mob_name  = 1
    args      = 'c1 c2 c3 eta1 eta3'
  [../]
  [./mult_ACBulkC_2]
    type = KKSMultiACBulkC
    variable  = lambda
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    cj_names  = 'c1 c2 c3'
    eta_i     = eta2
    args      = 'eta1 eta3'
    mob_name  = 1
  [../]
  [./mult_CoupledACint_2]
    type = SimpleCoupledACInterface
    variable = lambda
    v = eta2
    kappa_name = kappa
    mob_name = 1
  [../]
  
  [./mult_ACEla_3]
    type = CoupledAllenCahn
    variable = lambda
    v = eta3
    f_name = f_el_mat
    args = 'eta1 eta2 gr0 gr1 gr2 gr3 gr4 gr5'
    mob_name = 0
  [../]
  [./mult_ACBulkF_3]
    type = KKSMultiACBulkF
    variable  = lambda
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    gi_name   = g3
    eta_i     = eta3
    #wi        = 24             #double well height
    mob_name  = 0
    args      = 'c1 c2 c3 eta1 eta2'
  [../]
  [./mult_ACBulkC_3]
    type = KKSMultiACBulkC
    variable  = lambda
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    cj_names  = 'c1 c2 c3'
    eta_i     = eta3
    args      = 'eta1 eta2'
    mob_name  = 0
  [../]
  [./mult_CoupledACint_3]
    type = SimpleCoupledACInterface
    variable = lambda
    v = eta3
    kappa_name = kappa
    mob_name = 0
  [../]
 
  # Phase concentration constraints (KKS ASSUMPTION)
  [./chempot12]
    type = KKSPhaseChemicalPotential
    variable = c1
    cb       = c2
    fa_name  = F1
    fb_name  = F2
  [../]
  [./chempot23]
    type = KKSPhaseChemicalPotential
    variable = c2
    cb       = c3
    fa_name  = F2
    fb_name  = F3
  [../]
  [./phaseconcentration]
    type = KKSMultiPhaseConcentration
    variable = c3
    cj = 'c1 c2 c3'
    hj_names = 'h1 h2 h3'
    etas = 'eta1 eta2 eta3'
    c = c
  [../]
  
  # Kernels for Allen-Cahn equation of grain
  [./gr0]
    type = TimeDerivative
    variable = gr0
  [../]
  [./gr1]
    type = TimeDerivative
    variable = gr1
  [../]
  [./gr2]
    type = TimeDerivative
    variable = gr2
  [../]
  [./gr3]
    type = TimeDerivative
    variable = gr3
  [../]
  [./gr4]
    type = TimeDerivative
    variable = gr4
  [../]
  [./gr5]
    type = TimeDerivative
    variable = gr5
  [../]
  
  # Kernels for Allen-Cahn equation of grain
  [./in_stress_11]
    type = TimeDerivative
    variable = in_stress_11
  [../]
  [./in_stress_12]
    type = TimeDerivative
    variable = in_stress_12
  [../]
  [./in_stress_22]
    type = TimeDerivative
    variable = in_stress_22
  [../]
[]

[Executioner]
  type = Transient
  scheme = 'BDF2'
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm      ilu          nonzero'

  l_max_its = 30       #default 10000
  nl_max_its = 20      #default 50
  
  l_tol = 1e-4         #default 1e-05
  nl_abs_tol = 1e-03   #default 1e-50
  nl_rel_tol = 1e-03   #default 1e-08
  
  #num_steps = 20
  end_time = 0.1

  [./TimeStepper] 
    type = IterationAdaptiveDT 
    dt = 0.0001          
    cutback_factor = 0.8  #default 0.5
    growth_factor = 1.5   #default 2.0
    optimal_iterations = 10
    #timestep_limiting_postprocessor = dtnuc
  [../] 
[]

# [Adaptivity]
  # [./Indicators]
    # [./jump]
      # type = GradientJumpIndicator
      # variable = eta2
    # [../]
	# [./GJI_gr0]
      # type = GradientJumpIndicator
      # variable = gr0
    # [../]
    # [./GJI_gr1]
      # type = GradientJumpIndicator
      # variable = gr1
    # [../]
	# [./GJI_gr2]
      # type = GradientJumpIndicator
      # variable = gr2
    # [../]
    # [./GJI_gr3]
      # type = GradientJumpIndicator
      # variable = gr3
    # [../]
    # [./GJI_gr4]
      # type = GradientJumpIndicator
      # variable = gr4
    # [../]
	# [./GJI_gr5]
      # type = GradientJumpIndicator
      # variable = gr5
    # [../]
  # [../]
  # [./Markers]
    # [./nuc]
      # type = DiscreteNucleationMarker
      # map = map
    # [../]
    # [./grad_gr0]
      # type = ErrorFractionMarker
      # coarsen = 0.2
      # refine = 0.9
      # indicator = GJI_gr0
    # [../]
	# [./grad_gr1]
      # type = ErrorFractionMarker
      # coarsen = 0.2
      # refine = 0.9
      # indicator = GJI_gr1
    # [../]
    # [./grad_gr2]
      # type = ErrorFractionMarker
      # coarsen = 0.2
      # refine = 0.9
      # indicator = GJI_gr2
    # [../]
	# [./grad_gr3]
      # type = ErrorFractionMarker
      # coarsen = 0.2
      # refine = 0.9
      # indicator = GJI_gr3
    # [../]
	# [./grad_gr4]
      # type = ErrorFractionMarker
      # coarsen = 0.2
      # refine = 0.9
      # indicator = GJI_gr4
    # [../]
	# [./grad_gr5]
      # type = ErrorFractionMarker
      # coarsen = 0.2
      # refine = 0.9
      # indicator = GJI_gr5
    # [../]
	# [./grad]
      # type = ValueThresholdMarker
      # variable = jump
      # coarsen = 0.1
      # refine = 0.2
    # [../]
    # [./combo]
      # type = ComboMarker
      # markers = 'nuc grad_gr0 grad_gr1 grad_gr2 grad_gr3 grad_gr4 grad_gr5' # grad  
    # [../]
  # [../]
  # marker = combo
  # cycles_per_step = 1
  # recompute_markers_during_cycles = true
  # max_h_level = 1
# []


# Precondition using handcoded off-diagonal terms
#
[Preconditioning]
  [./full]
    type = SMP
    full = true
  [../]
[]


[Postprocessors]
  [./aver_vonmises]
    type = ElementAverageValue
    variable = vonmises
	execute_on = 'initial timestep_end'
  [../]
  [./aver_eta2]
    type = ElementAverageValue
    variable = eta2
	execute_on = 'initial timestep_end'
  [../]
  [./dt]
    type = TimestepSize
  [../]
  [./rate]
    type = DiscreteNucleationData
    value = RATE
    inserter = inserter
  [../]
  # [./dtnuc]
    # type = DiscreteNucleationTimeStep
    # inserter = inserter
    # p2nucleus = 0.03
    # dt_max = 0.01
  # [../]
  [./update]
    type = DiscreteNucleationData
    value = UPDATE
    inserter = inserter
  [../]
  [./count]
    type = DiscreteNucleationData
    value = COUNT
    inserter = inserter
  [../]
  [./nuc_insertions]
    type = DiscreteNucleationData
    inserter = inserter
    value = INSERTIONS
  [../]
  [./nuc_deletions]
    type = DiscreteNucleationData
    inserter = inserter
    value = DELETIONS
  [../]
[]

[Debug]               
  show_var_residual_norms = true 
[]

[Outputs]
  exodus = true
  [./csv]
    type = CSV
    execute_on = 'initial timestep_end'
  [../]
[]
