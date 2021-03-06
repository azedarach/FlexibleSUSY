Block MODSEL                 # Select model
    1   0                    # mSUGRA
    3   1                    # NMSSM
    6   0                    # flavour violation
   12   750		     # parameter input scale
Block FlexibleSUSY
    0   1.000000000e-04      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = two_scale, 1 = lattice)
    3   0                    # calculate SM pole masses
    4   2                    # pole mass loop order
    5   2                    # EWSB loop order
    6   2                    # beta-functions loop order
    7   2                    # threshold corrections loop order
    8   1                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   1                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   1                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   1                    # Higgs 2-loop corrections O(alpha_tau^2)
Block SOFTSUSY               # SOFTSUSY specific inputs
    1   1.000000000e-04      # tolerance
    2   2.000000000e+00      # up-quark mixing (=1) or down (=2)
    5   1.000000000E+00      # 2-loop running
    3   0.000000000E+00      # printout
    4   750                  # parameter output scale
    7   4.0		     # number of Higgs loops
   18   1.000000000E+00      # use soft Higgs masses as EWSB output
BLOCK SMINPUTS
    1	128.962	 	 # 1/alpha_em(MZ)
    2	0.000011663900000000002	 	 # GF
    3	0.1172	 	 # alpha_s(MZ)
    4	91.1876	 	 # MZ^pole
    5	4.2	 	 # mb(mb)^SM MS-bar
    6	172.9	 	 # mt^pole
    7	1.777	 	 # mtau^pole
    9	80.385	 	 # MW^pole
   11	0.00051099891	 # m_e
   13	1.10565836	 # m_mu
   21	0.00495	 	 # m_d
   22	0.0025	 	 # m_u
   23	0.1	 	 # m_s
   24	1.42	 	 # m_c
BLOCK EXTPAR
    0	750.	 	 # Qin
    1	120.	 	 # M1
    2	200.	 	 # M2
    3	1500.	 	 # M3
   11	1000.	 	 # At
   12	1000.	 	 # Ab
   13	1000.	 	 # Atau
   25	2.0	 	 # TanBeta
   31	1500.	 	 # ml1
   32	1500.	 	 # ml2
   33	1500.	 	 # ml3
   34	1500.	 	 # me1
   35	1500.	 	 # me2
   36	1500.	 	 # me3
   41	1500.	 	 # mq1
   42	1500.	 	 # mq2
   43	750.	 	 # mq3
   44	1500.	 	 # mu1
   45	1500.	 	 # mu2
   46	750.	 	 # mu3
   47	1500.	 	 # md1
   48	1500.	 	 # md2
   49	1500.	 	 # md3
   61	0.67	 	 # Lambda
   62	0.2	 	 # Kappa
   63	405.	 	 # ALambda
   64	-0.	 	 # AKappa
   65	200.	 	 # MuEff
