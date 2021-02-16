# -*- coding: utf-8 -*-

# Program to generate and solve the polynomial pedotransfer function from
# article (name or DOI here).

# You may choose the parameters and the program will build and solve the polynomial
# The parameters are availabe on "parameters" dictionary below.

# Some VARIABLES change from the journal article version.
# al   : alpha parameter
# lamb : l parameter related with conductivity function
# dep  : depth of the considered boudary layer


# ------------------------------------------------------------------------------------------------

parameters = {'al': 0.4624 ,
              'n': 2.108  ,
              'lamb': 0.601  ,
              'Ks': 41.98  ,
              'dep': 30  ,
              'q': 3.6  ,
              }



# ------------------------------------------------------------------------------------------------


from sympy.parsing.sympy_parser import parse_expr
from sympy import N

ctn = [0.8654654556811963, -0.05434726070081017, 0.08037116529839043, 0.008762800602361308,
       -0.04666291216858057, 0.2622780247098996, -0.09413191609868982, 0.06132681005269242, 
       -0.13323206575737548, -0.021039571391720215, 0.012148066143768961, -0.02987773097651346, 
       0.022671446739094357, 0.03335805429295968, 0.020147726195870116, 1.6279897446766827, 
       -1.4530569552105415, 0.7043355004062589, -0.541784083894665, 0.5902799985238507,
       0.09038263602784666, -0.24899715554877372, -0.10399563657752278, 0.21908354200635233,
       -0.1434519462819891, -0.030160399742816372, -0.042722279074579816, 0.12831152068623008,
       0.027436721757170452, -0.1565433350511172, -0.03806968196564138, -0.19466015404836676, 
       0.008474142129239706, -0.2914427928534872, 0.05235636114710501, -0.018086620387045928, 
       0.026784035410313163, -0.07537441951881563, 0.028860757171912103, -0.06277529823160982,
       0.003890244066777326, 0.02100705109908848, 0.030434602257902094, 0.35107292818654534,
       -0.03178639761011755, -0.6251196942516541, 0.1907936111626416, 0.03710196881236187, 
       0.09655013764417436, -0.0074343753768419465]

varia = ['al', 'al*al', 'al*al*al', 'al*al*1/n', 'al*1/n*1/n', 'al*1/n*Ks', 'al*1/n*q', 
         'al*1/n*dep', 'al*1/n*lamb', 'al*Ks*Ks', 'al*Ks*q', 'al*Ks*dep', 'al*dep*dep', 
         'al*lamb', '1/n*1/n', '1/n*1/n*1/n', '1/n*1/n*Ks', '1/n*1/n*q', '1/n*1/n*dep', 
         '1/n*1/n*lamb', '1/n*Ks', '1/n*Ks*Ks', '1/n*Ks*q', '1/n*Ks*dep', '1/n*Ks*lamb', 
         '1/n*q*q', '1/n*q*dep', '1/n*q*lamb', '1/n*dep*dep', '1/n*dep*lamb', '1/n*lamb',
         '1/n*lamb*lamb', 'Ks', 'Ks*Ks', 'Ks*Ks*q', 'Ks*Ks*dep', 'Ks*q', 'Ks*q*q', 'Ks*q*dep',
         'Ks*q*lamb', 'Ks*dep*dep', 'Ks*lamb', 'q', 'q*lamb', 'dep', 'dep*dep', 'dep*lamb',
         'lamb', 'lamb*lamb']


ctn05 = [0.8131647710603545, -0.1102031672431125, 0.09506884153342532, 0.011455932865381073,
        -0.04965007316753883, 0.28382412107956645, -0.1050035997168321, 0.07282709679978056,
        -0.1589842741168267, 0.014175516425969018, -0.03494374375313618, 0.02476540402469592,
        0.1149747290807096, 1.3989426440748427, -1.2155289954940884, 0.6761182696846543,
        -0.5271869263673907, 0.5580594840035178, -0.1959168298025158, -0.12298540966949786,
        0.22652698068648264, -0.14982734551570032, -0.09646615441904371, 0.13445733734620835,
        -0.1626018344970105, -0.3051626604946473, 0.06130721249420479, -0.021003997737688032,
        0.03281140582683211, -0.07069557592938963, 0.024012975895110422, -0.06854300165868961,
        0.019207769964987353, 0.33815206287610744, 0.02943941769713415, -0.4789542530153339, 
        0.1527851846068795]

varia05 = ['al', 'al*al', 'al*al*al', 'al*al*1/n', 'al*1/n*1/n', 'al*1/n*Ks', 'al*1/n*q', 
           'al*1/n*dep', 'al*Ks*Ks', 'al*Ks*q', 'al*Ks*dep', 'al*dep', '1/n*1/n', '1/n*1/n*1/n', 
           '1/n*1/n*Ks', '1/n*1/n*q', '1/n*1/n*dep', '1/n*Ks', '1/n*Ks*Ks', '1/n*Ks*q', 
           '1/n*Ks*dep', '1/n*q*q', '1/n*q*dep', '1/n*dep*dep', 'Ks', 'Ks*Ks', 'Ks*Ks*q', 
           'Ks*Ks*dep', 'Ks*q', 'Ks*q*q', 'Ks*q*dep', 'Ks*dep*dep', 'q', 'q*q*dep', 'dep', 
           'dep*dep']

def poli_m(ctn, varia):
    "Polynomial Maker"
    varia.insert(0, '1')
    poli = float(0)
    if len(varia) != len(ctn): print('!!! Problem, different number of predictors and coefficients !!!')
    
   # Remover variaveis    
    for i in range(len(varia)):
        
        varia[i] = varia[i].replace('al', 'log(al, 10)').replace('Ks', 'log(Ks, 10)').replace('dep',
        'log(dep, 10)').replace('q', 'log(q, 10)').replace('1/n','(1/n)')

        aa = parse_expr(varia[i])
        poli = poli + aa * ctn[i]
 
    return poli

if parameters['lamb'] == 0.5:
  var = (ctn05, varia05)
else:
  var = (ctn, varia)

a1 = poli_m(var[0], var[1]).subs(parameters)

ans = float(N(a1.subs(parameters)))
if ans > 1.0:
    ans = 1.0
elif ans < 0.0:
    ans = 0.01
print('The solution for the polinomial with proposed parameters is:')
print(ans)








