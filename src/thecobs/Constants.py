from scipy import constants

######################
# Declaring constants

ELECTRON_MASS_MEV = constants.m_e/constants.eV*constants.speed_of_light**2/1E6
ELECTRON_MASS_KEV = ELECTRON_MASS_MEV*1E3
ELECTRON_MASS_EV = ELECTRON_MASS_MEV*1E6
PROTON_MASS_EV = constants.m_p/constants.eV*constants.speed_of_light**2
PROTON_MASS_KEV = PROTON_MASS_EV/1E3
NUCLEON_MASS_EV = (constants.m_p+constants.m_n)/2.0/constants.eV*constants.speed_of_light**2
NUCLEON_MASS_KEV = NUCLEON_MASS_EV/1E3
NUCLEON_MASS_W = NUCLEON_MASS_EV/(1.0*ELECTRON_MASS_EV)
PROTON_MASS_W = PROTON_MASS_EV/(1.0*ELECTRON_MASS_EV)
AMU_MASS_KEV = constants.physical_constants['atomic mass constant energy equivalent in MeV'][0]*1E3
NATURAL_LENGTH = (constants.hbar*constants.speed_of_light/constants.e)/ELECTRON_MASS_EV
ALPHA = constants.alpha

####################
# Nucleon coupling constants

GA = 1.27
GM = 4.706
