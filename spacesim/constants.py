"""Space related constants.
"""
import scipy.constants as spconst
from datetime import datetime
                                        				# Description								Unit
                                        				# -------------------------------------------------------
M_EARTH : float = 5.97219e24                    		# Mass of Earth								kg
MU_EARTH: float = spconst.G * M_EARTH             		# Gravitational parameter of Earth			m^3 / s^2
R_EQAT_EARTH: float = 6378137.0							# Equatorial radius of Earth				m
R_POLE_EARTH: float = 6356752.0							# Polar radius of Earth						m
R_EARTH: float = R_EQAT_EARTH                  			# Nominal radius of Earth					m
e2_EARTH: float = 6.69437999014e-3						# Eccentricty of earth ellipsoid squared	-
ROT_V_EARTH: float = 7.2921159e-5               		# Rotation rate of Earth					rad / s
J2000: datetime = datetime(2000, 1, 1, 12, 0, 0) 		# J2000 epoch								UTC
J2000_ERA: float = 4.89496121272365 			  		# Earth rotation angle at J2000 epoch		rad
sidereal_day_length: float = 86164.1					# Sidereal day length						seconds
solar_day_length: float = 86400.0						# Solar day length							seconds
J2_EARTH: float = 1.082626925638815e-3					# J2 coefficient of Earth					-
TU: float = 806.8116									# Earth Time Unit							seconds
C_ZERO_K: float = 273.15								# Zero degrees Celsius in Kelvin			K
D_EARTH: float = 86164.10035							# Sidereal rotation period of the Earth		seconds