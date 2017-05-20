#!/usr/bin/python3.5
#RK4 flight integration

from math import *

#useful constants
KerbinMu = 3.5316e12 #gravitational parameter
KerbinR = 6e5
KerbinDay = 21549.425
Kerbing0 = KerbinMu/KerbinR**2
KerbinAtm = 101.325
KerbinAtmScale1 = 6930.5
KerbinAtmScale2 = 20240.0

class Engine:
  def __init__(self, vacIsp, slIsp, maxthrustvac):
    self.ispv = vacIsp
    self.isps = slIsp
    self.maxff = maxthrustvac/vacIsp/Kerbing0

  def thrust(self, press = 0, thrustlevel = 1.0, thrustlim = 1.0):
    if press > 0:
      return Kerbing0*thrustlevel*thrustlim*self.maxff*(self.ispv + (self.isps-self.ispv)*press/KerbinAtm)
    else:
      return Kerbing0*thrustlevel*thrustlim*self.maxff*self.ispv

class EngineMount:
  def __init__(self, Engn, attAngle, thrustlim = 1.0):
    "Input attAngle in deg, self.angle in radians"
    self.engine = Engn
    self.angle = attAngle/180.0*3.1415927
    self.englimit = thrustlim

EngineDic = {'Reliant': Engine(310,265,240e3), 'Swivel': Engine(320,250,215e3), 'Thud': Engine(305,275,120e3), 'Spark': Engine(320,270,20e3), 'Twitch': Engine(290,250,16e3), 'Terrier': Engine(345,85,60e3), 'Vector': Engine(315,295,1.0e6), 'Skipper': Engine(320,280,650e3), 'Mainsail': Engine(310,285,1.5e6), 'TwinBoar': Engine(300,280,2.0e6), 'Poodle': Engine(350,90,250e3), 'Rhino': Engine(340,205,2.0e6), 'Mammoth': Engine(315,295,4.0e6), 'VestaVR1': Engine(335,260,90e3)}

def KerbinGravConst(h,lat=0):
  "Calculate gravity acceleration at Kerbin, height h in meters, latitude in degrees"
  return KerbinMu/(KerbinR+h)**2 - (2*3.1415927/KerbinDay)**2*KerbinR*cos(lat/180.0*3.1415927)

def KerbinPressure(h):
  if h < 70000:
    return KerbinAtm*exp(-h/KerbinAtmScale1 - (h/KerbinAtmScale2)**2)
  else:
    return 0

def DragForce(p, v, dragcoeff):
  return dragcoeff*p*29/305/8.314*v**2/2
	
def VesselMass(mwet,mdry,fuelflow,t):
  return max(mdry,mwet-fuelflow*t)
  
def TotalThrust(EngMntList,press,thrustlevel = 1.0):
  "First output value is thrust in Newtons, second is effective Isp in seconds"
  sumthrust = 0
  sumff = 0
  for em in EngMntList:
    sumthrust += em.engine.thrust(press, thrustlevel, em.englimit)*cos(em.angle)
    sumff += em.engine.maxff * thrustlevel
  return [ sumthrust, sumthrust/sumff/Kerbing0 ]

def dvexdh(EngMntList, hasl, press, pressasl):
  "Returns du/dh at height Hasl above sea level"
  if press == 0:
    return 0
  hscale = hasl/log(pressasl/press) #assuming p=p0exp(-h/hscale)
  vexh = TotalThrust(EngMntList, press)[1]
  vexsl = TotalThrust(EngMntList, pressasl)[1]
  dvexdp = (vexh - vexsl)/(press - pressasl)*Kerbing0
  
  return -press/hscale*dvexdp

def intlog(x):
  "Returns indefinite integral of ln(x)"
  return x*(log(x) - 1)
  
def RK4integrate(h0, v0, mwet, mdry, dragcoeff, dt, maxtime, enginelist, thrustlevel, lat=0, fn="log", OutputAsCSV=False):
  if OutputAsCSV:
    filesuffix = ".csv"
    delimiter = "; "
  else:
    filesuffix = ".txt"
    delimiter = "  "    
  mass = mwet
  tt = 0
  h1 = h0
  v1 = v0
  
  dgdh = -2*KerbinMu/(KerbinR + h0)**3 #gradient of free fall acceleration
  thrustisp = TotalThrust(enginelist, KerbinPressure(h0), thrustlevel)
  effthrust = thrustisp[0]
  exhvel = thrustisp[1]*Kerbing0
  fuelflow = effthrust/exhvel
  tburn = (mwet-mdry)/fuelflow
  
  drag = DragForce(KerbinPressure(h0), v0, dragcoeff)
  geff = KerbinGravConst(h0,lat) - drag/(3*mass)
  maxa = effthrust/mass
  dudh = dvexdh(enginelist, h0, KerbinPressure(h0), KerbinAtm)

  v_lin_term = maxa - geff
  v_sq_term = maxa**2/(exhvel*2) + 1.0/6.0*fuelflow/mass*dudh*v0 - dgdh*v0 / 3.0
  
  print(repr(dgdh))

  ts = (-v_lin_term + sqrt( v_lin_term**2 - 4*v0*v_sq_term ) ) / (2*v_sq_term)
  
  #dv/dt = (v_e[0] + v_e[1]t + v_e[2]t^2)/(v_d - t) + A*(v0 - v0/ts*t)^2 - (g + dg/dh*v0*t - v0/ts*dg/dh*t^2/2)
  v_e = [ exhvel, dudh*v0, -dudh*v0/(2*ts) ]
  v_d = mass/fuelflow
  
  dv_frac_const = -(v_e[1] + v_e[2]*v_d)
  dv_frac_lin = -v_e[2] - dgdh*v0
  dv_frac_e = v_e[0] - v_d * dv_frac_const
  
  hs = -v0*ts + (KerbinGravConst(h0,lat) - drag/(2*mass) - dv_frac_const)*ts**2/2 - dv_frac_e * v_d *(1 + intlog(1 - ts/v_d)) - dv_frac_lin/6*ts**3
  
  outfile = open(fn+filesuffix,"wt")
  
  outfile.write("Guess of t_s, h_s:\n")
  outfile.write("t_s = " + '{0:.2f}'.format(ts) + " s" + delimiter + "h_s = " + '{0:.1f}'.format(hs) + " m\n")
  
  outfile.write("Time,s" + delimiter + "Distance traveled" + delimiter + " H, m" + delimiter + "VSpeed, m/s" + delimiter + "Isp, s\n")
  
  while (tt <= maxtime) and (v1 < 0):
    outfile.write('{0:.3f}'.format(tt) + delimiter + '{}'.format(h0-h1) + delimiter + '{}'.format(h1) + delimiter + '{}'.format(v1) + delimiter + repr(TotalThrust(enginelist, KerbinPressure(h1))[1]) + "\n")
		
    k1a = (TotalThrust(enginelist, KerbinPressure(h1), thrustlevel)[0]*max(0,copysign(1,tburn-tt))  - DragForce(KerbinPressure(h1), v1, dragcoeff)*copysign(1,v1))/VesselMass(mwet,mdry,fuelflow,tt) - KerbinGravConst(h1,lat)
    k1v = v1
    
    h2 = h1 + k1v*dt*0.5
    v2 = v1 + k1a*dt*0.5
    
    k2a = (TotalThrust(enginelist, KerbinPressure(h2), thrustlevel)[0]*max(0,copysign(1,tburn-tt-dt*0.5)) - DragForce(KerbinPressure(h2), v2, dragcoeff)*copysign(1,v2))/VesselMass(mwet,mdry,fuelflow,tt+dt*0.5) - KerbinGravConst(h2,lat)
    k2v = v2
    
    h3 = h1 + k2v*dt*0.5
    v3 = v1 + k2a*dt*0.5
    
    k3a = (TotalThrust(enginelist, KerbinPressure(h3), thrustlevel)[0]*max(0,copysign(1,tburn-tt-dt*0.5)) - DragForce(KerbinPressure(h3), v3, dragcoeff)*copysign(1,v3))/VesselMass(mwet,mdry,fuelflow,tt+dt*0.5) - KerbinGravConst(h3,lat)
    k3v = v3
    
    h4 = h1 + k3v*dt
    v4 = v1 + k3a*dt
    
    k4a = (TotalThrust(enginelist, KerbinPressure(h4), thrustlevel)[0]*max(0,copysign(1,tburn-tt-dt)) - DragForce(KerbinPressure(h4), v4, dragcoeff)*copysign(1,v4))/VesselMass(mwet,mdry,fuelflow,tt+dt) - KerbinGravConst(h4,lat)
    k4v = v4
    
    h1 = h1 + dt/6.0*(k1v+2*k2v+2*k3v+k4v)
    v1 = v1 + dt/6.0*(k1a+2*k2a+2*k3a+k4a)
    
    tt=tt+dt
  
  outfile.close()
  return 0

EngineList = [EngineMount(EngineDic['Thud'],22.0719225), EngineMount(EngineDic['Thud'],22.0719225), EngineMount(EngineDic['Thud'],22.0719225), EngineMount(EngineDic['Thud'],22.0719225)]

UserThrustLevel = 0.75
InitialHeight = 3500 #above sea level
InitialVspeed = -60
InitialMass = 23000
DryMass = 15800
DragArea = 9.8
Timestep = 0.01
MaxIntegrationTime = 1000
InitLatitude = 90

LogFile = "droplog_1_" + '{0:.0f}'.format(InitialHeight) + "_" + '{0:.0f}'.format(-InitialVspeed)

RK4integrate(InitialHeight, InitialVspeed, InitialMass, DryMass, DragArea, Timestep, MaxIntegrationTime, EngineList, UserThrustLevel, InitLatitude, LogFile,True)
