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

def RK4integrate(h0, v0, mwet, mdry, dragcoeff, dt, maxtime, enginelist, thrustlevel, lat=0, fn="log.txt"):
  mass = mwet
  tt = 0
  h1 = h0
  v1 = v0
  thrustisp = TotalThrust(enginelist, KerbinPressure(h1), thrustlevel)
  effthrust = thrustisp[0]
  exhvel = thrustisp[1]*Kerbing0
  fuelflow = effthrust/exhvel
  tburn = (mwet-mdry)/fuelflow
  
  drag = DragForce(KerbinPressure(h0), v0, dragcoeff)
  geff = KerbinGravConst(h0,lat) - drag/(3*mass)
  maxa = effthrust/mass

#initial guess for t_s, h_s  
  ts = (geff - maxa + sqrt( (geff - maxa)**2 - 2*v0*maxa**2/exhvel ) ) * exhvel/maxa**2
  
  hs = exhvel*(mass - fuelflow*ts) / fuelflow * log( mass / (mass - fuelflow*ts) ) - (v0 + exhvel)*ts + ( KerbinGravConst(h0,lat) - drag/(2*mass) )*ts**2/2

#need to correct for Isp changes with height  
  h_ave = h0 - 2.0/3.0*hs
  
  avethrustisp = TotalThrust(enginelist, KerbinPressure(h_ave), thrustlevel)
  effthrust = avethrustisp[0]
  exhvel = avethrustisp[1]*Kerbing0
  geff = KerbinGravConst(h_ave,lat) - drag/(3*mass)

#correction for t_s, h_s
  ts = ( v0 - exhvel * log( 1 - fuelflow * ts / mass ) ) / geff
  
  hs = exhvel*(mass - fuelflow*ts) / fuelflow * log( mass / (mass - fuelflow*ts) ) - (v0 + exhvel)*ts + ( KerbinGravConst(h_ave,lat) - drag/(2*mass) )*ts**2/2
  
  outfile = open(fn,"wt")
  
  outfile.write("Guess of t_s, h_s:\n")
  outfile.write("t_s = " + '{0:.2f}'.format(ts) + " s;  h_s = " + '{0:.1f}'.format(hs) + " m\n")
  
  outfile.write("Time,s // Distance traveled // H, m // VSpeed, m/s // Isp, s\n")
  
  while (tt <= maxtime) and (v1 < 0):
    outfile.write('{0:.3f}'.format(tt) + "  " + '{}'.format(h0-h1) + "  " + '{}'.format(h1) + "  " + '{}'.format(v1) + "  " + repr(TotalThrust(enginelist, KerbinPressure(h1))[1]) + "\n")
		
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
InitialVspeed = -160
InitialMass = 23000
DryMass = 15800
DragArea = 9.8
Timestep = 0.01
MaxIntegrationTime = 1000
InitLatitude = 90

LogFile = "droplog_" + '{0:.0f}'.format(InitialHeight) + "_" + '{0:.0f}'.format(-InitialVspeed) + ".txt"

RK4integrate(InitialHeight, InitialVspeed, InitialMass, DryMass, DragArea, Timestep, MaxIntegrationTime, EngineList, UserThrustLevel, InitLatitude, LogFile)
