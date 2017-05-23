global Kerbing0 to 9.80665.

function limiter {
  parameter val, minval, maxval.
  return min( maxval, max( val, minval ) ).
}

function teleprint {
  parameter pid.
  print "Error: " + round(pid:error,2) at (0,5).
  print "PTerm: " + round(pid:pterm,2) at (0,6).
  print "ITerm: " + round(pid:iterm,2) at (0,7).
  print "DTerm: " + round(pid:dterm,2) at (0,8).
  print "Output: " + round(pid:output,2) at (0,9).
}

function intlog {
  parameter x.
  return x*(ln(x) - 1).
}

function BrakeDistance {
  parameter hasl.
  parameter vspd.
  parameter vt.
  parameter tl.
  parameter geff to Kerbing0.
  parameter p to 1. // in Kerbin atmospheres
  parameter drag to 0.
  
  local maxa to vt[0]*tl/mass.
  
  local hscale to hasl/ln(slp/p).
  local vexh to vt[1].
  local slvexh to vthrustsl[1].
  local dudh to p*(vexh - slvexh)/((slp - p)*hscale).
  local dgdh to -2*muhere/(body:radius + hasl)^3.
  local ff to vt[0]*tl/vexh.

  local v_lin_term to maxa - geff + drag / (3.0*mass).
  local v_sq_term to maxa*maxa/(vexh*2) + 1.0/6.0*ff/mass*dudh*vspd - dgdh*vspd / 3.0.
  
  local ts to (-v_lin_term + sqrt( v_lin_term*v_lin_term - 4*vspd*v_sq_term ) ) / (2*v_sq_term).
  
  local v_e to list( vexh, dudh*vspd, -dudh*vspd/(2*ts) ).
  local v_d to mass/ff.

  local dv_frac_const to -(v_e[1] + v_e[2]*v_d).
  local dv_frac_lin to -v_e[2] - dgdh*vspd.
  local dv_frac_e to v_e[0] - v_d * dv_frac_const.
  
  local hs to -vspd*ts + (geff - drag/(2*mass) - dv_frac_const)*ts*ts/2 - dv_frac_e * v_d *(1 + intlog(1 - ts/v_d)) - dv_frac_lin/6*ts^3.
  return hs.
}

function updvthrust {
  parameter vthrust.
  parameter p.
  list engines in el.
  local tt to 0.
  local tff to 0.
  for eng in el {
    local tmax to eng:availablethrustat(p).
    if tmax > 0 {
      set tt to tt + tmax*vdot(eng:facing:vector,ship:facing:vector).
      set tff to tff + tmax / (eng:ispat(p) * Kerbing0).
    }
  }
  set vthrust[0] to tt.
  set vthrust[1] to tt/tff.
  print "Max thrust: " + tt at (0,20).
}
  
//PID parameters
local Kp to 0.2.
local Ki to 0.005.
local Kd to 0.02.

//control parameters
local wantedthrust to 0.95.
local sds to -1. //safe descent speed
local safeheight to 4.
local drag to 0.
global muhere to body:mu.

local mode to 1.
on brakes set mode to 2.

lock steering to up.
lock throttle to 0.
wait until status <> "prelaunch".
stage.
wait until mode > 1.

global slp to body:atm:sealevelpressure.
local verticalthrust to list(0,0).
updvthrust(verticalthrust,ship:sensors:pres).
global vthrustsl to list(0,0).
updvthrust(vthrustsl, slp).
global vexsl to vthrustsl[1].
local vv_old to ship:verticalspeed.
local geff to muhere / ship:orbit:position:sqrmagnitude - ship:groundspeed * ship:groundspeed / ship:orbit:position:mag.
local t_old to time:seconds.
local h to alt:radar.
local hbrake to 0.
wait 0. // wait until next frame.

until mode > 2 {
  set tm to time:seconds.
  set vv to ship:verticalspeed.
  set h to alt:radar.
  set dt to tm - t_old.
  set geff to muhere / body:position:sqrmagnitude. // - ship:groundspeed * ship:groundspeed / body:position:mag.
  set drag to ((vv - vv_old) / dt + geff)*mass.
  set pressure to ship:sensors:pres * constant:kpatoatm.
  updvthrust(verticalthrust,pressure).
  set hbrake to BrakeDistance(h, vv, verticalthrust, wantedthrust, geff, pressure, drag).
  print "Effective Isp: " + round(verticalthrust[1]) at (0,21).
  print "dt: " + round(dt, 3) at (0,22).
  print "Pressure: " + round(pressure,2) + "    " at (0,23).
  print "Drag: " + round(drag) + "    " at (0,24).
  print "Suicide burn distance: " + round(hbrake) + "    " at (0,25).
  if h < safeheight + hbrake + abs(vv)*0.25 set mode to 3.
  set vv_old to ship:verticalspeed.
  set t_old to time:seconds.
  wait 0.
}

lock steering to srfretrograde.
local thrustpid to pidloop(Kp, Ki, Kd, -1, 1).
set thrustpid:setpoint to hbrake.

until vv > (geff-verticalthrust[0]/mass)/4 {
  set thrustpid:setpoint to hbrake.
  lock throttle to limiter(wantedthrust + thrustpid:update(time:seconds, alt:radar - safeheight), 0.7, 1).
  set tm to time:seconds.
  set vv to ship:verticalspeed.
  set h to alt:radar.
  set dt to tm - t_old.
  set geff to muhere / body:position:sqrmagnitude. // - ship:groundspeed * ship:groundspeed / body:position:mag.
  set drag to ((vv - vv_old) / dt + geff)*mass - verticalthrust[0]*throttle.
  set pressure to ship:sensors:pres * constant:kpatoatm.
  updvthrust(verticalthrust,pressure).
  set hbrake to BrakeDistance(h, vv, verticalthrust, wantedthrust, geff, pressure, drag).
  print "dt: " + round(dt, 3) at (0,22).
  print "Pressure: " + round(pressure,2) + "    " at (0,23).
  print "Drag: " + round(drag) + "    " at (0,24).
  print "Suicide burn distance: " + round(hbrake) + "    " at (0,25).
  teleprint(thrustpid).
  
  if (vxcl(velocity:surface, up:vector):mag < 0.2) or (hbrake < abs(vv*0.5)) { lock steering to up. }
  else lock steering to srfretrograde.
  
  set vv_old to ship:verticalspeed.
  set t_old to time:seconds.
  wait 0.
}
lock throttle to 0.
lock steering to up.
wait until status = "landed".
set ship:control:pilotmainthrottle to 0.
