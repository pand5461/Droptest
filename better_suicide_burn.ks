clearscreen.
core:doevent("open terminal").
deletepath("0:/log.txt").

local tlevel0 to 0.95.
local tlevel1 to 0.95.
local tlevel to 0.
local minthrott to 0.8.

local g0 to 9.80665.

function ThrustIsp {
  list engines in el.
  local vex to 1.
  local ff to 0.
  local tt to 0.
  local te to 0.
  for e in el {
    set te to e:availablethrust.
    set ff to ff + te/max(e:visp,0.01).
    set tt to tt + te*vdot(facing:vector,e:facing:vector).
  }
  if tt<>0 set vex to g0*tt/ff.
  return list(tt, vex).
}

function ThrustIspat {
  parameter p.
  list engines in el.
  local vex to 1.
  local ff to 0.
  local tt to 0.
  local te to 0.
  for e in el {
    set te to e:availablethrustat(p).
    set ff to ff + te/max(e:ispat(p),0.01).
    set tt to tt + te*vdot(facing:vector,e:facing:vector).
  }
  if tt<>0 set vex to g0*tt/ff.
  return list(tt, vex).
}

wait 0.5.
lock throttle to 0.
until ship:availablethrust > 0 {
  stage.
  wait 0.2.
}

wait until verticalspeed < -10.

lock steering to srfretrograde.
local h0 to alt:radar + vdot(ship:partsnamed("landingLeg1-2")[0]:position - ship:controlpart:position, up:vector) - 1.8.
local h_ground to altitude - h0.
local g_ground to body:mu / (body:radius + h_ground)^2.
local tvex_ground to ThrustIspat(body:atm:altitudepressure(h_ground)).
local t_ground to tvex_ground[0].
local maxflow to t_ground / tvex_ground[1].
set t_ground to tlevel1 * t_ground.
local fuelflow to (tlevel0 + tlevel1) * 0.5 * maxflow.
local ff_ground to tlevel1 * maxflow.

local vv0 to verticalspeed.
local vv1 to 0.
local hstop to 0.
local ghere to body:mu / body:position:sqrmagnitude.
local init_acc to ThrustIsp()[0] * tlevel0 / mass - ghere.
local tstop to -vv0 / init_acc.
local final_mass to mass - fuelflow * tstop.
local final_acc to t_ground / final_mass - g_ground.
local ave_acc to (2 * final_acc + init_acc) / 3.
local tstop1 to -vv0 / ave_acc.
local t0 to 0.
local t1 to 0.
local drag_acc to 0.
local j_ground to 0.

until h0 - hstop < -10 * vv0 {
  set vv0 to verticalspeed.
  set t0 to time:seconds.
  
  local niter to 0.
  until abs(tstop - tstop1) < 0.05 or niter >= 10 {
    set tstop to tstop1.
    set final_mass to mass - fuelflow * tstop.
    set final_acc to t_ground / final_mass - g_ground.
    set j_ground to ff_ground / final_mass * final_acc.
    set tstop1 to ((2 * final_acc + init_acc) - sqrt((2 * final_acc + init_acc)^2 + 6 * vv0 * j_ground)) / j_ground.
    set niter to niter + 1.
  }

  set tstop to 0.5 * (tstop + tstop1).

  set hstop to ((init_acc - final_acc) * tstop / 12 - vv0 * 0.5) * tstop.
  print "Stopping alt. above ground: " + round(h0 - hstop) + "        " at (0,26).
  print "Stop time: " + round(tstop, 2) + "        " at (0, 25).
  wait 0.
  set vv1 to verticalspeed.
  set t1 to time:seconds.
  if body:atm:exists {
    set drag_acc to (vv1 - vv0) / (t1 - t0) + ghere.
  }
  set h0 to alt:radar + vdot(ship:partsnamed("landingLeg1-2")[0]:position - ship:controlpart:position, up:vector) - 1.8.
  set ghere to body:mu / body:position:sqrmagnitude.
  set init_acc to ThrustIsp()[0] * tlevel0 / mass - ghere.
  set final_mass to mass - fuelflow * tstop.
  set final_acc to t_ground / final_mass.
  set j_ground to ff_ground / final_mass * final_acc.
  set final_acc to final_acc - g_ground.
  set tstop1 to ((2 * final_acc + init_acc) - sqrt((2 * final_acc + init_acc)^2 + 6 * vv1 * j_ground)) / j_ground.
}

set vv0 to vv1.
set t0 to t1.
clearscreen.
local t_now to ThrustIsp()[0].
local min_tlevel to tlevel1.
until tstop <= 0.1 {
  wait 0.
  set vv1 to verticalspeed.
  set t1 to time:seconds.
  local dt to t1 - t0.
  if body:atm:exists {
    local expt to 0.5 * body:atm:altitudepressure(altitude) * (body:atm:altitudepressure(altitude - h0) - body:atm:altitudepressure(altitude)) / constant:e.
    print "Drag exponent: " + expt at (0, 20).
    set drag_acc to (t_ground / t_now)^expt * min(ghere, (0.9 * drag_acc + max(0.1 * ((vv1 - vv0) / (t1 - t0) + ghere - min(throttle, 1) * t_now / mass), 0))).
  }
  set t_now to ThrustIsp()[0].
  set vv0 to vv1.
  set h0 to alt:radar + vdot(ship:partsnamed("landingLeg1-2")[0]:position - ship:controlpart:position, up:vector) - 1.8.
  set ghere to body:mu / body:position:sqrmagnitude.
  set h_ground to altitude - alt:radar.
  set g_ground to body:mu / (body:radius + h_ground)^2.
  set t_ground to tlevel1 * ThrustIspat(body:atm:altitudepressure(h_ground))[0].
  local a_total to throttle * t_now / (mass + fuelflow * dt) - g_ground + drag_acc.
  set h0 to h0 + (vv1 + 0.5 * a_total * dt) * dt.
  set vv1 to vv1 + a_total * dt.
  set final_acc to 1.
  set tstop to max(0.01, -2 * h0 / vv1 ).
  set dt to 1.
  if tstop < 1 {
    set dt to 0.
    set init_acc to 0.5 * vv1^2 / h0.
  }
  local niter to 0.
  local jts to 0.
  log tstop + "   " + init_acc to "0:/log.txt".
  local loopstart to time:seconds.
  until abs(dt) < 0.02 or niter >= 8 {
    set final_mass to mass - fuelflow * tstop.
    set final_acc to t_ground / final_mass.
    set j_ground to ff_ground / final_mass * final_acc.
    set final_acc to final_acc - g_ground.
    set jts to j_ground * tstop.
    local f to 4 * h0 + (vv1 + (jts / 6 - final_acc) * tstop) * tstop.
    local fdot to vv1 + tstop * (tstop * (fuelflow / final_mass * jts / 3 + j_ground * (0.5 - fuelflow / ff_ground)) - 2 * final_acc).
    set dt to f / fdot.
    set tstop to tstop - dt.
    set init_acc to 7 * final_acc - jts - 12 * (vv1 + 3 * h0 / tstop) / tstop.
    log tstop + "   " + init_acc to "0:/log.txt".
    set niter to niter + 1.
  }
  log "=================" to "0:/log.txt".
  log "Tstop converged in " + round(time:seconds - loopstart, 2) to "0:/log.txt".
  log "Throttle: " + throttle to "0:/log.txt".

  print "Stop time: " + round(tstop, 2) + "        " at (0, 25).
  print "Stopping alt. above ground: " + round(h0 - 0.25 * tstop * (tstop * (final_acc - jts / 6) - vv1)) + "     " at (0,26).
  
  if abs(dt) < 0.5 and tstop > 0 {
    local j_init to 2 * (final_acc - init_acc) / tstop - j_ground.
    set tlevel to (init_acc + ghere - drag_acc) * mass / t_now.
    local tj_init to mass * sqrt(max(0, j_init) / (maxflow * t_now)).
    log "Req'd throttle : " + tlevel to "0:/log.txt".
    if tj_init > tlevel set tlevel to 0.2 * (4 * tlevel + min(1, tj_init)).
    log "Drag acc: " + drag_acc to "0:/log.txt".
    log "Req'd throttle from j: " + tj_init to "0:/log.txt".
    if throttle > 0 or tlevel > tlevel0 {
      lock throttle to max(minthrott, min(1, tlevel)).
      print "Throttle level: " + round(throttle, 2) + "    " at (0, 23).
    }
  }
  if throttle > 0 set fuelflow to (throttle + tlevel1) * 0.5 * maxflow.
  set t0 to t1.
}
lock throttle to 0.
set ship:control:pilotmainthrottle to 0.
lock steering to up.

wait 1.
unlock throttle. unlock steering.
wait until ship:status = "landed" or ship:status = "splashed".
print "The ship has " + ship:status + ".".
