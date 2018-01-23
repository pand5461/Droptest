// Wanted thrust levels at start and end of burn
local tlevel0 to 0.9.
local tlevel1 to 0.9.

wait 0.5.
lock throttle to 0.
until ship:availablethrust > 0 {
  stage.
  wait 0.2.
}

list engines in el.
wait until verticalspeed < -10.

lock steering to srfretrograde.
local h_ground to latlng(latitude, longitude):terrainheight. // heightof ground directly underneath
local g_ground to body:mu / (body:radius + h_ground)^2. // gravity at surface
local Isp to el[0]:ispat(body:atm:altitudepressure(h_ground)). // specific impulse of all active engines at terrain level
local maxflow to el[0]:maxthrustat(0) / (el[0]:visp * 9.80665). // total fuel flow at max throttle
local t_ground to ship:availablethrustat(body:atm:altitudepressure(h_ground)). // thrust at ground level at max throttle
local fuelflow to (tlevel0 + tlevel1) * 0.5 * maxflow. // average fuel flow during the burn
local init_acc to availablethrust * tlevel0 / mass. // starting thrust acceleration

local h0 to altitude - 5.
local vv0 to verticalspeed.
local vv1 to 0.
local hstop to alt:radar.
local ghere to body:mu / body:position:sqrmagnitude. // gravity at current altitude
local tstop to -vv0 / (init_acc - ghere). // initial estimate of time-to-stop
local final_acc to t_ground * tlevel1 / (mass - fuelflow * tstop) - g_ground. // estimate of acceleration at touchdown (incl. gravity)
local ave_acc to (2 * (init_acc - ghere) + final_acc) / 3. // estimate of average acceleration through burn
local tstop1 to -vv0 / ave_acc. // corrected estimate of time-to-stop
local t0 to 0. // used to estimate drag
local t1 to 0.
local drag_acc to 0.
local j0 to 6 * (vv0 / tstop + final_acc) / tstop. // wanted jerk (time derivative of acceleration) at burn start

until hstop - h_ground < -vv0 * 0.2 {
  set vv0 to verticalspeed.
  set t0 to time:seconds.
  // refine time-to-stop estimate
  until abs(tstop - tstop1) < 0.05 {
    set tstop to tstop1.
    set final_acc to t_ground * tlevel1 / (mass - fuelflow * tstop) - g_ground.
    set j0 to 6 * (vv0 / tstop + final_acc) / tstop.
    set ave_acc to final_acc - j0 * tstop / 6.
    set tstop1 to -vv0 / ave_acc.
  }

  set tstop to 0.5 * (tstop + tstop1).

  local snap to -j0 / tstop. // wanted snap (time derivative of jerk)
  
  set hstop to h0 + (vv0 + ((init_acc + drag_acc - ghere) * 0.5 + (j0 / 6 + snap * tstop / 24) * tstop) * tstop) * tstop. // wicked maths to calculate altitude ASL at stop if burn starts right now
  print "Stopping alt. above ground: " + round(hstop - h_ground) + "        " at (0,26).
  print "Stop time: " + round(tstop, 2) + "        " at (0, 25).
  wait 0.
  set vv1 to verticalspeed.
  set t1 to time:seconds.
  set init_acc to availablethrust * tlevel0 / mass.
  set h0 to altitude - 5.5.
  set ghere to body:mu / body:position:sqrmagnitude.
  set final_acc to t_ground * tlevel1 / (mass - fuelflow * tstop) - g_ground.
  set ave_acc to (2 * (init_acc + drag_acc - ghere) + final_acc) / 3.
  set tstop1 to -vv1 / ave_acc.
}
// exit from loop must be roughly when suicide burn is to start (a bit higher in atmospheres, because drag is ignored)

set vv0 to vv1.
set t0 to t1.
until vv0 > -0.5 {
  wait 0.
  // estimate drag
  set vv1 to verticalspeed.
  if body:atm:exists {
    set t1 to time:seconds.
    set drag_acc to  (vv1 - vv0) / (t1 - t0) + ghere - min(throttle, 1) * availablethrust / mass.
    set t0 to t1.
  }
  set vv0 to vv1.
  set h0 to alt:radar - 5.5.
  set ghere to body:mu / body:position:sqrmagnitude.
  set final_acc to t_ground * tlevel1 / (mass - fuelflow * tstop) - g_ground.
  // another way to calculate stopping time, now based on current altitude and speed
  local discrim to max(0, (vv1 / 4)^2 + final_acc * h0).
  set tstop1 to max( 0, 2 / final_acc * (vv1 / 4 + sqrt(discrim))).
  until abs(tstop - tstop1) < 0.05 {
    set tstop to tstop1.
    set final_acc to t_ground * tlevel1 / (mass - fuelflow * tstop) - g_ground.
    set discrim to max(0, (vv1 / 4)^2 + final_acc * h0).
    set tstop1 to max( 0, 2 / final_acc * (vv1 / 4 + sqrt(discrim))).
  }
  set tstop to 0.5 * (tstop + tstop1).
  print "Stop time: " + round(tstop, 2) + "        " at (0, 25).
  
  local j0 to 6 * (vv1 / tstop + final_acc) / tstop.
  // after calculating jerk, we can calculate the acceleration needed to stop at ground level, assuming constant snap
  set init_acc to ghere - drag_acc + final_acc - j0 * tstop * 0.5.
  set tlevel to init_acc * mass / availablethrust.
  if tlevel >= tlevel0 or throttle > 0 lock throttle to min(1,tlevel).
  if throttle > 0 set fuelflow to (throttle + tlevel1) * 0.5 * maxflow.
}
lock throttle to 0.
lock steering to up.

wait until ship:status = "landed" or ship:status = "splashed".
set ship:control:pilotmainthrottle to 0.
unlock all.
