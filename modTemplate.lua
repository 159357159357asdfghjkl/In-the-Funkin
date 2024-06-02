print('LOADED TEMPLATE')

TEMPLATE = {}

-- Taro's janky "TEMPLATE 1" implementation for FNF
-- shoutouts to Kade for letting me do this
-- psyched by a invisible man
-- it's a 2D modchart system so this makes it don't support zIndex
-- use dofile or addLuaScript to call the modchart template in order to write mods
-- use psych engine 0.7 or higher version
-- DO NOT CHANGE ALL THE FUNCTIONS

-- EASING EQUATIONS

---------------------------------------------------------------------------------------
----------------------DON'T TOUCH IT KIDDO---------------------------------------------
---------------------------------------------------------------------------------------
			
-- Adapted from
-- Tweener's easing functions (Penner's Easing Equations)
-- and http://code.google.com/p/tweener/ (jstweener javascript version)
--

--[[
Disclaimer for Robert Penner's Easing Equations license:

TERMS OF USE - EASING EQUATIONS

Open source under the BSD License.

Copyright © 2001 Robert Penner
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    * Neither the name of the author nor the names of contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
]]

-- For all easing functions:
-- t = elapsed time
-- b = begin
-- c = change == ending - beginning
-- d = duration (total time)

local pow = math.pow
local sin = math.sin
local cos = math.cos
local pi = math.pi
local sqrt = math.sqrt
local abs = math.abs
local asin  = math.asin

function linear(t, b, c, d)
  return c * t / d + b
end

function inQuad(t, b, c, d)
  t = t / d
  return c * pow(t, 2) + b
end

function outQuad(t, b, c, d)
  t = t / d
  return -c * t * (t - 2) + b
end

function inOutQuad(t, b, c, d)
  t = t / d * 2
  if t < 1 then
    return c / 2 * pow(t, 2) + b
  else
    return -c / 2 * ((t - 1) * (t - 3) - 1) + b
  end
end

function outInQuad(t, b, c, d)
  if t < d / 2 then
    return outQuad (t * 2, b, c / 2, d)
  else
    return inQuad((t * 2) - d, b + c / 2, c / 2, d)
  end
end

function inCubic (t, b, c, d)
  t = t / d
  return c * pow(t, 3) + b
end

function outCubic(t, b, c, d)
  t = t / d - 1
  return c * (pow(t, 3) + 1) + b
end

function inOutCubic(t, b, c, d)
  t = t / d * 2
  if t < 1 then
    return c / 2 * t * t * t + b
  else
    t = t - 2
    return c / 2 * (t * t * t + 2) + b
  end
end

function outInCubic(t, b, c, d)
  if t < d / 2 then
    return outCubic(t * 2, b, c / 2, d)
  else
    return inCubic((t * 2) - d, b + c / 2, c / 2, d)
  end
end

function inQuart(t, b, c, d)
  t = t / d
  return c * pow(t, 4) + b
end

function outQuart(t, b, c, d)
  t = t / d - 1
  return -c * (pow(t, 4) - 1) + b
end

function inOutQuart(t, b, c, d)
  t = t / d * 2
  if t < 1 then
    return c / 2 * pow(t, 4) + b
  else
    t = t - 2
    return -c / 2 * (pow(t, 4) - 2) + b
  end
end

function outInQuart(t, b, c, d)
  if t < d / 2 then
    return outQuart(t * 2, b, c / 2, d)
  else
    return inQuart((t * 2) - d, b + c / 2, c / 2, d)
  end
end

function inQuint(t, b, c, d)
  t = t / d
  return c * pow(t, 5) + b
end

function outQuint(t, b, c, d)
  t = t / d - 1
  return c * (pow(t, 5) + 1) + b
end

function inOutQuint(t, b, c, d)
  t = t / d * 2
  if t < 1 then
    return c / 2 * pow(t, 5) + b
  else
    t = t - 2
    return c / 2 * (pow(t, 5) + 2) + b
  end
end

function outInQuint(t, b, c, d)
  if t < d / 2 then
    return outQuint(t * 2, b, c / 2, d)
  else
    return inQuint((t * 2) - d, b + c / 2, c / 2, d)
  end
end

function inSine(t, b, c, d)
  return -c * cos(t / d * (pi / 2)) + c + b
end

function outSine(t, b, c, d)
  return c * sin(t / d * (pi / 2)) + b
end

function inOutSine(t, b, c, d)
  return -c / 2 * (cos(pi * t / d) - 1) + b
end

function outInSine(t, b, c, d)
  if t < d / 2 then
    return outSine(t * 2, b, c / 2, d)
  else
    return inSine((t * 2) -d, b + c / 2, c / 2, d)
  end
end

function inExpo(t, b, c, d)
  if t == 0 then
    return b
  else
    return c * pow(2, 10 * (t / d - 1)) + b - c * 0.001
  end
end

function outExpo(t, b, c, d)
  if t == d then
    return b + c
  else
    return c * 1.001 * (-pow(2, -10 * t / d) + 1) + b
  end
end

function inOutExpo(t, b, c, d)
  if t == 0 then return b end
  if t == d then return b + c end
  t = t / d * 2
  if t < 1 then
    return c / 2 * pow(2, 10 * (t - 1)) + b - c * 0.0005
  else
    t = t - 1
    return c / 2 * 1.0005 * (-pow(2, -10 * t) + 2) + b
  end
end

function outInExpo(t, b, c, d)
  if t < d / 2 then
    return outExpo(t * 2, b, c / 2, d)
  else
    return inExpo((t * 2) - d, b + c / 2, c / 2, d)
  end
end

function inCirc(t, b, c, d)
  t = t / d
  return(-c * (sqrt(1 - pow(t, 2)) - 1) + b)
end

function outCirc(t, b, c, d)
  t = t / d - 1
  return(c * sqrt(1 - pow(t, 2)) + b)
end

function inOutCirc(t, b, c, d)
  t = t / d * 2
  if t < 1 then
    return -c / 2 * (sqrt(1 - t * t) - 1) + b
  else
    t = t - 2
    return c / 2 * (sqrt(1 - t * t) + 1) + b
  end
end

function outInCirc(t, b, c, d)
  if t < d / 2 then
    return outCirc(t * 2, b, c / 2, d)
  else
    return inCirc((t * 2) - d, b + c / 2, c / 2, d)
  end
end

function inElastic(t, b, c, d, a, p)
  if t == 0 then return b end

  t = t / d

  if t == 1  then return b + c end

  if not p then p = d * 0.3 end

  local s

  if not a or a < abs(c) then
    a = c
    s = p / 4
  else
    s = p / (2 * pi) * asin(c/a)
  end

  t = t - 1

  return -(a * pow(2, 10 * t) * sin((t * d - s) * (2 * pi) / p)) + b
end

-- a: amplitud
-- p: period
function outElastic(t, b, c, d, a, p)
  if t == 0 then return b end

  t = t / d

  if t == 1 then return b + c end

  if not p then p = d * 0.3 end

  local s

  if not a or a < abs(c) then
    a = c
    s = p / 4
  else
    s = p / (2 * pi) * asin(c/a)
  end

  return a * pow(2, -10 * t) * sin((t * d - s) * (2 * pi) / p) + c + b
end

-- p = period
-- a = amplitud
function inOutElastic(t, b, c, d, a, p)
  if t == 0 then return b end

  t = t / d * 2

  if t == 2 then return b + c end

  if not p then p = d * (0.3 * 1.5) end
  if not a then a = 0 end

  local s

  if not a or a < abs(c) then
    a = c
    s = p / 4
  else
    s = p / (2 * pi) * asin(c / a)
  end

  if t < 1 then
    t = t - 1
    return -0.5 * (a * pow(2, 10 * t) * sin((t * d - s) * (2 * pi) / p)) + b
  else
    t = t - 1
    return a * pow(2, -10 * t) * sin((t * d - s) * (2 * pi) / p ) * 0.5 + c + b
  end
end

-- a: amplitud
-- p: period
function outInElastic(t, b, c, d, a, p)
  if t < d / 2 then
    return outElastic(t * 2, b, c / 2, d, a, p)
  else
    return inElastic((t * 2) - d, b + c / 2, c / 2, d, a, p)
  end
end

function inBack(t, b, c, d, s)
  if not s then s = 1.70158 end
  t = t / d
  return c * t * t * ((s + 1) * t - s) + b
end

function outBack(t, b, c, d, s)
  if not s then s = 1.70158 end
  t = t / d - 1
  return c * (t * t * ((s + 1) * t + s) + 1) + b
end

function inOutBack(t, b, c, d, s)
  if not s then s = 1.70158 end
  s = s * 1.525
  t = t / d * 2
  if t < 1 then
    return c / 2 * (t * t * ((s + 1) * t - s)) + b
  else
    t = t - 2
    return c / 2 * (t * t * ((s + 1) * t + s) + 2) + b
  end
end

function outInBack(t, b, c, d, s)
  if t < d / 2 then
    return outBack(t * 2, b, c / 2, d, s)
  else
    return inBack((t * 2) - d, b + c / 2, c / 2, d, s)
  end
end

function outBounce(t, b, c, d)
  t = t / d
  if t < 1 / 2.75 then
    return c * (7.5625 * t * t) + b
  elseif t < 2 / 2.75 then
    t = t - (1.5 / 2.75)
    return c * (7.5625 * t * t + 0.75) + b
  elseif t < 2.5 / 2.75 then
    t = t - (2.25 / 2.75)
    return c * (7.5625 * t * t + 0.9375) + b
  else
    t = t - (2.625 / 2.75)
    return c * (7.5625 * t * t + 0.984375) + b
  end
end

function inBounce(t, b, c, d)
  return c - outBounce(d - t, 0, c, d) + b
end

function inOutBounce(t, b, c, d)
  if t < d / 2 then
    return inBounce(t * 2, 0, c, d) * 0.5 + b
  else
    return outBounce(t * 2 - d, 0, c, d) * 0.5 + c * .5 + b
  end
end

function outInBounce(t, b, c, d)
  if t < d / 2 then
    return outBounce(t * 2, b, c / 2, d)
  else
    return inBounce((t * 2) - d, b + c / 2, c / 2, d)
  end
end

function scale(x, l1, h1, l2, h2)
	return (((x) - (l1)) * ((h2) - (l2)) / ((h1) - (l1)) + (l2))
end

function math.clamp(val,min,max)
	if val < min then return min end
	if val > max then return max end
	return val
end

---------------------------------------------------------------------------------------
----------------------END DON'T TOUCH IT KIDDO-----------------------------------------
---------------------------------------------------------------------------------------

beat = 0
ARROW_SIZE = 112

--all of our mods, with default values
modList = {
	beat = 0,
	flip = 0,
	invert = 0,
	drunk = 0,
	tipsy = 0,
	adrunk = 0, --non conflict accent mod
	atipsy = 0, --non conflict accent mod
	movex = 0,
	movey = 0,
	amovex = 0,
	amovey = 0,
	reverse = 0,
	split = 0,
	cross = 0,
	dark = 0,
	stealth = 0,
	alpha = 1,
	confusion = 0,
	dizzy = 0,
	wave = 0,
	brake = 0,
	hidden = 0,
	hiddenoffset = 0,
	alternate = 0,
	camx = 0,
	camy = 0,
	rotationz = 0,
	camwag = 0,
	xmod = 1, --scrollSpeed
    cmod = -1,
	drawsize = 10 --beatcutoff
}

--column specific mods
for i=0,3 do
	modList['movex'..i] = 0
	modList['movey'..i] = 0
	modList['amovex'..i] = 0
	modList['amovey'..i] = 0
	modList['dark'..i] = 0
	modList['stealth'..i] = 0
	modList['confusion'..i] = 0
	modList['reverse'..i] = 0
	modList['xmod'..i] = 1 --column specific scrollSpeed multiplier
    modList['cmod'..i] = -1
end

activeMods = {{},{}}

for pn=1,2 do
	for k,v in pairs(modList) do
		activeMods[pn][k] = v
	end
end

storedMods = {{},{}}
targetMods = {{},{}}
isTweening = {{},{}}
tweenStart = {{},{}}
tweenLen = {{},{}}
tweenCurve = {{},{}}
tweenEx1 = {{},{}}
tweenEx2 = {{},{}}
modnames = {}

function definemod(t)
	local k,v = t[1],t[2]
	if not v then v = 0 end
	for pn=1,2 do
		storedMods[pn][k] = v
		targetMods[pn][k] = v
		isTweening[pn][k] = false
		tweenStart[pn][k] = 0
		tweenLen[pn][k] = 0
		tweenCurve[pn][k] = linear
		tweenEx1[pn][k] = nil
		tweenEx2[pn][k] = nil
		if pn == 1 then
			--print('registered modifier: '..k)
			table.insert(modnames,k)
		end
	end
end

function TEMPLATE.InitMods()
	for pn=1,2 do
		for k,v in pairs(activeMods[pn]) do
			definemod{k,v}
		end
	end
end

function TEMPLATE.setup()
	--sort tables, optimization step
	function modtable_compare(a,b)
		return a[1] < b[1]
	end
	
	if table.getn(event) > 1 then
		table.sort(event, modtable_compare)
	end
	if table.getn(mods) > 1 then
		table.sort(mods, modtable_compare)
	end
end



function receptorAlpha(iCol,pn)
	local alp = 1
	
	local m = activeMods[pn]
	
	if m.alpha ~= 1 then
		alp = alp*m.alpha
	end
	if m.dark ~= 0 or m['dark'..iCol] ~= 0 then
		alp = alp*(1-m.dark)*(1-m['dark'..iCol])
	end
	
	return alp
end

function arrowAlpha(fYOffset, iCol,pn)
	local alp = 1
	
	local m = activeMods[pn]
	
	if m.alpha ~= 1 then
		alp = alp*m.alpha
	end
	if m.stealth ~= 0 or m['stealth'..iCol] ~= 0 then
		alp = alp*(1-m.stealth)*(1-m['stealth'..iCol])
	end
	if m.hidden ~= 0 then
		if fYOffset < m.hiddenoffset and fYOffset >= m.hiddenoffset-200 then
			local hmult = -(fYOffset-m.hiddenoffset)/200
			alp = alp*(1-hmult)*m.hidden
		elseif fYOffset < m.hiddenoffset-100 then
			alp = alp*(1-m.hidden)
		end
	end
	return alp
end

function getReverseForCol(iCol,pn)
	local m = activeMods[pn]
	local val = 0
	
	val = val+m.reverse+m['reverse'..iCol]
	
	if m.split ~= 0 and iCol > 1 then val = val+m.split end
	if m.cross ~= 0 and iCol == 1 or iCol == 2 then val = val+m.cross end
	if m.alternate ~= 0 and iCol % 2 == 1 then val = val+m.alternate end
	
	return val
end

function getYAdjust(fYOffset, iCol, pn)
	
	local m = activeMods[pn]
	
	local yadj = 0
	if m.wave ~= 0 then
		yadj = yadj + m.wave * 20*math.sin( (fYOffset+250)/76 )
	end
	
	if m.brake ~= 0 then

		local fEffectHeight = 500;
		local fScale = scale( fYOffset, 0, fEffectHeight, 0, 1 )
		local fNewYOffset = fYOffset * fScale; 
		local fBrakeYAdjust = m.brake * (fNewYOffset - fYOffset)
		
		fBrakeYAdjust = math.clamp( fBrakeYAdjust, -400, 400 )
		yadj = yadj+fBrakeYAdjust;
	
	end
	
	fYOffset = fYOffset+yadj
	
	return fYOffset
end

function arrowEffects(fYOffset, iCol, pn)
    local m = activeMods[pn]
	
    local xpos, ypos, rotz = 0, 0, 0
	
	if m['confusion'..iCol] ~= 0 or m.confusion ~= 0 then
		rotz = rotz + m['confusion'..iCol] + m.confusion
	end
	if m.dizzy ~= 0 then
		rotz = rotz + m.dizzy*fYOffset
	end
	
    if m.drunk ~= 0 then
        xpos = xpos + m.drunk * ( math.cos( getSongPosition()*0.001 + iCol*(0.2) + 1*(0.2) + fYOffset*(10)/720) * ARROW_SIZE*0.5 )
    end
    if m.tipsy ~= 0 then
        ypos = ypos + m.tipsy * ( math.cos( getSongPosition()*0.001 *(1.2) + iCol*(2.0) + 1*(0.2) ) * ARROW_SIZE*0.4 )
    end
    if m.adrunk ~= 0 then
        xpos = xpos + m.adrunk * ( math.cos( getSongPosition()*0.001 + iCol*(0.2) + 1*(0.2) + fYOffset*(10)/720) * ARROW_SIZE*0.5 )
    end
    if m.atipsy ~= 0 then
        ypos = ypos + m.atipsy * ( math.cos( getSongPosition()*0.001 *(1.2) + iCol*(2.0) + 1*(0.2) ) * ARROW_SIZE*0.4 )
    end
	
	if m['movex'..iCol] ~= 0 or m.movex ~= 0 then
		xpos = xpos + m['movex'..iCol] + m.movex
	end
	if m['amovex'..iCol] ~= 0 or m.amovex ~= 0 then
		xpos = xpos + m['amovex'..iCol] + m.amovex
	end
	if m['movey'..iCol] ~= 0 or m.movey ~= 0 then
		ypos = ypos + m['movey'..iCol] + m.movey
	end
	if m['amovey'..iCol] ~= 0 or m.amovey ~= 0 then
		ypos = ypos + m['amovey'..iCol] + m.amovey
	end
	
	if m['reverse'..iCol] ~= 0 or m.reverse ~= 0 or m.split ~= 0 or m.cross ~= 0 or m.alternate ~= 0 then
		ypos = ypos + getReverseForCol(iCol,pn) * 450
	end
	
	if m.flip ~= 0 then
		local fDistance = ARROW_SIZE * 2 * (1.5 - iCol);
		xpos = xpos + fDistance * m.flip;
	end

	if m.invert ~= 0 then
		local fDistance = ARROW_SIZE * (iCol%2 == 0 and 1 or -1);
		xpos = xpos + fDistance * m.invert;
	end
	
	if m.beat ~= 0 then
			
		local fBeatStrength = m.beat;
		
		local fAccelTime = 0.3;
		local fTotalTime = 0.7;
		
		-- If the song is really fast, slow down the rate, but speed up the
		-- acceleration to compensate or it'll look weird.
		fBeat = beat + fAccelTime;
		
		local bEvenBeat = false;
		if math.floor(fBeat) % 2 ~= 0 then
			bEvenBeat = true;
		end
		
		fBeat = fBeat-math.floor( fBeat );
		fBeat = fBeat+1;
		fBeat = fBeat-math.floor( fBeat );
		
		if fBeat<fTotalTime then
		
			local fAmount = 0;
			if fBeat < fAccelTime then
				fAmount = scale( fBeat, 0.0, fAccelTime, 0.0, 1.0);
				fAmount = fAmount*fAmount;
			else 
				--fBeat < fTotalTime
				fAmount = scale( fBeat, fAccelTime, fTotalTime, 1.0, 0.0);
				fAmount = 1 - (1-fAmount) * (1-fAmount);
			end

			if bEvenBeat then
				fAmount = fAmount*-1;
			end

			local fShift = 40.0*fAmount*math.sin( ((fYOffset/30.0)) + (math.pi/2) );
			
			xpos = xpos + fBeatStrength * fShift
			
		end
	
	end
	
    return xpos, ypos, rotz
    
end




defaultPositions = {{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0}}



-- events

mods,curmod = {},1
perframe = {}
event,curevent = {},1
songStarted = false

function ease(t)
	table.insert(mods,t)
end

function set(t)
	table.insert(mods,{t[1],0,linear,t[2],t[3],pn=t.pn})
end

function callback(t)
	table.insert(event,t)
end


function beats(t)
	table.insert(perframe,t)
end

function me(t)
    --deprecated function, aka ease
	table.insert(mods,t)
end
function mpf(t)
    --deprecated aka beats
	table.insert(perframe,t)
end
function m2(t)
	--deprecated aka callback
	table.insert(event,t)
end

function TEMPLATE.songStart()
    
    downscroll = false

	for i=0,7 do
        defaultPositions[i+1].x = getPropertyFromGroup("strumLineNotes",i,"x")+(i < 4 and 20 or 0)
        defaultPositions[i+1].y = getPropertyFromGroup("strumLineNotes",i,"y")

        --print(i .. ": " .. defaultPositions[i+1].x .. " " .. defaultPositions[i+1].y)
    end
	
	--fuck it, it's mods. You don't get a say here.
	--(this is done to prevent a lot of bugs and weird cases)
	storedScrollSpeed = 1.8
	--storedScrollSpeed = scrollSpeed
	
	for pn=1,2 do
		activeMods[pn].xmod = storedScrollSpeed
	end
	
	songStarted = true
	
end

function TEMPLATE.update(elapsed)
    beat = (getSongPosition() / 1000) * (curBpm/60)
    beatcutoff = activeMods[1].drawsize
	
	--------------------------------------------------------------
	-- modified version of exschwasion's template 1 ease reader
	-- format changed to make it more mirin-like
	-- v[1]=startBeat v[2]=len/end v[3]=curve v[4]=newval v[5]=modname, v.pn=player
	-- len is now implied, but v.timing='len' or v.timing='end' for specifics
	-- v.startVal for specifying new start val
	-- v.ex1, v.ex2 for the extra params of elastic and back eases
	--------------------------------------------------------------
	
	while curmod <= #mods and beat > mods[curmod][1] do
		local v = mods[curmod]
		
		local mn = v[5]
		local dur = v[2]
		if v.timing and v.timing == 'end' then
			dur = v[2]-v[1]
		end
		
		--print('launch attack '..mn..' at beat '..v[1])
		
		if v.plr and not v.pn then v.pn = v.plr end
		
		for pn=1,2 do
			if not v.pn or pn == v.pn then
				tweenStart[pn][mn] = v[1]
				tweenLen[pn][mn] = dur
				tweenCurve[pn][mn] = v[3]
				if v.startVal then
					storedMods[pn][mn] = v.startVal
				else
					storedMods[pn][mn] = activeMods[pn][mn]
				end
				targetMods[pn][mn] = v[4]
				tweenEx1[pn][mn] = v.ex1
				tweenEx2[pn][mn] = v.ex2
				isTweening[pn][mn] = true
			end
		end
		
		curmod = curmod+1
	end
	
	for pn=1,2 do
		for _,v in pairs(modnames) do
			
			if isTweening[pn][v] then
				local curtime = beat - tweenStart[pn][v]
				local duration = tweenLen[pn][v]
				local startstrength = storedMods[pn][v]
				local diff = targetMods[pn][v] - startstrength
				local curve = tweenCurve[pn][v]
				local strength = curve(curtime, startstrength, diff, duration, tweenEx1[pn][v], tweenEx2[pn][v])
				activeMods[pn][v] = strength
				if beat > tweenStart[pn][v]+duration then
					isTweening[pn][v] = false
				end
			else
				activeMods[pn][v] = targetMods[pn][v]
			end
			
		end
	end
	
	----------------------------------------
	-- do this stuff every frame --
	----------------------------------------
	if #perframe>0 then
		for i=1,#perframe do
			local a = perframe[i]
			if beat > a[1] and beat < a[2] then
				a[3](beat)
			end
		end
	end
	
	-----------------------------------------
	-- event queue --event,curevent = {},1 --
	-----------------------------------------
	while curevent <= #event and beat>=event[curevent][1] do
		if event[curevent][3] or beat < event[curevent][1]+2 then
			event[curevent][2]()
		end
		curevent = curevent+1;
	end
	
	---------------------------------------
	-- ACTUALLY APPLY THE RESULTS OF THE ABOVE CALCULATIONS TO THE NOTES
	---------------------------------------

	setCamNotes(activeMods[1].camx,activeMods[1].camy,activeMods[1].rotationz + activeMods[1].camwag * math.sin(beat*math.pi))
	
	if songStarted then
		for pn=1,2 do
			local xmod = activeMods[pn].xmod
			for col=0,3 do
			
				local c = (pn-1)*4 + col
				local xp, yp, rz = arrowEffects(0, col, pn)
				local alp = receptorAlpha(col,pn)

				--print('Areceptor '..c..' is '..tostring(receptor))
			
				local defaultx, defaulty = defaultPositions[c+1].x, defaultPositions[c+1].y
				setPropertyFromGroup("strumLineNotes",c,"x",defaultx + xp)
                setPropertyFromGroup("strumLineNotes",c,"y",defaulty + yp)
                setPropertyFromGroup("strumLineNotes",c,"angle",rz)
                setPropertyFromGroup("strumLineNotes",c,"alpha",alp)

				
				--local scrollSpeed = xmod * activeMods[pn]['xmod'..col] * (1 - 2*getReverseForCol(col,pn))
				--setLaneScrollspeed(c,scrollSpeed)
				
				--print('Breceptor '..c..' is '..tostring(receptor))
				
			end
		end
		
		--for i=1,getNumberOfNotes() do
		--	local note = _G['note_'..i]

		for v = 0, getProperty("notes.length")-1 do
			if getPropertyFromGroup('notes',v,"alive") then
				
				--print(tostring(note)..' sus '..tostring(note.isSustain))
				
				local pn = 1
				if getPropertyFromGroup('notes',v,"mustPress") then pn = 2 end
				
				
				local xmod = activeMods[pn].xmod
				
				local isSus = getPropertyFromGroup('notes',v,"isSustainNote")
				--local isParent = note.isParent
				local col = getPropertyFromGroup('notes',v,"noteData")
				local c = (pn-1)*4 + col
				
				local targTime = getPropertyFromGroup('notes',v,"strumTime")
				
				local defaultx, defaulty = defaultPositions[c+1].x, defaultPositions[c+1].y

				local function getcmod()
				local daSpeed = activeMods[pn]['xmod'..col]
				if daSpeed < 0 then
					daSpeed = activeMods[pn].cmod < 0 and getProperty('SONG.speed') or activeMods[pn].cmod
					return daSpeed
				end
					return daSpeed;
				end
				local scrollSpeed = xmod * activeMods[pn]['xmod'..col] * getcmod() * getPropertyFromGroup("notes",v,"multSpeed") * (1 - 2*getReverseForCol(col,pn))
				
				local off = (1 - 2*getReverseForCol(col,pn))

				local ypos = getYAdjust(defaulty - (getSongPosition() - targTime),col,pn) * scrollSpeed * 0.45 - off
				
				local xa, ya, rz = arrowEffects(ypos-defaulty, col, pn)
				local alp = arrowAlpha(ypos-defaulty, col, pn)
				
				if getPropertyFromGroup('notes',v,"isSustainNote") --[[and not note.isParent]] then
					local ypos2 = getYAdjust(defaulty - ((getSongPosition()+.1) - targTime),col,pn) * scrollSpeed * 0.45 - off
					local xa2, ya2 = arrowEffects(ypos2-defaulty, col, pn)

					--if scrollSpeed >= 0 then
					setPropertyFromGroup("notes",v,"angle",math.deg(math.atan2(((ypos2 + ya2)-(ypos + ya))*100,(xa2-xa)*100) + math.pi/2))
					--else
					--	note.angle = 180+math.deg(math.atan2(((ypos2 + ya2)-(ypos + ya))*100,(xa2-xa)*100) + math.pi/2)
					--end
				else
					setPropertyFromGroup("notes",v,"angle",rz)
				end
            	setPropertyFromGroup("notes",v,"x",defaultx + xa + (getPropertyFromGroup('notes',v,"isSustainNote") and 35 or 0))
            	setPropertyFromGroup("notes",v,"y",ypos + ya)
            	setPropertyFromGroup("notes",v,"alpha",alp)
				
				
			end
			
		end
		
	end
	
	
end


function onCreatePost()
	TEMPLATE.InitMods()

	--WRITE MODS HERE! 
	

	run()
	
	--must be called at END of start
	TEMPLATE.setup()
	
end
function run()
	if songName == 'detected' then --example
	-------START-------
	local disable = false
	
	local ALLOW_REVERSE = true

	function flicker()
	if (getSongPosition() * 0.001 * 60) % 2 < 1 then
		return -1
	else
		return 1
	end
	end
	
	--some basic helpers
	function hide(t)
		local bt,tpn = t[1],t.pn
		for i=0,3 do
			me{bt+i*.125-1,.5,outExpo,-70,'movey'..i,pn=tpn}
			me{bt+i*.125-.5,1.25,inExpo,650,'movey'..i,pn=tpn}
			set{bt+i*.125+1.75,1,'stealth',pn=tpn}
			set{bt+i*.125+1.75,1,'dark',pn=tpn}
		end
	end
	function unhide(t)
		local bt,tpn = t[1],t.pn
		for i=0,3 do
			set{bt+i*.125-2,0,'stealth',pn=tpn}
			set{bt+i*.125-2,0,'dark',pn=tpn}
			me{bt+i*.125-2,1,outExpo,-70,'movey'..i,pn=tpn}
			me{bt+i*.125-1,1,inExpo,50,'movey'..i,pn=tpn}
			me{bt+i*.125-0,1.25,outElastic,0,'movey'..i,pn=tpn}
		end
	end
	--wiggle(beat,num,div,ease,amt,mod)
	function wig(t)
		local b,num,div,ea,am,mo = t[1],t[2],t[3],t[4],t[5],t[6]
		local f = 1
		for i=0,num do
			local smul = i==0 and 1 or 0
			local emul = i==num and 0 or 1
			
			me{b+i*(1/div),1/div,ea,startVal = am*smul*f, am*emul*-f,mo,pn=t.pn}
			
			f = f*-1
		end
	end
	--simple mod 2
	local function sm2(tab)
		local b,len,eas,amt,mods,intime = tab[1],tab[2],tab[3],tab[4],tab[5],tab.intime
		if not intime then intime = .1 end
		if intime <= 0 then intime = .001 end
		me{b-intime,intime,linear,amt,mods,pn=tab.pn}
		me{b,len-intime,eas,0,mods,pn=tab.pn}
	end
	
		definemod{'colspacing',112}
		definemod{'spacing',620}
		definemod{'addx',0}
		definemod{'waveamp',0}
		
		definemod{'halo',0}
		
	
	
	
	
	
	set{0,1.8,'xmod'}
	set{0,32,'drawsize'}
	
	if not disable then
		me{0,2,inCirc,1,'dizzy'}
		me{0,2,inOutQuad,1,'camwag'}
		me{0,4,outCubic,2,'tipsy',pn=1}
		me{12,4,inCubic,1,'tipsy',pn=1}
		
		me{16,4,outCubic,1,'tipsy',pn=2}
		me{28,4,inCubic,0,'tipsy'}
		
		
		
		--kade "sway"
		mpf{32,64,function(beat)
			
			local pn = 1
			
			if beat > 48 then
				pn = 2
			end
			
			for col = 0,3 do
				activeMods[pn]['movex'..col] = 32 * math.sin((beat + col*0.25) * math.pi)
				local mu = col%2 == 0 and 1 or -1
				activeMods[pn]['movey'..col] = 32 * mu * math.tan((beat) * 0.25 * math.pi)
			end
			
		end}
		
		local f = 1
		for i=32,63,4 do
			--spinny
			local spn = 1
			if i >= 48 then
				spn = 2
			end
			me{i,2,outCubic,startVal=360*f,0,'confusion',pn=spn}
			f = f*-1
		end
		
		
		
		me{62,4,inOutCubic,.4,'alpha',pn=1}
		me{94,4,inOutCubic,1,'alpha',pn=1}
		
		set{91.5,3,'beat'}
		set{95.5,0,'beat'}
		
		
		--NEW wrapping
		--WITHOUT FLICKER
		me{64,8,inOutCubic,1280/8,'colspacing'}
		me{64,8,inOutCubic,(1280/8)*4,'spacing'}
		me{64,8,inOutCubic,-176,'addx'}
		me{64,8,inOutCubic,1,'waveamp'}
		
		me{92,4,inCubic,112,'colspacing'}
		me{92,4,inCubic,620,'spacing'}
		me{92,4,inCubic,0,'addx'}
		me{92,4,inCubic,0,'waveamp'}
		
		mpf{64,95,function(beat)
			
			local scrollerpos = beat-64
			
			for pn=1,2 do
			
				for col = 0,3 do
					local cpos = col*-112
					if pn == 2 then cpos = cpos-620 end
					local c = (pn-1)*4+col
					
					local colspacing = activeMods[pn].colspacing
					local spacing = activeMods[pn].spacing
					local addx = activeMods[pn].addx
					local waveamp = activeMods[pn].waveamp
					
					local newpos = (col*colspacing + (pn-1)*spacing + scrollerpos*1280/8)%1280 + addx
					
					activeMods[pn]['movex'..col] = cpos + newpos
				
					local ang = 2*math.pi*(c/8)
					activeMods[pn]['reverse'..col] = waveamp * .1 * math.sin(beat*math.pi + ang)
					
				end
				
			end
			
		end}
		
		
		local wso = 0
		
		hide{(96-wso)-1,pn=2}
		unhide{(96-wso)+16,pn=2}
		hide{(96-wso)+16,pn=1}
		unhide{(96-wso)+32,pn=1}
		hide{(96-wso)+32,pn=2}
		unhide{(96-wso)+48,pn=2}
		hide{(96-wso)+48,pn=1}
		unhide{(96-wso)+64,pn=1}
		
		me{(96-wso)-3,3,linear,1.6,'xmod'}
		--me{(96-wso)-3,3,linear,1,'drunk'}
		--set{(96-wso)-3,1,'drunk'}
		
		me{(96-wso)-4,3,linear,50,'movey'}
		me{(96-wso),3,linear,-.1,'flip'}
		
		me{(96-wso)+28,4,inCubic,0,'flip'}
		
		for i=0,1 do
			me{(96-wso)+8*i,4,linear,448*1.2,'movex',pn=1}
			me{(96-wso)+4+8*i,4,linear,0,'movex',pn=1}
			me{(96-wso)+16+8*i,4,linear,-448*1.2,'movex',pn=2}
			me{(96-wso)+16+4+8*i,4,linear,0,'movex',pn=2}
		end
		
		mpf{(96-wso),(96-wso)+64,function(beat)
			
			local pn=1
			if (beat-(96-wso)) % 32 > 16 then pn = 2 end
			
			local pingpong = beat - math.floor(beat)
			if beat%2 > 1 then
				pingpong = 1-pingpong
			end
			
			local wdir = -1 * (pn*2-3)
			if (beat-(96-wso)) % 8 > 4 then wdir = wdir*-1 end
			
			if beat < (96-wso)+32 then
				for col=0,3 do
					activeMods[pn]['amovey'..col] = math.min(-112 * wdir * math.sin(beat*math.pi + col*math.pi),0)
				end
				activeMods[pn].invert = pingpong*1.2
			else
				for col=0,3 do
					activeMods[pn]['amovey'..col] = -56*1.2 * wdir * math.sin(beat*math.pi + col*math.pi)
				end
				activeMods[pn].invert = 0.6 - 0.6*math.cos(beat*math.pi)
			end
			
		end}
		
		
		local kick = {
			{208.000,0,1},
			{208.750,0,1},
			{209.500,0,1},
			{210.250,0,1},
			{211.000,3,1},
			{211.500,3,1},
			{212.000,0,1},
			{212.750,0,1},
			{213.500,0,1},
			{214.250,0,1},
			{215.000,3,1},
			{215.500,3,1},
			{216.000,0,1},
			{216.750,0,1},
			{217.500,0,1},
			{218.250,0,1},
			{219.000,3,1},
			{219.500,3,1},
		}
		
		for pn=1,2 do
			local f1 = 1
			local f2 = 1
			for _,v in pairs(kick) do
				if v[2] == 0 then
					me{v[1],.75,outCubic,startVal=2*f1*(pn*2-3),0,'adrunk',plr=pn}
					f1 = f1*-1
				else
					me{v[1],.75,outCubic,startVal=2*f2*(pn*2-3),0,'atipsy',plr=pn}
					f2 = f2*-1
				end
			end
		end
		
		local ddo = 0
		
		me{160-ddo-3,3,inOutCubic,1.4,'xmod'}
		set{160-ddo-.3,2,'beat'}
		set{208-ddo-.7,0,'beat'}

		me{160-ddo-2,4,inOutCubic,.2,'alpha',pn=1}
		me{208-ddo-2,4,inOutCubic,.5,'alpha',pn=1}
		me{219-ddo,2,inOutCubic,1,'alpha',pn=1}
		me{223-ddo,2,inOutCubic,.5,'alpha',pn=1}
		
		for pn=1,2 do
			me{160-ddo-2,4,inOutCubic,-310*(pn*2-3),'amovex',plr=pn}
			me{208-ddo-2,2,inCubic,0,'amovex',plr=pn}
			me{208-ddo,12,linear,-310*(pn*2-3),'amovex',plr=pn}
		end
		
		for i=160-ddo,208-ddo-1 do
			me{i,1,linear,1.7*112,startVal=0,'amovey'}
		end
		set{208-ddo,0,'amovey'}
		
		me{174-ddo,4,inOutCubic,180,'rotationz'}
		me{174-ddo,4,inOutCubic,-560,'camy'}
		me{190-ddo,4,inOutCubic,90,'rotationz'}
		me{190-ddo,4,inOutCubic,-280,'camy'}
		me{206-ddo,4,inOutCubic,0,'rotationz'}
		me{206-ddo,4,inOutCubic,0,'camy'}
		
		me{174-ddo,4,inOutCubic,-180,'confusion'}
		me{190-ddo,4,inOutCubic,-90,'confusion'}
		me{206-ddo,4,inOutCubic,0,'confusion'}
		
		me{174-ddo,4,inOutCubic,1,'flip'}
		me{190-ddo,4,inOutCubic,0,'flip'}
		
		me{216,4,linear,0,'alpha',pn=2}
		
		wig{220,16,4,outCubic,1,'drunk',pn=1}
		
		
		set{224,1,'alpha',pn=2}
		set{224,0,'stealth',pn=2}
		set{224,1,'dark',pn=2}
		me{234,3,linear,0,'dark',pn=2}
		
		
		
		local pieo = 0
		me{224-pieo,28,linear,3,'wave'}
		me{224-pieo,8,linear,1,'hidden'}
		set{224-pieo,50,'hiddenoffset'}
		me{240-pieo,12,linear,400,'hiddenoffset'}
		me{252-pieo,8,inOutCubic,2,'wave'}
		me{252-pieo,8,inOutCubic,.5,'drunk'}
		me{252-pieo,8,inOutCubic,200,'hiddenoffset'}
		
		for pn=1,2 do
			set{224-pieo,-310*(pn*2-3),'amovex',pn=pn}
		end
		set{224-pieo,.3,'alpha',pn=1}
		
		--kade "sway" again
		mpf{256-pieo,284-pieo,function(beat)
			
			for pn=1,2 do
				for col = 0,3 do
					activeMods[pn]['movex'..col] = 32 * math.sin((beat + col*0.25) * math.pi)
					local mu = col%2 == 0 and 1 or -1
					activeMods[pn]['movey'..col] = 32 * mu * math.tan((beat) * 0.25 * math.pi)
				end
			end
			
		end}
		
		me{280-pieo,8,inOutCubic,0,'wave'}
		me{280-pieo,8,inOutCubic,0,'drunk'}
		me{276-pieo,8,inOutCubic,0,'hidden'}
		set{288,0,'hiddenoffset'}
		
		local f = 1
		for i=256-pieo,283-pieo,4 do
			for pn = 1,2 do
				me{i,2,outCubic,startVal=360*f,0,'confusion',pn=spn}
			end
			f = f*-1
		end
		
		for pn=1,2 do
			me{284-pieo,4,inCubic,0,'amovex',plr=pn}
		end
		
		--NEW wrapping?
		--WITH NO FLICKER
		
		local fwo = 0
		
		me{288-fwo,8,inOutCubic,1280/8,'colspacing'}
		me{288-fwo,8,inOutCubic,(1280/8)*4,'spacing'}
		me{288-fwo,8,inOutCubic,-176,'addx'}
		me{288-fwo,8,inOutCubic,1,'waveamp'}
		
		me{316-fwo,4,inCubic,112,'colspacing'}
		me{316-fwo,4,inCubic,620,'spacing'}
		me{316-fwo,4,inCubic,0,'addx'}
		me{316-fwo,4,inCubic,0,'waveamp'}
		
		mpf{288-fwo,320-fwo,function(beat)
			
			local scrollerpos = beat-(288-fwo)
			
			for pn=1,2 do
			
				for col = 0,3 do
					local cpos = col*-112
					if pn == 2 then cpos = cpos-620 end
					local c = (pn-1)*4+col
					
					local colspacing = activeMods[pn].colspacing
					local spacing = activeMods[pn].spacing
					local addx = activeMods[pn].addx
					local waveamp = activeMods[pn].waveamp
					
					local newpos = (col*colspacing + (pn-1)*spacing + scrollerpos*1280/8)%1280 + addx
					
					activeMods[pn]['movex'..col] = cpos + newpos
				
					local ang = 2*math.pi*(c/8)
					activeMods[pn]['reverse'..col] = waveamp * .1 * math.sin(beat*math.pi + ang)
					
				end
				
			end
			
		end}
		
		
		--[[
		me{288-fwo,8,outExpo,startVal=620,448,'spacing'}
		me{288-fwo,8,outExpo,startVal=0,224,'addx'}
		me{288-fwo,8,outExpo,startVal=0,1,'waveamp'}
		
		me{316-fwo,4,inCubic,0,'addx'}
		me{316-fwo,4,inCubic,0,'waveamp'}
		me{316-fwo,4,inCubic,620,'spacing'}
		
		mpf{288-fwo,320-fwo,function(beat)
			
			local f = flicker()
			if f == -1 then f = 0 end
			
			local scrollerpos = beat-288 
			
			for pn=1,2 do
				local cpos = 0
				if pn == 2 then cpos = -620 end
				
				local spacing = activeMods[pn].spacing
				local addx = activeMods[pn].addx
				local waveamp = activeMods[pn].waveamp
				
				local newpos = (spacing * (pn - 1 + f*2) + scrollerpos * spacing * 0.25)%(spacing*4) - spacing*2 + addx
				
				activeMods[pn].movex = cpos + newpos
				for col = 0,3 do
					local ang = 2*math.pi*((pn-1)*4 + col)/8
					activeMods[pn]['reverse'..col] = waveamp * .1 * math.sin(beat*math.pi + ang)
				end
				
			end
			
		end}
		]]
		
		set{315.5,3,'beat'}
		set{319.5,0,'beat'}
		
		
		
		
		local wso2 = 0
		
		hide{(320-wso2)-1,pn=2}
		unhide{(320-wso2)+16,pn=2}
		hide{(320-wso2)+16,pn=1}
		unhide{(320-wso2)+32,pn=1}
		
		me{(320-wso2)-3,3,linear,1.7,'xmod'}
		--me{(320-wso2)-3,3,linear,1,'drunk'}
		--set{(320-wso2)-3,1,'drunk'}
		
		me{(320-wso2)-2,2,linear,1,'alpha',pn=1}
		
		me{(320-wso2)-4,3,linear,50,'movey'}
		me{(320-wso2),3,linear,-.1,'flip'}
		
		me{(320-wso2)+28,4,inCubic,0,'flip'}
		
		for i=0,1 do
			me{(320-wso2)+8*i,4,linear,448*1.2,'movex',pn=1}
			me{(320-wso2)+4+8*i,4,linear,0,'movex',pn=1}
			me{(320-wso2)+16+8*i,4,linear,-448*1.2,'movex',pn=2}
			me{(320-wso2)+16+4+8*i,4,linear,0,'movex',pn=2}
		end
		
		mpf{(320-wso2),(320-wso2)+64,function(beat)
			
			local pn=1
			if (beat-(320-wso2)) % 32 > 16 then pn = 2 end
			
			local pingpong = beat - math.floor(beat)
			if beat%2 > 1 then
				pingpong = 1-pingpong
			end
			
			local wdir = -1 * (pn*2-3)
			if (beat-(320-wso2)) % 8 > 4 then wdir = wdir*-1 end
			
			if beat < (320-wso2)+32 then
				for col=0,3 do
					activeMods[pn]['amovey'..col] = math.min(-112 * wdir * math.sin(beat*math.pi + col*math.pi),0)
				end
				activeMods[pn].invert = pingpong*1.2
			else
				for col=0,3 do
					activeMods[pn]['amovey'..col] = -56*1.2 * wdir * math.sin(beat*math.pi + col*math.pi)
				end
				activeMods[pn].invert = 0.6 - 0.6*math.cos(beat*math.pi)
			end
			
		end}
		
		for i=320-wso2,351-wso2,8 do
			local plr = 2
			if (i-(320-wso2)) % 32 < 16 then plr = 1 end
			if ALLOW_REVERSE then
				me{i+3,2,inOutExpo,1,'reverse',pn=plr}
				me{i+7,2,inOutExpo,0,'reverse',pn=plr}
			end
		end
		
		local bso = 0
		
		me{352-bso-1,2,linear,.5,'alpha',pn=1}
		
		--sm2 = b,len,eas,amt,mods,intime
		for i=352-bso,383-bso,4 do
		
			sm2{i,2,outCubic,2,'drunk'}
			sm2{i,2,outCubic,1,'brake'}
			
			for pn=1,2 do
			
				me{i+1,1,inOutExpo,-310*(pn*2-3),'amovex',plr=pn}
				if i ~= 380-bso then
					me{i+1,1,inOutCubic,1,'invert',plr=pn}
					me{i+2.5,1,inOutExpo,0,'amovex',plr=pn}
					me{i+2.5,1,inOutCubic,0,'invert',plr=pn}
					
					for c=0,3 do
						local mu = c%2 == 0 and 1 or -1
						me{i+1,.5,outCubic,.1*mu,'reverse'..c,plr=pn}
						me{i+1.5,.5,inCubic,0,'reverse'..c,plr=pn}
						
						me{i+2.5,.5,outCubic,-.1*mu,'reverse'..c,plr=pn}
						me{i+3,.5,inCubic,0,'reverse'..c,plr=pn}
					end
				end
				
			end
		end
		
		
		local ddo2 = 0
			
		me{384-ddo2-3,3,inOutCubic,1.4,'xmod'}
		set{384-ddo2-.3,2,'beat'}
		set{416-ddo2-.7,0,'beat'}

		me{384-ddo2-2,4,inOutCubic,.25,'alpha',pn=1}
		me{416-ddo2-2,4,inOutCubic,.45,'alpha',pn=1}
		
		for pn=1,2 do
			set{384-ddo2,-310*(pn*2-3),'amovex',pn=pn}
			me{416-ddo2-2,2,inCubic,0,'amovex',pn=pn}
			me{428-ddo2,12,linear,-310*(pn*2-3),'amovex',pn=pn}
		end
		
		for i=384-ddo2,416-ddo2-1 do
			me{i,1,linear,1.7*112,startVal=0,'amovey'}
		end
		set{416-ddo2,0,'amovey'}
		
		set{384-ddo2,8,'camwag'}
		me{414-ddo2,2,linear,0,'camwag'}
		
		me{391-ddo2,2,inOutCubic,180,'rotationz'}
		me{391-ddo2,2,inOutCubic,-560,'camy'}
		
		me{399-ddo2,2,inOutCubic,0,'rotationz'}
		me{399-ddo2,2,inOutCubic,0,'camy'}
		
		me{407-ddo2,2,inOutCubic,-90,'rotationz'}
		me{407-ddo2,2,inOutCubic,-280,'camy'}
		
		me{415-ddo2,2,inOutCubic,0,'rotationz'}
		me{415-ddo2,2,inOutCubic,0,'camy'}
		
		me{391-ddo2,2,inOutCubic,-180,'confusion'}
		me{399-ddo2,2,inOutCubic,0,'confusion'}
		me{407-ddo2,2,inOutCubic,90,'confusion'}
		me{415-ddo2,2,inOutCubic,0,'confusion'}
		
		me{391-ddo2,2,inOutCubic,1,'flip'}
		me{399-ddo2,2,inOutCubic,0,'flip'}
		
		me{415-ddo2,2,inOutCubic,1.7,'xmod'}
		me{436-ddo2,8,linear,0,'alpha',plr=1}
		me{448-ddo2,2,linear,.3,'alpha',plr=1}
		
		
		wig{440-ddo2,16,4,outCubic,1,'drunk'}
		wig{444-ddo2,16,4,outCubic,2,'drunk'}
		
		set{416-ddo2,32,'drawsize'}
		
		if ALLOW_REVERSE then
			for i=0,1 do
				me{416+8*i-ddo2,4,linear,1,'reverse'}
				me{420+8*i-ddo2,4,linear,0,'reverse'}
				me{432+4*i-ddo2,2,linear,1,'reverse'}
				me{434+4*i-ddo2,2,linear,0,'reverse'}
			end
		end
		
		set{448-ddo2,32,'drawsize'}
		
		kick2 = {
			{416.000,0,1},
			{416.750,0,1},
			{417.500,0,1},
			{418.250,0,1},
			{419.000,3,1},
			{419.500,3,1},
			{420.000,0,1},
			{420.750,0,1},
			{421.500,0,1},
			{422.250,0,1},
			{423.000,3,1},
			{423.500,3,1},
			{424.000,0,1},
			{424.750,0,1},
			{425.500,0,1},
			{426.250,0,1},
			{427.000,3,1},
			{427.500,3,1},
			{428.000,0,1},
			{428.750,0,1},
			{429.500,0,1},
			{430.250,0,1},
			{431.000,3,1},
			{431.500,3,1},
			{432.000,0,1},
			{432.500,0,1},
			{433.000,0,1},
			{433.500,0,1},
			{434.000,0,1},
			{434.500,0,1},
			{435.000,0,1},
			{435.500,0,1},
			{436.000,0,1},
			{436.500,0,1},
			{437.000,0,1},
			{437.500,0,1},
			{438.000,0,1},
			{438.500,0,1},
			{439.000,0,1},
			{439.500,0,1},
		}
		
		for pn=1,2 do
			local f1 = 1
			local f2 = 1
			for _,v in pairs(kick2) do
				if v[2] == 0 then
					me{v[1]-ddo2,.75,outCubic,startVal=2*f1*(pn*2-3),0,'adrunk',plr=pn}
					f1 = f1*-1
				else
					me{v[1]-ddo2,.75,outCubic,startVal=2*f2*(pn*2-3),0,'atipsy',plr=pn}
					f2 = f2*-1
				end
			end
		end
		
		
		--NEW wrapping?
		--NO FLICKER
		
		local fwo2 = 0
		
		set{448-fwo2,.4,'alpha',pn=1}
		set{448-fwo2,32,'drawsize'}
		me{472-fwo2,4,inExpo,.5,'flip'}
		
		for pn=1,2 do
			set{448-fwo2,0,'amovex',pn=pn}
			set{476-fwo2,-310*(pn*2-3),'amovex',pn=pn}
		end
		
		me{448-fwo2,4,outCubic,startVal=112,1280/8,'colspacing'}
		--me{448-fwo2,8,inOutCubic,startVal=0,0,'colspacing'}
		me{448-fwo2,4,outCubic,startVal=0,(1280/8)*4,'spacing'}
		me{448-fwo2,4,outCubic,startVal=310,-176,'addx'}
		me{448-fwo2,8,inOutCubic,1,'waveamp'}
		
		me{472-fwo2,4,inCubic,112,'colspacing'}
		me{472-fwo2,4,inCubic,0,'spacing'}
		me{472-fwo2,4,inCubic,-310,'addx'}
		me{472-fwo2,4,inCubic,0,'waveamp'}
		
		mpf{448-fwo2,476-fwo2,function(beat)
			
			local scrollerpos = beat-(448-fwo2)
			
			for pn=1,2 do
			
				for col = 0,3 do
					local cpos = col*-112 + 0
					if pn == 2 then cpos = cpos-620 end
					local c = (pn-1)*4+col
					
					local colspacing = activeMods[pn].colspacing
					local spacing = activeMods[pn].spacing
					local addx = activeMods[pn].addx
					local waveamp = activeMods[pn].waveamp
					
					local newpos = (col*colspacing + (pn-1)*spacing + scrollerpos*1280/8)%1280 + addx
					
					activeMods[pn]['movex'..col] = cpos + newpos
				
					local ang = 2*math.pi*(c/8)
					activeMods[pn]['reverse'..col] = waveamp * .1 * math.sin(beat*math.pi + ang)
					
				end
				
			end
			
		end}
		
		if ALLOW_REVERSE then
			for i=448-fwo2,463-fwo2,8 do
				me{i+3,2,inOutCubic,1,'reverse'}
				me{i+7,2,inOutCubic,0,'reverse'}
			end
			for i=464-fwo2,471-fwo2,4 do
				me{i+1.5,1,inOutCubic,1,'reverse'}
				me{i+3.5,1,inOutCubic,0,'reverse'}
			end
		end
		
		me{475-fwo2,2,linear,1,'alpha',plr=1}
		me{476-fwo2,8,outExpo,1,'halo'}
		
		mpf{476-fwo2,488-fwo2,function(beat)
			
			for pn=1,2 do
			
				local halo = activeMods[pn].halo
				
				for col=0,3 do
					
					local ang = 2*math.pi*((pn-1)*4 + col)/8 + ((beat-476)*0.5*math.pi)
					
					activeMods[pn]['movex'..col] = 280 * halo * math.sin(ang)
					activeMods[pn]['movey'..col] = 70 * (1-(activeMods[pn].amovey/300)) * halo * math.cos(ang)
					
				end
				
			end
			
		end}
		
			me{484-fwo2,4,inCubic,750,'amovey'}

		end
	end
	-------END-------
end
function onSongStart()
    
    TEMPLATE.songStart()
	
	--for i=0,7,1 do
	--	print('default position '..i..' = '..defaultPositions[i+1].x)
	--end
	
end

function onUpdatePost(elapsed)
	TEMPLATE.update(elapsed)
end

return 0