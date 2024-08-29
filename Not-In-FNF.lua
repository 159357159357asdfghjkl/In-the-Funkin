print('LOADED TEMPLATE')

TEMPLATE = {}

--[[
			READ MEEEEEEEE!
	This version of modchart template is modified from TaroNuke's TEMPLATE 1, the file has many differences.
	By comparing this version with my original modified version, i removed stuff from other modchart tools,
	and add a ton of new stuff. e.g. it can use many mods from one ease
			now i will keep updating this shit instead of updating the old one
]]



-- Taro's janky "TEMPLATE 1" implementation for FNF
-- shoutouts to Kade for letting me do this

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

function instant()
	return 1
end

function scale(x, l1, h1, l2, h2)
	return (((x) - (l1)) * ((h2) - (l2)) / ((h1) - (l1)) + (l2))
end

function math.clamp(val,min,max)
	if val < min then return min end
	if val > max then return max end
	return val
end

function lerp(a, b, t)
    return a + t * (b-a)
end
function square(angle)
	local fAngle = angle % (math.pi * 2)
		--Hack: This ensures the hold notes don't flicker right before they're hit.
		if fAngle < 0.01 then
		    fAngle = fAngle + math.pi * 2
		end
	return fAngle >= math.pi and -1.0 or 1.0;
end

function triangle( angle )
	local fAngle= angle % math.pi * 2.0
	if fAngle < 0.0 then
		fAngle= fAngle+math.pi * 2.0
	end
	local result= fAngle * (1 / math.pi)
	if result < .5 then
		return result * 2.0
	elseif result < 1.5 then
		return 1.0 - ((result - .5) * 2.0)
	else
		return -4.0 + (result * 2.0)
	end
	
end

function quantize(f,interval)
return int((f+interval/2)/interval)*interval
end

function rotationXYZ( rX, rY, rZ )
    local PI=math.pi
	rX = rX*(PI/180)
	rY = rY*(PI/180)
	rZ = rZ*(PI/180)

	local cX = math.cos(rX)
	local sX = math.sin(rX)
	local cY = math.cos(rY)
	local sY = math.sin(rY)
	local cZ = math.cos(rZ)
	local sZ = math.sin(rZ)

	 return {
	 	m00=cZ*cY, m01=cZ*sY*sX+sZ*cX, m02=cZ*sY*cX+sZ*(-sX), m03=0,
	 	m10=(-sZ)*cY, m11=(-sZ)*sY*sX+cZ*cX, m12=(-sZ)*sY*cX+cZ*(-sX), m13=0,
	 	m20=-sY, m21=cY*sX, m22=cY*cX, m23=0,
	 	m30=0, m31=0, m32=0, m33=1
	  }
end

function selectTanType(angle,is_csc)
if is_csc ~= 0 then
return fastCsc(angle)
else
return fastTan(angle)
end
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
	drunkspeed = 0,
	drunkoffset = 0,
	drunkperiod = 0,
	tipsy = 0,
	tipsyspeed = 0,
	tipsyoffset = 0,
	tandrunk = 0,
	tandrunkspeed = 0,
	tandrunkoffset = 0,
	tandrunkperiod = 0,
	drunky = 0,
	drunkyspeed = 0,
	drunkyoffset = 0,
	drunkyperiod = 0,
	tandrunky = 0,
	tandrunkyspeed = 0,
	tandrunkyoffset = 0,
	tandrunkyperiod = 0,
	tipsyx = 0,
	tipsyxspeed = 0,
	tipsyxoffset = 0,
	tantipsy = 0,
	tantipsyspeed = 0,
	tantipsyoffset = 0,
	tantipsyx = 0,
	tantipsyxspeed = 0,
	tantipsyxoffset = 0,
	adrunk = 0, --non conflict accent mod
	atipsy = 0, --non conflict accent mod
	bumpyx = 0,
	bumpyxoffset = 0,
	bumpyxperiod = 0,
	bumpyy = 0,
	bumpyyoffset = 0,
	bumpyyperiod = 0,
	beaty = 0,
	sawtooth = 0,
	sawtoothperiod = 0,
	digital = 0,
	digitalsteps = 0,
	digitaloffset = 0,
	digitalperiod = 0,
	tandigital = 0,
	tandigitalsteps = 0,
	tandigitaloffset = 0,
	tandigitalperiod = 0,
	square = 0,
	squareoffset = 0,
	squareperiod = 0,
	bounce = 0,
	bounceoffset = 0,
	bounceperiod = 0,
	bouncey = 0,
	bounceyoffset = 0,
	bounceyperiod = 0,
	xmode = 0,
	tiny = 0,
	zigzag = 0,
	zigzagoffset = 0,
	zigzagperiod = 0,
	attenuatex = 0,
	attenuatey = 0,
	swap = 0,
	parabolax = 0,
	parabolay = 0,
	tornado = 0,
	scale = 0,
	zoom = 0,
	big = 0,
	scalex = 0,
	scaley = 0,
	squish = 0,
	stretch = 0,
	pulseinner = 0,
	pulseouter = 0,
	pulseoffset = 0,
	pulseperiod = 0,
	shrinkmult = 0,
	shrinklinear = 0,
	movex = 0,
	movey = 0,
	movez = 0,
	amovex = 0,
	amovez = 0,
	amovey = 0,
	reverse = 0,
	split = 0,
	cross = 0,
	centered = 0,
	dark = 0,
	stealth = 0,
	alpha = 1,
	randomvanish = 0,
	blink = 0,
	confusion = 0,
	dizzy = 0,
	wave = 0,
	waveperiod = 0,
	brake = 0,
	boost = 0,
	boomerang = 0,
	expand = 0,
	expandperiod = 0,
	tanexpand = 0,
	tanexpandperiod = 0,
	hidden = 0,
	hiddenoffset = 0,
	sudden = 0,
	suddenoffset = 0,
	alternate = 0,
	inside = 0,
	outside = 0,
	camx = 0,
	camy = 0,
	camalpha = 1,
	camzoom = 0,
	rotationz = 0,
	camwag = 0,
	xmod = 1, --scrollSpeed
	cosecant = 0,
	vibrate = 0,
	bumpy = 0,
	bumpyoffset = 0,
	bumpyperiod = 0,
	tanbumpy = 0,
	tanbumpyoffset = 0,
	tanbumpyperiod = 0,
	attenuatez = 0,
	tornadoz = 0,
	parabolaz = 0,
	sawtoothz = 0,
	sawtoothzperiod = 0,
	squarez = 0,
	squarezoffset = 0,
	squarezperiod = 0,
	bouncez = 0,
	bouncezoffset = 0,
	bouncezperiod = 0,
	digitalz = 0,
	digitalzsteps = 0,
	digitalzoffset = 0,
	digitalzperiod = 0,
	tandigitalz = 0,
	tandigitalzsteps = 0,
	tandigitalzoffset = 0,
	tandigitalzperiod = 0,
	zigzagz = 0,
	zigzagzoffset = 0,
	zigzagzperiod = 0,
	beatz = 0,
	drunkz = 0,
	drunkzspeed = 0,
	drunkzoffset = 0,
	drunkzperiod = 0,
	tandrunkz = 0,
	tandrunkzspeed = 0,
	tandrunkzoffset = 0,
	tandrunkzperiod = 0,
	tipsyz = 0,
	tipsyzspeed = 0,
	tipsyzoffset = 0,
	tantipsyz = 0,
	tantipsyzspeed = 0,
	tantipsyzoffset = 0,
	confusionx = 0,
	confusionxoffset = 0,
	confusiony = 0,
	confusionyoffset = 0,
	confusionoffset = 0,
	twirl = 0,
	roll = 0,
	receptorscroll = 0,
}

--column specific mods
for i=0,3 do
	modList['movex'..i] = 0
	modList['movey'..i] = 0
	modList['movez'..i] = 0
	modList['amovex'..i] = 0
	modList['amovey'..i] = 0
	modList['amovez'..i] = 0
	modList['dark'..i] = 0
	modList['stealth'..i] = 0
	modList['confusion'..i] = 0
	modList['confusionx'..i] = 0
	modList['confusionxoffset'..i] = 0
	modList['confusiony'..i] = 0
	modList['confusionyoffset'..i] = 0
	modList['confusionoffset'..i] = 0
	modList['reverse'..i] = 0
	modList['tiny'..i] = 0
	modList['scale'..i] = 0
	modList['zoom'..i] = 0
	modList['scalex'..i] = 0
	modList['scaley'..i] = 0
	modList['squish'..i] = 0
	modList['stretch'..i] = 0
	modList['bumpy'..i] = 0
	modList['xmod'..i] = 1 --column specific scrollSpeed multiplier
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
	if m.sudden ~= 0 then
		if fYOffset > m.suddenoffset and fYOffset <= m.suddenoffset+200 then
			local hmult = -(fYOffset-m.suddenoffset)/200
			alp = alp*(1-hmult)*m.sudden
		elseif fYOffset > m.suddenoffset+100 then
			alp = alp*(1-m.sudden)
		end
	end
	if m.blink ~= 0 then
		local time = getSongPosition()/1000
		local f = math.sin(time*10)
		f=quantize(f,0.3333)
		alp = alp + scale( f, 0, 1, -1, 0 );
	end
	if m.randomvanish ~= 0 then
		local fRealFadeDist = 80;
		alp = alp + scale( math.abs(fYOffset-360), fRealFadeDist, 2*fRealFadeDist, -1, 0 ) * m.randomvanish;
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
	if m.centered ~= 0 then val = scale( m.centered, 0, 1, val, 0 ) end

	return val
end

function getYAdjust(fYOffset, iCol, pn)
	
	local m = activeMods[pn]
	
	local yadj = 0
	local fScrollSpeed = 1
	if m.wave ~= 0 then
		yadj =yadj + m.wave * 20 *math.sin( fYOffset/((m.waveperiod*38)+38) );
	end
	
	if m.brake ~= 0 then

		local fEffectHeight = 500;
		local fScale = scale( fYOffset, 0, fEffectHeight, 0, 1 )
		local fNewYOffset = fYOffset * fScale; 
		local fBrakeYAdjust = m.brake * (fNewYOffset - fYOffset)
		
		fBrakeYAdjust = math.clamp( fBrakeYAdjust, -400, 400 )
		yadj = yadj+fBrakeYAdjust;
	
	end

	if m.boost ~= 0 then

		local fEffectHeight = 500;
		local fNewYOffset = fYOffset * 1.5 / ((fYOffset+fEffectHeight/1.2)/fEffectHeight); 
		local fAccelYAdjust = m.boost * (fNewYOffset - fYOffset)
		
		fAccelYAdjust = math.clamp( fAccelYAdjust, -400, 400 )
		yadj = yadj+fAccelYAdjust;
	
	end
	
    fYOffset = fYOffset+yadj
    
	if m.boomerang ~= 0 then
		fYOffset = ((-1*fYOffset*fYOffset/500) + 1.5*fYOffset)*m.boomerang
	end
	if m.expand ~= 0 then
	local last = 0
	local time = getSongPosition() / 1000
	local expandSeconds = 0
    expandSeconds = expandSeconds + (time - last);
    expandSeconds = expandSeconds % ((math.pi * 2) / (m.expandperiod + 1));
    last = time
	    local fExpandMultiplier = scale(math.cos(expandSeconds * 3 * (m.expandperiod + 1)), -1, 1, 0.75, 1.75);
      fScrollSpeed = fScrollSpeed * scale(m.expand, 0, 1, 1, fExpandMultiplier);
	end
	if m.tanexpand ~= 0 then
	local last = 0
	local time = getSongPosition() / 1000
	local tanExpandSeconds = 0
    tanExpandSeconds = tanExpandSeconds + (time - last);
    tanExpandSeconds = tanExpandSeconds % ((math.pi * 2) / (m.tanexpandperiod + 1));
    last = time
	    local fTanExpandMultiplier = scale(selectTanType(tanExpandSeconds * 3 * (m.tanexpandperiod + 1),m.cosecant), -1, 1, 0.75, 1.75);
      fScrollSpeed = fScrollSpeed * scale(m.tanexpand, 0, 1, 1, fTanExpandMultiplier);
	end
	fYOffset = fYOffset * fScrollSpeed
	return fYOffset
end

function getZoom(fYOffset, iCol, pn)
	local m = activeMods[pn]
    local fZoom = 1
	local fPulseInner = 1.0
	if m.pulseinner ~= 0 or m.pulseouter ~= 0 then
		fPulseInner = ((m.pulseinner*0.5)+1)
		if fPulseInner == 0 then
			fPulseInner = 0.01
		end
	end
    if m.pulseinner ~= 0 or m.pulseouter ~= 0 then
    local sine = math.sin(((fYOffset+(100.0*m.pulseoffset))/(0.4*(ARROW_SIZE+(m.pulseperiod*ARROW_SIZE)))))

		fZoom = fZoom*((sine*(m.pulseouter*0.5))+fPulseInner)
	end
	if m.shrinkmult ~= 0 and fYOffset >= 0 then
	fZoom = fZoom * (1/(1+(fYOffset*(m.shrinkmult/100.0))))
    end
    if m.shrinklinear ~= 0 and fYOffset >= 0 then
    fZoom = fZoom + (fYOffset*(0.5*m.shrinklinear/ARROW_SIZE))
    end
	local fTinyPercent = m.tiny
	if fTinyPercent ~= 0 then
		fTinyPercent = math.pow( 0.5, fTinyPercent )
		fZoom = fZoom * fTinyPercent
	end
	if m['tiny'..iCol] ~= 0 then
		fTinyPercent = math.pow( 0.5,  m['tiny'..iCol] )
		fZoom = fZoom * fTinyPercent
	end
	    fZoom = fZoom + m.zoom + m['zoom'..iCol]
	return fZoom
end

--im fucked
function receptorRotation(fYOffset,iCol,pn)
    local fRotationX, fRotationY, fRotationZ = 0, 0, 0
    local m = activeMods[pn]
    if m['confusionx'..iCol]~= 0 then
		fRotationX = fRotationX +m['confusionx'..iCol] * 180.0/math.pi;
	end
	if m.confusionxoffset ~= 0 then
		fRotationX = fRotationX +m.confusionxoffset * 180.0/math.pi;
	end
	if m['confusionxoffset'..iCol] ~= 0 then
		fRotationX = fRotationX +m['confusionxoffset'..iCol] * 180.0/math.pi;
	end
	if m.confusionx ~= 0 then
		local fConfRotation = beat
		local PI = math.pi
		fConfRotation = fConfRotation * m.confusionx
		fConfRotation = fConfRotation % ( 2*PI )
		fConfRotation = fConfRotation*(-180/PI)
		fRotationX = fRotationX + fConfRotation;
	end
	if m['confusiony'..iCol] ~= 0 then
		fRotationY = fRotationY + (m['confusiony'..iCol] * 180.0/math.pi)
	end
	if m['confusionyoffset'..iCol] ~= 0 then
		fRotationY = fRotationY + (m['confusionyoffset'..iCol] * 180.0/math.pi)
	end
	if m.confusionyoffset ~= 0 then
		fRotationY = fRotationY + m.confusionyoffset * 180.0/math.pi
	end
	if m.confusiony ~= 0 then
		local fConfRotation = beat
		local PI=math.pi
		fConfRotation = fConfRotation * m.confusiony
		fConfRotation = fConfRotation%(2*PI )
		fConfRotation = fConfRotation *(-180/PI)
		fRotationY = fRotationY + fConfRotation
	end
	if m['confusion'..iCol] ~= 0 then
		fRotationZ = fRotationZ+m['confusion'..iCol] * 180.0/math.pi
	end
	if m.confusionoffset ~= 0 then
		fRotationZ = fRotationZ+m.confusionoffset * 180.0/math.pi
	end
	if m['confusionoffset'..iCol] ~= 0 then
		fRotationZ = fRotationZ+m['confusionoffset'..iCol] * 180.0/math.pi
	end
	if m.confusion ~= 0 then
		local fConfRotation = beat
		local PI=math.pi
		fConfRotation = fConfRotation * m.confusion
		fConfRotation = fConfRotation%( 2*PI );
		fConfRotation = fConfRotation*(-180/PI)
		fRotationZ = fRotationZ*fConfRotation;
	end
    return fRotationX, fRotationY, fRotationZ
end

function arrowRotation(fYOffset,iCol,pn,noteBeat)
    local fRotationX, fRotationY, fRotationZ = 0, 0, 0
    local m = activeMods[pn]
    if m['confusionx'..iCol]~= 0 then
		fRotationX = fRotationX +m['confusionx'..iCol] * 180.0/math.pi;
	end
	if m.confusionxoffset ~= 0 then
		fRotationX = fRotationX +confusionxoffset * 180.0/math.pi;
	end
	if m['confusionxoffset'..iCol] ~= 0 then
		fRotationX = fRotationX +m['confusionxoffset'..iCol] * 180.0/math.pi;
	end
	if m.confusionx ~= 0 then
		local fConfRotation = beat
		fConfRotation = fConfRotation * m.confusionx
		fConfRotation = fConfRotation % ( 2*PI )
		fConfRotation = fConfRotation*(-180/PI)
		fRotationX = fRotationX + fConfRotation;
	end
	if m.roll ~= 0 then
		fRotationX = fRotationX + (m.roll * fYOffset/2)
	end
	if m['confusiony'..iCol] ~= 0 then
		fRotationY = fRotationY + (m['confusiony'..iCol] * 180.0/math.pi)
	end
	if m['confusionyoffset'..iCol] ~= 0 then
		fRotationY = fRotationY + (m['confusionyoffset'..iCol] * 180.0/math.pi)
	end
	if m.confusionyoffset ~= 0 then
		fRotationY = fRotationY + m.confusionyoffset * 180.0/math.pi
	end
	if m.confusiony ~= 0 then
		local fConfRotation = beat
		local PI=math.pi
		fConfRotation = fConfRotation * m.confusiony
		fConfRotation = fConfRotation%(2*PI )
		fConfRotation = fConfRotation *(-180/PI)
		fRotationY = fRotationY + fConfRotation
	end
	if m.twirl ~= 0 then
		fRotationY = fRotationY + (m.twirl * fYOffset/2)
	end

	if m.dizzy ~= 0 then
	local fSongBeat = beat
	local PI = math.pi
	local fDizzyRotation = noteBeat - fSongBeat
		fDizzyRotation = fDizzyRotation * m.dizzy
		fDizzyRotation = fDizzyRotation % (2*PI)
		fDizzyRotation = fDizzyRotation*(180/PI)
		fRotationZ = fRotationZ+fDizzyRotation
	end
		if m['confusion'..iCol] ~= 0 then
		fRotationZ = fRotationZ+m['confusion'..iCol] * 180.0/math.pi
	end
	if m.confusionoffset ~= 0 then
		fRotationZ = fRotationZ+m.confusionoffset * 180.0/math.pi
	end
	if m['confusionoffset'..iCol] ~= 0 then
		fRotationZ = fRotationZ+m['confusionoffset'..iCol] * 180.0/math.pi
	end
	if m.confusion ~= 0 then
		local fConfRotation = beat
		local PI=math.pi
		fConfRotation = fConfRotation * m.confusion
		fConfRotation = fConfRotation%( 2*PI );
		fConfRotation = fConfRotation*(-180/PI)
		fRotationZ = fRotationZ*fConfRotation;
	end
    return fRotationX, fRotationY, fRotationZ
end

function getScale(fYOffset, iCol, pn, sx, sy)
    local x = sx
    local y = sy
    local m = activeMods[pn]
    x = x + m.scalex + m['scalex'..iCol] + m.scale + m['scale'..iCol] + m.big
    y = y + m.scaley + m['scaley'..iCol] + m.scale + m['scale'..iCol] + m.big
    local angle = 0;
    local stretch = m.stretch + m['stretch'..iCol]
    local squish = m.squish + m['squish'..iCol]
    local stretchX = lerp(1, 0.5, stretch);
    local stretchY = lerp(1, 2, stretch);

    local squishX = lerp(1, 2, squish);
    local squishY = lerp(1, 0.5, squish);
	x = x * ((math.sin(angle * math.pi / 180) * squishY) + (math.cos(angle * math.pi / 180) * squishX));
	x = x * ((math.sin(angle * math.pi / 180) * stretchY) + (math.cos(angle * math.pi / 180) * stretchX));

	y = y * ((math.cos(angle * math.pi / 180) * stretchY) + (math.sin(angle * math.pi / 180) * stretchX));
	y = y * ((math.cos(angle * math.pi / 180) * squishY) + (math.sin(angle * math.pi / 180) * squishX));
		return x,y
end
function arrowEffects(fYOffset, iCol, pn)
    local m = activeMods[pn]
	
    local xpos, ypos, rotz, zpos = 0, 0, 0, 0
	
	if m['confusion'..iCol] ~= 0 or m.confusion ~= 0 ~= 0 then
		rotz = rotz + m['confusion'..iCol] + m.confusion
	end
	if m.dizzy ~= 0 then
		rotz = rotz + m.dizzy*fYOffset
	end
    if m.drunk ~= 0 then
        xpos = xpos + m.drunk * math.cos(getSongPosition()*0.001 * (1 + m.drunkspeed) + iCol * ((m.drunkoffset * 0.2) + 0.2) + fYOffset * ((m.drunkperiod * 10) + 10) / screenHeight) * ARROW_SIZE * 0.5;
    end
    if m.drunkz ~= 0 then
        zpos = zpos + m.drunkz * math.cos(getSongPosition()*0.001 * (1 + m.drunkzspeed) + iCol * ((m.drunkzoffset * 0.2) + 0.2) + fYOffset * ((m.drunkzperiod * 10) + 10) / screenHeight) * ARROW_SIZE * 0.5;
    end
    if m.drunky ~= 0 then
        ypos = ypos + m.drunky * math.cos(getSongPosition()*0.001 * (1 + m.drunkyspeed) + iCol * ((m.drunkyoffset * 0.2) + 0.2) + fYOffset * ((m.drunkyperiod * 10) + 10) / screenHeight) * ARROW_SIZE * 0.5;
    end
    if m.tandrunk ~= 0 then
    xpos = xpos + m.tandrunk * selectTanType(getSongPosition()*0.001 * (1 + m.tandrunkspeed) + iCol * ((m.tandrunkoffset * 0.2) + 0.2) + fYOffset * ((m.tandrunkperiod * 10) + 10) / screenHeight,m.cosecant) * ARROW_SIZE * 0.5;
    end
    if m.tandrunkz ~= 0 then
    zpos = zpos + m.tandrunkz * selectTanType(getSongPosition()*0.001 * (1 + m.tandrunkzspeed) + iCol * ((m.tandrunkzoffset * 0.2) + 0.2) + fYOffset * ((m.tandrunkzperiod * 10) + 10) / screenHeight,m.cosecant) * ARROW_SIZE * 0.5;
    end
    if m.tandrunky ~= 0 then
    ypos = ypos + m.tandrunky * selectTanType(getSongPosition()*0.001 * (1 + m.tandrunkyspeed) + iCol * ((m.tandrunkyoffset * 0.2) + 0.2) + fYOffset * ((m.tandrunkyperiod * 10) + 10) / screenHeight,m.cosecant) * ARROW_SIZE * 0.5;
    end
    if m.tipsy ~= 0 then
        ypos = ypos + m.tipsy * math.cos( getSongPosition() * 0.001 * ((m.tipsyspeed * 1.2) + 1.2) + (iCol * ((m.tipsyoffset * 1.8) + 1.8)))* ARROW_SIZE * 0.4
    end
    if m.tantipsy ~= 0 then
        ypos = ypos + m.tantipsy * selectTanType( getSongPosition() * 0.001 * ((m.tantipsyspeed * 1.2) + 1.2) + (iCol * ((m.tantipsyoffset * 1.8) + 1.8)),m.cosecant)* ARROW_SIZE * 0.4
    end
    if m.tipsyz ~= 0 then
        zpos = zpos + m.tipsyz * math.cos( getSongPosition() * 0.001 * ((m.tipsyzspeed * 1.2) + 1.2) + (iCol * ((m.tipsyzoffset * 1.8) + 1.8)))* ARROW_SIZE * 0.4
    end
    if m.tantipsyz ~= 0 then
        zpos = zpos + m.tantipsyz * selectTanType( getSongPosition() * 0.001 * ((m.tantipsyzspeed * 1.2) + 1.2) + (iCol * ((m.tantipsyzoffset * 1.8) + 1.8)),m.cosecant)* ARROW_SIZE * 0.4
    end
    if m.tipsyx ~= 0 then
        xpos = xpos + m.tipsyx * math.cos( getSongPosition() * 0.001 * ((m.tipsyxspeed * 1.2) + 1.2) + (iCol * ((m.tipsyxoffset * 1.8) + 1.8))) * ARROW_SIZE * 0.4
    end
     if m.tantipsyx ~= 0 then
        xpos = xpos + m.tantipsyx * selectTanType( getSongPosition() * 0.001 * ((m.tantipsyxspeed * 1.2) + 1.2) + (iCol * ((m.tantipsyxoffset * 1.8) + 1.8)),m.cosecant)* ARROW_SIZE * 0.4
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
	if m['movez'..iCol] ~= 0 or m.movez ~= 0 then
	    zpos = zpos + m['movez'..iCol] + m.movez
	end
	if m['amovez'..iCol] ~= 0 or m.amovez ~= 0 then
		zpos = zpos + m['amovez'..iCol] + m.amovez
	end
	
	if m['reverse'..iCol] ~= 0 or m.reverse ~= 0 or m.split ~= 0 or m.cross ~= 0 or m.alternate ~= 0 or m.centered ~= 0 then
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

	if m.beaty ~= 0 then
			
		local fBeatStrength = m.beaty;
		
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
			
			ypos = ypos + fBeatStrength * fShift
			
		end
	
	end

	if m.sawtooth ~= 0 then
		xpos = xpos + (m.sawtooth*ARROW_SIZE) * ((0.5 / (m.sawtoothperiod+1) * fYOffset) / ARROW_SIZE - math.floor((0.5 / (m.sawtoothperiod+1) * fYOffset) / ARROW_SIZE) );
	end

	if m.digital ~= 0 then
		xpos = xpos + (m.digital * ARROW_SIZE * 0.5) * round((m.digitalsteps+1) * math.sin(math.pi * (fYOffset + (1.0 * m.digitaloffset ) ) / (ARROW_SIZE + (m.digitalperiod * ARROW_SIZE) )) )/(m.digitalsteps+1);
	end
	if m.tandigital ~= 0 then
		xpos = xpos + (m.tandigital * ARROW_SIZE * 0.5) * round((m.tandigitalsteps+1) * selectTanType(math.pi * (fYOffset + (1.0 * m.tandigitaloffset ) ) / (ARROW_SIZE + (m.tandigitalperiod * ARROW_SIZE) ),m.cosecant) )/(m.tandigitalsteps+1);
	end
	if m.bumpyx ~= 0 then
		xpos = xpos + m.bumpyx * 40*math.sin((fYOffset+(100.0*m.bumpyxoffset))/((m.bumpyxperiod*16.0)+16.0));
	end

	if m.square ~= 0 then
		local fResult = square( (math.pi * (fYOffset+(1.0*(m.squareoffset))) / (ARROW_SIZE+(m.squareperiod*ARROW_SIZE))) );
		xpos = xpos + (m.square * ARROW_SIZE * 0.5) * fResult;
	end
    if m.bumpyy ~= 0 then
		ypos = ypos + m.bumpyy * 40*math.sin((fYOffset+(100.0*m.bumpyyoffset))/((m.bumpyyperiod*16.0)+16.0));
	end
	if m.bounce ~= 0 then
		local fBounceAmt = math.abs( math.sin( ( (fYOffset + (1.0 * (m.bounceoffset) ) ) / ( 60 + m.bounceperiod*60) ) ) );
		xpos = xpos + m.bounce * ARROW_SIZE * 0.5 * fBounceAmt;
	end

	if m.bouncey ~= 0 then
		local fBounceAmt = math.abs( math.sin( ( (fYOffset + (1.0 * (m.bounceyoffset) ) ) / ( 60 + m.bounceyperiod*60) ) ) );
		ypos = ypos + m.bouncey * ARROW_SIZE * 0.5 * fBounceAmt;
	end
	
	if m.xmode ~= 0 then
		xpos = xpos + m.xmode * (pn == 2 and -fYOffset or fYOffset)
	end

	if m.tiny ~= 0 then
		local fTinyPercent = m.tiny
		fTinyPercent = math.min( math.pow(0.5, fTinyPercent), 1.0 );
		xpos = xpos * fTinyPercent
	end

	if m.zigzag ~= 0 then
		local fResult = triangle( (math.pi * (1/(m.zigzagperiod+1)) * 
		((fYOffset+(100.0*(m.zigzagoffset)))/ARROW_SIZE) ) );
	    
		xpos = xpos + (m.zigzag*ARROW_SIZE/2) * fResult;
	end

    if m.swap ~= 0 then
        xpos = xpos + screenWidth / 2 * m.swap * (pn == 2 and -1 or 1)
    end
    
    if m.parabolax ~= 0 then
        xpos = xpos + m.parabolax * (fYOffset / ARROW_SIZE) * (fYOffset / ARROW_SIZE)
    end
    
    if m.parabolay ~= 0 then
        ypos = ypos + m.parabolay * (fYOffset / ARROW_SIZE) * (fYOffset / ARROW_SIZE)
    end
    if m.inside ~= 0 then
        xpos = xpos +math.sin(0 + (fYOffset*0.004))*(ARROW_SIZE* (iCol % 2 == 0 and 1 or -1) * m.inside*0.5);
    end
    if m.attenuatex ~= 0 then
    local fXOffset =getPropertyFromGroup('strumLineNotes',(pn==2 and iCol+4 or iCol),"x")
   xpos = xpos +m.attenuatex * (fYOffset / ARROW_SIZE) * (fYOffset / ARROW_SIZE) * (fXOffset / ARROW_SIZE);
    end
    if m.attenuatey ~= 0 then
    local fXOffset =getPropertyFromGroup('strumLineNotes',(pn==2 and iCol+4 or iCol),"x")
    ypos = ypos +m.attenuatey * (fYOffset / ARROW_SIZE) * (fYOffset / ARROW_SIZE) * (fXOffset / ARROW_SIZE);
    end
    if m.tornado ~= 0 then
		local iTornadoWidth = 2

		local iStartCol = iCol - iTornadoWidth;
		local iEndCol = iCol + iTornadoWidth;
		iStartCol = math.clamp( iStartCol, 0, 4-1 );
		iEndCol = math.clamp( iEndCol, 0, 4-1 );

		local fMinX = 3.402823466*(10^38)
		local fMaxX = 1.175494351*(10^-38)
		
		-- TODO: Don't index by PlayerNumber.

		for i=iStartCol,iEndCol do
		
			fMinX = math.min( fMinX, getPropertyFromGroup('strumLineNotes',(pn==2 and i+4 or i),"x") );
			fMaxX = math.max( fMaxX, getPropertyFromGroup('strumLineNotes',(pn==2 and i+4 or i),"x") );
	end

		local fRealPixelOffset = getPropertyFromGroup('strumLineNotes',(pn==2 and iCol+4 or iCol),"x")
		local fPositionBetween = scale( fRealPixelOffset, fMinX, fMaxX, -1, 1 );
		local fRads = math.acos( fPositionBetween );
		fRads = fRads + (fYOffset * 6 / screenHeight)
		
		local fAdjustedPixelOffset = scale( math.cos(fRads), -1, 1, fMinX, fMaxX );

		xpos = xpos + (fAdjustedPixelOffset - fRealPixelOffset) * m.tornado
    end
    if m.outside ~= 0 then
        local multIn = iCol%4 <= 1 and -1 or 1
        xpos =xpos+math.pow(math.max((fYOffset*0.01)-1, 0)*m.outside, 2) * multIn
    end
    if m.vibrate ~= 0 then
		xpos = xpos + (math.random() - 0.5) * m.vibrate * 20;
		ypos = ypos + (math.random() - 0.5) * m.vibrate * 20;
    end
    if m.bumpy ~= 0 then
		zpos = zpos + m.bumpy * 40*math.sin((fYOffset+(100.0*m.bumpyoffset))/((m.bumpyperiod*16.0)+16.0))
    end
	if m['bumpy'..iCol] ~= 0 then
	    zpos = zpos + m['bumpy'..iCol] * 40*math.sin((fYOffset+(100.0*m.bumpyoffset))/((m.bumpyperiod*16.0)+16.0))
	end
	if m.tanbumpy ~= 0 then
		zpos = zpos + m.tanbumpy * 40*selectTanType((fYOffset+(100.0*m.tanbumpyoffset))/((m.tanbumpyperiod*16.0)+16.0),m.cosecant)
    end
    if m.attenuatez ~= 0 then
    local fXOffset =getPropertyFromGroup('strumLineNotes',(pn==2 and iCol+4 or iCol),"x")
    zpos = zpos +m.attenuatez * (fYOffset / ARROW_SIZE) * (fYOffset / ARROW_SIZE) * (fXOffset / ARROW_SIZE);
    end
    if m.tornadoz ~= 0 then
		local iTornadoWidth = 2

		local iStartCol = iCol - iTornadoWidth;
		local iEndCol = iCol + iTornadoWidth;
		iStartCol = math.clamp( iStartCol, 0, 4-1 );
		iEndCol = math.clamp( iEndCol, 0, 4-1 );

		local fMinX = 3.402823466*(10^38)
		local fMaxX = 1.175494351*(10^-38)
		
		-- TODO: Don't index by PlayerNumber.

		for i=iStartCol,iEndCol do
		
			fMinX = math.min( fMinX, getPropertyFromGroup('strumLineNotes',(pn==2 and i+4 or i),"x") );
			fMaxX = math.max( fMaxX, getPropertyFromGroup('strumLineNotes',(pn==2 and i+4 or i),"x") );
	end

		local fRealPixelOffset = getPropertyFromGroup('strumLineNotes',(pn==2 and iCol+4 or iCol),"x")
		local fPositionBetween = scale( fRealPixelOffset, fMinX, fMaxX, -1, 1 );
		local fRads = math.acos( fPositionBetween );
		fRads = fRads + (fYOffset * 6 / screenHeight)
		
		local fAdjustedPixelOffset = scale( math.cos(fRads), -1, 1, fMinX, fMaxX );

		zpos = zpos + (fAdjustedPixelOffset - fRealPixelOffset) * m.tornadoz
    end
    if m.parabolaz ~= 0 then
        zpos = zpos + m.parabolaz * (fYOffset / ARROW_SIZE) * (fYOffset / ARROW_SIZE)
    end
	if m.sawtoothz ~= 0 then
		zpos = zpos + (m.sawtoothz*ARROW_SIZE) * ((0.5 / (m.sawtoothzperiod+1) * fYOffset) / ARROW_SIZE - math.floor((0.5 / (m.sawtoothzperiod+1) * fYOffset) / ARROW_SIZE) );
	end

	if m.digitalz ~= 0 then
		zpos = zpos + (m.digitalz * ARROW_SIZE * 0.5) * round((m.digitalzsteps+1) * math.sin(math.pi * (fYOffset + (1.0 * m.digitalzoffset ) ) / (ARROW_SIZE + (m.digitalzperiod * ARROW_SIZE) )) )/(m.digitalzsteps+1);
	end
	if m.tandigitalz ~= 0 then
		zpos = zpos + (m.tandigitalz * ARROW_SIZE * 0.5) * round((m.tandigitalzsteps+1) * selectTanType(math.pi * (fYOffset + (1.0 * m.tandigitalzoffset ) ) / (ARROW_SIZE + (m.tandigitalzperiod * ARROW_SIZE) ),m.cosecant) )/(m.tandigitalzsteps+1);
	end
	if m.squarez ~= 0 then
		local fResult = square( (math.pi * (fYOffset+(1.0*(m.squarezoffset))) / (ARROW_SIZE+(m.squarezperiod*ARROW_SIZE))) );
		zpos = zpos + (m.squarez * ARROW_SIZE * 0.5) * fResult;
	end

	if m.bouncez ~= 0 then
		local fBounceAmt = math.abs( math.sin( ( (fYOffset + (1.0 * (m.bouncezoffset) ) ) / ( 60 + m.bouncezperiod*60) ) ) );
		zpos = zpos + m.bouncez * ARROW_SIZE * 0.5 * fBounceAmt;
	end

	if m.zigzagz ~= 0 then
		local fResult = triangle( (math.pi * (1/(m.zigzagzperiod+1)) * ((fYOffset+(100.0*(m.zigzagzoffset)))/ARROW_SIZE) ) );
	    
		zpos = zpos + (m.zigzagz*ARROW_SIZE/2) * fResult;
	end
	if m.beatz ~= 0 then
			
		local fBeatStrength = m.beatz;
		
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
			
			zpos = zpos + fBeatStrength * fShift
			
		end
	
	end
    return xpos, ypos, rotz, zpos
    
end




defaultPositions = {{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0}}
defaultscale = {{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0}}


-- events

mods,curmod = {},1
perframefunc = {}
event,curevent = {},1
songStarted = false

-- some of mirin functions
function set(t)
	table.insert(t,2,0)
	table.insert(t,3,instant)
	ease(t)
end
function ease(t)
	table.insert(mods,t)
end
function perframe(t)
	table.insert(perframefunc,t)
end
function func(t)
	table.insert(event,t)
end

function TEMPLATE.songStart()
    
    downscroll = false

	for i=0,7 do
        defaultPositions[i+1].x = getPropertyFromGroup("strumLineNotes",i,"x")
                defaultPositions[i+1].y = getPropertyFromGroup("strumLineNotes",i,"y")
        defaultscale[i+1].x = getPropertyFromGroup("strumLineNotes",i,"scale.x")
        defaultscale[i+1].y = getPropertyFromGroup("strumLineNotes",i,"scale.y")
        
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
	luaDebugMode = true

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
		for i = 4, #v, 2 do
		local mn = v[i + 1]
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
				targetMods[pn][mn] = v[i]
				tweenEx1[pn][mn] = v.ex1
				tweenEx2[pn][mn] = v.ex2
				isTweening[pn][mn] = true
			end
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
	if #perframefunc>0 then
		for i=1,#perframefunc do
			local a = perframefunc[i]
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

	setCamNotes(activeMods[1].camx,activeMods[1].camy,activeMods[1].rotationz + activeMods[1].camwag * math.sin(beat*math.pi),activeMods[1].camalpha,activeMods[1].camzoom)
	

			
	if songStarted then
		for pn=1,2 do
			local xmod = activeMods[pn].xmod
			for col=0,3 do
				local c = (pn-1)*4 + col
				local xp, yp, rz, zp = arrowEffects(0, col, pn)
				local alp = receptorAlpha(col,pn)

				--print('Areceptor '..c..' is '..tostring(receptor))
			
				local defaultx, defaulty = defaultPositions[c+1].x, defaultPositions[c+1].y
				    setPropertyFromGroup('strumLineNotes', c, 'x', defaultx + xp)
			setPropertyFromGroup('strumLineNotes', c, 'y', defaulty + yp)
                setPropertyFromGroup("strumLineNotes",c,"angle",rz)
                setPropertyFromGroup("strumLineNotes",c,"alpha",alp)
			
    --[[local rotx,roty,rotz = receptorRotation(0,col,pn)
    local rotation=rotationXYZ(rotx, roty, rotz)
    local anglepos={x=rotation.m00+rotation.m01+rotation.m02+rotation.m03,
    y=rotation.m10+rotation.m11+rotation.m12+rotation.m13,
    z=rotation.m20+rotation.m21+rotation.m22+rotation.m23}
    setPropertyFromGroup('strumLineNotes', c, 'x',getPropertyFromGroup('strumLineNotes', c, 'x')+anglepos.x)
	setPropertyFromGroup('strumLineNotes', c, 'y', getPropertyFromGroup("strumLineNotes",c,"y")+anglepos.y)
			zp=zp+anglepos.z]]
			local scalex, scaley = getScale(0, col, pn, defaultscale[c+1].x, defaultscale[c+1].y)

	local zNear,zFar = 0,100
	local zRange = zNear - zFar
	local fov = 90
	local tanHalfFOV = math.tan(math.rad(fov/2))

			local origin={x=getPropertyFromGroup("strumLineNotes",c,"x") - (screenWidth/2),y=getPropertyFromGroup("strumLineNotes",c,"y") - (screenHeight/2),z=zp}

			local pos={x=origin.x,y=origin.y,z=origin.z/1000-1}
			local X = pos.x*(1/tanHalfFOV)/-pos.z+(screenWidth/2)
			local Y = pos.y/(1/tanHalfFOV)/-pos.z+(screenHeight/2)
			setPropertyFromGroup('strumLineNotes', c, 'x', X)
			setPropertyFromGroup('strumLineNotes', c, 'y', Y)

			local scale = -pos.z
			scale = 1 / scale
setPropertyFromGroup('strumLineNotes', c, 'scale.x', scalex * scale)
			setPropertyFromGroup('strumLineNotes', c, 'scale.y', scaley * scale)
    local zoom = getZoom(0,col,pn)

setPropertyFromGroup("strumLineNotes",c,"scale.x",getPropertyFromGroup("strumLineNotes",c,"scale.x")*zoom)
setPropertyFromGroup("strumLineNotes",c,"scale.y",getPropertyFromGroup("strumLineNotes",c,"scale.y")*zoom)

			
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
				local scrollSpeeds = xmod * activeMods[pn]['xmod'..col] * (1 - 2*getReverseForCol(col,pn)) * scrollSpeed
				
				local off = (1 - 2*getReverseForCol(col,pn))

				local ypos = getYAdjust(defaulty - (getSongPosition() - targTime),col,pn) * scrollSpeeds * 0.45 - off + ARROW_SIZE / 4
					local zoom = getZoom(ypos-defaulty,col,pn)
				local xa, ya, rz, za = arrowEffects(ypos-defaulty, col, pn)
				local alp = arrowAlpha(ypos-defaulty, col, pn)
			    local scalex, scaley = getScale(ypos-defaulty, col, pn, defaultscale[c+1].x, defaultscale[c+1].y)
				if getPropertyFromGroup('notes',v,"isSustainNote") --[[and not note.isParent]] then
					local ypos2 = getYAdjust(defaulty - ((getSongPosition()+.1) - targTime),col,pn) * scrollSpeeds * 0.45 - off + ARROW_SIZE / 2

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
            	setPropertyFromGroup("notes",v,"y",ypos + ya + (getPropertyFromGroup('notes',v,"isSustainNote") and 35 or 0))
            	setPropertyFromGroup("notes",v,"alpha",alp)



    local zNear,zFar = 0,100
	local zRange = zNear - zFar
	local fov = 90
	local tanHalfFOV = math.tan(math.rad(fov/2))
			local origin={x=getPropertyFromGroup('notes',v,"x") - (screenWidth/2),y= getPropertyFromGroup('notes',v,"y") - (screenHeight/2),z=za}
			local pos={x=origin.x,y=origin.y,z=origin.z/1000-1}

			local X = pos.x*(1/tanHalfFOV)/-pos.z+(screenWidth/2)
			local Y = pos.y/(1/tanHalfFOV)/-pos.z+(screenHeight/2)
			setPropertyFromGroup('notes', v, 'x', X)
			setPropertyFromGroup('notes', v, 'y', Y)

			local scale = -pos.z
			scale = 1 / scale
			local scalenewy=getPropertyFromGroup('notes',v,"isSustainNote") and 1 or scaley
			setPropertyFromGroup('notes', v, 'scale.x', scalex * scale)
			setPropertyFromGroup('notes', v, 'scale.y', scalenewy * scale)
			
			setPropertyFromGroup("notes",v,"scale.x",(getPropertyFromGroup('notes',v,"scale.x")*zoom))

  			setPropertyFromGroup("notes",v,"scale.y",(getPropertyFromGroup('notes',v,"scale.y")*(getPropertyFromGroup('notes',v,"isSustainNote") and 1 or zoom)))

  
			end
			
		end
		
	end
	updateCommand(elapsed,beat)
end

function onCreatePost()
	TEMPLATE.InitMods()

	--WRITE MODS HERE! 
	
	initCommand()
	
	--must be called at END of start
	TEMPLATE.setup()
	
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

--callbacks
function initCommand()
	local me = ease
	local m2 = func
	local mpf = perframe
    local function hide(t)
		local bt,tpn = t[1],t.pn
		for i=0,3 do
			me{bt+i*.125-1,.5,outExpo,-70,'movey'..i,pn=tpn}
			me{bt+i*.125-.5,1.25,inExpo,650,'movey'..i,pn=tpn}
			set{bt+i*.125+1.75,1,'stealth',pn=tpn}
			set{bt+i*.125+1.75,1,'dark',pn=tpn}
		end
	end
	local function unhide(t)
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
	local function wig(t)
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
	local function mod_bounce2(t)
		local pn = t.pn
		local beat,length,ease1,ease2,amt,mod,pn = t[1],t[2],t[3],t[4],t[5],t[6],t.plr
		me{beat, (length/2), ease1, amt, mod, plr=pn}
		me{beat+(length/2), (length/2), ease2, -amt, mod, plr=pn}
	end
	me{0,1,linear,5,'pulseouter',1,'drunk'}
end

function updateCommand(elapsed,beat)

end
--end callbacks

return 0
