print('LOADED TEMPLATE')

TEMPLATE = {}

--[[
			READ MEEEEEEEE!
	The template is modified from TaroNuke's TEMPLATE 1,
	now i will keep updating this shit instead of updating the old one,
	only for psych engine (support version >= 0.6.3)
	i prefer to use windows psych engine 0.7.3
	it's like psych modcharting tool lua style
	Credits:
		Original Template : TaroNuke
		Modified : UntilYouAreGone / asdf1234 / 159357159357asdfghjkl
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

function inCubic(t, b, c, d)
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

local sine_table_size = 1024
local sine_index_mod = sine_table_size * 2
local sine_table_index_mult = sine_index_mod / (math.pi * 2)
local sine_table = {}
local table_is_inited = false

function fastSin(x)
	if table_is_inited == false then
		for i = 0,sine_table_size-1 do
			local angle = i * math.pi / sine_table_size
			sine_table[i] = math.sin(angle)
		end
		table_is_inited = true
	end
	if x == 0 then return 0 end
	local index = x * sine_table_index_mult
	while index < 0 do
		index = index + sine_index_mod
	end
	local first_index = math.floor(index)
	local second_index = (first_index + 1) % sine_index_mod
	local remainder = index - first_index
	first_index = first_index % sine_index_mod
	local first = 0.0
	local second = 0.0
	if first_index >= sine_table_size then
		first = -sine_table[first_index - sine_table_size]
	else
		first = sine_table[first_index]
	end
	if second_index >= sine_table_size then
		second = -sine_table[second_index - sine_table_size]
	else
		second = sine_table[second_index]
	end
	local result = remainder * (second - first) + first
	return result
end

function fastCos(x)
	return fastSin(x+math.pi*0.5)
end

function fastTan(x)
	return fastSin(x) / fastCos(x)
end

function fastCot(x)
	return fastCos(x) / fastSin(x)
end

function fastSec(x)
	return 1 / fastCos(x)
end

function fastCsc(x)
	return 1 / fastSin(x)
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

function quantize(f,fRoundInterval)
return math.floor((f + fRoundInterval / 2) / fRoundInterval) * fRoundInterval
end

function RotationXYZ( vec, rX, rY, rZ )
    local PI=math.pi
	rX = rX*(PI/180)
	rY = rY*(PI/180)
	rZ = rZ*(PI/180)

	local cX = fastCos(rX)
	local sX = fastSin(rX)
	local cY = fastCos(rY)
	local sY = fastSin(rY)
	local cZ = fastCos(rZ)
	local sZ = fastSin(rZ)

	local mat = {
	 	m00=cZ*cY, m01=cZ*sY*sX+sZ*cX, m02=cZ*sY*cX+sZ*(-sX), m03=0,
	 	m10=(-sZ)*cY, m11=(-sZ)*sY*sX+cZ*cX, m12=(-sZ)*sY*cX+cZ*(-sX), m13=0,
	 	m20=-sY, m21=cY*sX, m22=cY*cX, m23=0,
	 	m30=0, m31=0, m32=0, m33=1
	  }
	local vec = {
		x = vec.x * mat.m00 + vec.y * mat.m10 + vec.z * mat.m20,
		y = vec.x * mat.m01 + vec.y * mat.m11 + vec.z * mat.m21,
		z = vec.x * mat.m02 + vec.y * mat.m12 + vec.z * mat.m22
	}
	return vec
end

function round(val)
	return math.floor(val+.5)
end

function selectTanType(angle,is_csc,is_sec,is_cot)
if is_csc ~= 0 then
return fastCsc(angle)
elseif is_sec ~= 0 then
return fastSec(angle)
elseif is_cot ~= 0 then
return fastCot(angle)
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
	beatmult = 0,
	beatoffset = 0,
	beatperiod = 0,
	flip = 0,
	invert = 0,
	divide = 0,
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
	beatymult = 0,
	beatyoffset = 0,
	beatyperiod = 0,
	sawtooth = 0,
	sawtoothoffset = 0,
	sawtoothperiod = 0,
	digital = 0,
	digitalsteps = 0,
	digitaloffset = 0,
	digitalperiod = 0,
	tandigital = 0,
	tandigitalsteps = 0,
	tandigitaloffset = 0,
	tandigitalperiod = 0,
	digitaly = 0,
	digitalysteps = 0,
	digitalyoffset = 0,
	digitalyperiod = 0,
	tandigitaly = 0,
	tandigitalysteps = 0,
	tandigitalyoffset = 0,
	tandigitalyperiod = 0,
	square = 0,
	squareoffset = 0,
	squareperiod = 0,
	squarey = 0,
	squareyoffset = 0,
	squareyperiod = 0,
	bounce = 0,
	bounceoffset = 0,
	bounceperiod = 0,
	bouncey = 0,
	bounceyoffset = 0,
	bounceyperiod = 0,
	xmode = 0,
	tiny = 0,
	tinyx = 0,
	tinyy = 0,
	tinyz = 0,
	zigzag = 0,
	zigzagoffset = 0,
	zigzagperiod = 0,
	attenuatex = 0,
	attenuatey = 0,
	attenuatez = 0,
	swap = 0,
	parabolax = 0,
	parabolay = 0,
	tornado = 0,
	tornadooffset = 0,
	tornadoperiod = 0,
	tantornado = 0,
	tantornadooffset = 0,
	tantornadoperiod = 0,
	tantornadoz = 0,
	tantornadozoffset = 0,
	tantornadozperiod = 0,
	tornadoy = 0,
	tornadoyoffset = 0,
	tornadoyperiod = 0,
	zoom = 0,
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
	camx = 0,
	camy = 0,
	camalpha = 1,
	camzoom = 0,
	camangle = 0,
	camwag = 0,
	xmod = 1, --scrollSpeed
	cotangent = 0,
	secant = 0,
	cosecant = 0,
	vibrate = 0,
	vibratez = 0,
	bumpy = 0,
	bumpyoffset = 0,
	bumpyperiod = 0,
	tanbumpy = 0,
	tanbumpyoffset = 0,
	tanbumpyperiod = 0,
	attenuatez = 0,
	tornadoz = 0,
	tornadozoffset = 0,
	tornadozperiod = 0,
	parabolaz = 0,
	sawtoothz = 0,
	sawtoothzoffset = 0,
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
	beatzmult = 0,
	beatzoffset = 0,
	beatzperiod = 0,
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
	confusionoffset = 0,
	rotatex = 0,
	rotatey = 0,
	rotatez = 0,
	rotatetype = 0,
	rotateoriginx = 0, -- optional value
	rotateoriginy = 0, -- optional value
	scalex = 1,
	scaley = 1,
	scalez = 1,
	scale = 1,
	movew = 1,
	amovew = 1,
	waveoffset = 0,
	cmod = -1,
	mmod = 10,
	pingpong = 0,
	pingpongtype = 0,
	zoomx = 1,
	zoomy = 1,
	skew = 0,
	skewx = 0,
	skewy = 0,
	rotationx = 0,
	rotationy = 0,
	rotationz = 0,
	unboundedreverse = 1,
	incominganglex = 0,
	incomingangley = 0,
	incominganglez = 0,
	drawsize = 1,
	stealthtype = 1, -- bool
	stealthpastreceptors = 1, -- bool
	smooth = 0,
	centered = 0,
	infinite = 0,
	glitch = 0,
	hidenoteflash = 0,
	longboy = 0,
}

--column specific mods
for i=0,3 do
	modList['movex'..i] = 0
	modList['movey'..i] = 0
	modList['movez'..i] = 0
	modList['movew'..i] = 1
	modList['amovex'..i] = 0
	modList['amovey'..i] = 0
	modList['amovez'..i] = 0
	modList['amovew'..i] = 1
	modList['dark'..i] = 0
	modList['stealth'..i] = 0
	modList['tinyx'..i] = 0
	modList['tinyy'..i] = 0
	modList['tinyz'..i] = 0
	modList['confusion'..i] = 0
	modList['confusionoffset'..i] = 0
	modList['incominganglex'..i] = 0
	modList['incomingangley'..i] = 0
	modList['incominganglez'..i] = 0
	modList['reverse'..i] = 0
	modList['tiny'..i] = 0
	modList['zoom'..i] = 0
	modList['scalex'..i] = 1
	modList['scaley'..i] = 1
	modList['scalez'..i] = 1
	modList['scale'..i] = 1
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
tweenFlip = {{},{}}
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

local FADE_DIST_Y = 40
function GetHiddenSudden(m)
	return m.hidden * m.sudden
end
function GetHiddenEndLine(m)
	return screenHeight / 2 +
		FADE_DIST_Y * scale( GetHiddenSudden(m), 0., 1., -1.0, -1.25 ) +
		screenHeight / 2 * m.hiddenoffset;
end

function GetHiddenStartLine(m)
	return screenHeight / 2 +
		FADE_DIST_Y * scale( GetHiddenSudden(m), 0., 1., 0.0, -0.25 ) +
		screenHeight / 2 * m.hiddenoffset;
end

function GetSuddenEndLine(m)
	return screenHeight / 2 +
		FADE_DIST_Y * scale( GetHiddenSudden(m), 0., 1., -0.0, 0.25 ) +
		screenHeight / 2 * m.suddenoffset;
end

function GetSuddenStartLine(m)
	return screenHeight / 2 +
		FADE_DIST_Y * scale( GetHiddenSudden(m), 0., 1., 1.0, 1.25 ) +
		screenHeight / 2 * m.suddenoffset;
end

function receptorAlpha(iCol,pn)
	local alp = 1

	local m = activeMods[pn]
	if m.glitch ~= 0 then
		local time = getSongPosition()/1000
		local f = fastSin(time * 50)
		alp = alp + (f-1)
	end
	if m.alpha ~= 1 then
		alp = alp*m.alpha
	end
	if m.dark ~= 0 or m['dark'..iCol] ~= 0 then
		alp = alp*(1-m.dark)*(1-m['dark'..iCol])
	end

	return alp
end

function arrowAlpha(fYOffset, iCol,pn,ypos)
	local alp = 1
	local m = activeMods[pn]
	local fYPos = 0
	if m.stealthtype == 1 then
		fYPos = fYOffset
	else
		fYPos = ypos
	end

	if fYPos < 0 and m.stealthpastreceptors == false then
		return 1
	end

	if m.alpha ~= 1 then
		alp = alp*m.alpha
	end
	if m.stealth ~= 0 or m['stealth'..iCol] ~= 0 then
		alp = alp - m.stealth - m['stealth'..iCol]
	end
	if m.hidden ~= 0 then
		local fHiddenVisibleAdjust = scale( fYPos, GetHiddenStartLine(m), GetHiddenEndLine(m), 0, -1 );
		fHiddenVisibleAdjust = math.clamp( fHiddenVisibleAdjust, -1, 0 );
		alp = alp + m.hidden * fHiddenVisibleAdjust;
	end
	if m.sudden ~= 0 then
		local fSuddenVisibleAdjust = scale( fYPos, GetSuddenStartLine(m), GetSuddenEndLine(m), -1, 0 );
		fSuddenVisibleAdjust = math.clamp( fSuddenVisibleAdjust, -1, 0 );
		alp = alp + m.sudden * fSuddenVisibleAdjust;
	end
	if m.blink ~= 0 then
		local time = getSongPosition()/1000
		local f = fastSin(time*10)
		f=quantize(f,0.3333)
		alp = alp + scale( f, 0, 1, -1, 0 );
	end

	if m.randomvanish ~= 0 then
		local fRealFadeDist = 80;
		alp = alp + scale( math.abs(fYPos-360), fRealFadeDist, 2*fRealFadeDist, -1, 0 ) * m.randomvanish;
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
	if not m.unboundedreverse ~= 0 then
	if val > 2 then val = val % 2 end
	if val > 1 then val = scale(val, 1, 2, 1, 0) end
	end
	return val
end

function getYAdjust(fYOffset, iCol, pn)

	local m = activeMods[pn]

	local yadj = 0
	local fScrollSpeed = 1
	if m.wave ~= 0 then
		yadj =yadj + m.wave * 20 *fastSin( (m.waveoffset+fYOffset)/((m.waveperiod*38)+38) );
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
	    local fExpandMultiplier = scale(fastCos(expandSeconds * 3 * (m.expandperiod + 1)), -1, 1, 0.75, 1.75);
      fScrollSpeed = fScrollSpeed * scale(m.expand, 0, 1, 1, fExpandMultiplier);
	end
	if m.tanexpand ~= 0 then
	local last = 0
	local time = getSongPosition() / 1000
	local tanExpandSeconds = 0
    tanExpandSeconds = tanExpandSeconds + (time - last);
    tanExpandSeconds = tanExpandSeconds % ((math.pi * 2) / (m.tanexpandperiod + 1));
    last = time
	    local fTanExpandMultiplier = scale(selectTanType(tanExpandSeconds * 3 * (m.tanexpandperiod + 1),m.cosecant,m.secant,m.cotangent), -1, 1, 0.75, 1.75);
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
    local sine = fastSin(((fYOffset+(100.0*m.pulseoffset))/(0.4*(ARROW_SIZE+(m.pulseperiod*ARROW_SIZE)))))

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

function getScale(fYOffset, iCol, pn, sx, sy, isNote)
    local x = sx
    local y = sy
	local z = 0
    local m = activeMods[pn]
	x = x * m.scalex * m['scalex'..iCol] * m.scale * m['scale'..iCol]
    y = y * m.scaley * m['scaley'..iCol] * m.scale * m['scale'..iCol]
	z = z * m.scalez * m['scalez'..iCol] * m.scale * m['scale'..iCol]

    local stretch = m.stretch + m['stretch'..iCol]
    local squish = m.squish + m['squish'..iCol]
    local stretchX = lerp(1, 0.5, stretch);
    local stretchY = lerp(1, 2, stretch);

    local squishX = lerp(1, 2, squish);
    local squishY = lerp(1, 0.5, squish);
	x = x * math.pow( 0.5, m.tinyx ) * math.pow( 0.5, m['tinyx'..iCol] )
	y = y * math.pow( 0.5, m.tinyy ) * math.pow( 0.5, m['tinyy'..iCol] )
	z = z * math.pow( 0.5, m.tinyz ) * math.pow( 0.5, m['tinyz'..iCol] )

	x = x * squishX
	x = x * stretchX

	y = y * stretchY
	y = y * squishY
		return x,y
end

function cameraEffects(pn)
    local m = activeMods[pn]
	local xpos, ypos, rotz, alpha, zoom, zoomx, zoomy, skewx, skewy = 0, 0, 0, 0, 0, 0, 0, 0, 0
	xpos = xpos + activeMods[1].camx
	ypos = ypos + activeMods[1].camy
	rotz = rotz + activeMods[1].camangle + activeMods[1].camwag * fastSin(beat*math.pi)
	alpha = alpha + activeMods[1].camalpha
	zoom = zoom + activeMods[1].camzoom
	zoomx = zoomx + activeMods[1].zoomx
	zoomy = zoomy + activeMods[1].zoomy
	skewx = skewx + activeMods[1].skewx + activeMods[1].skew
	skewy = skewy + activeMods[1].skewy + activeMods[1].skew
	return xpos, ypos, rotz, alpha, zoom, zoomx, zoomy, skewx, skewy
end

function arrowEffects(fYOffset, iCol, pn, withreverse)
    local m = activeMods[pn]

    local xpos, ypos, rotz, zpos = 0, 0, 0, 0

	if m['confusion'..iCol] ~= 0 or m.confusion ~= 0 or m.confusionoffset ~= 0 or m['confusionoffset'..iCol]~= 0 then
		rotz = rotz + m['confusion'..iCol] + m.confusion + m['confusionoffset'..iCol] + m.confusionoffset
	end
	if m.dizzy ~= 0 then
		rotz = rotz + m.dizzy*fYOffset
	end
    if m.drunk ~= 0 then
        xpos = xpos + m.drunk * fastCos(getSongPosition()*0.001 * (1 + m.drunkspeed) + iCol * ((m.drunkoffset * 0.2) + 0.2) + fYOffset * ((m.drunkperiod * 10) + 10) / screenHeight) * ARROW_SIZE * 0.5;
    end
    if m.drunkz ~= 0 then
        zpos = zpos + m.drunkz * fastCos(getSongPosition()*0.001 * (1 + m.drunkzspeed) + iCol * ((m.drunkzoffset * 0.2) + 0.2) + fYOffset * ((m.drunkzperiod * 10) + 10) / screenHeight) * ARROW_SIZE * 0.5;
    end
    if m.drunky ~= 0 then
        ypos = ypos + m.drunky * fastCos(getSongPosition()*0.001 * (1 + m.drunkyspeed) + iCol * ((m.drunkyoffset * 0.2) + 0.2) + fYOffset * ((m.drunkyperiod * 10) + 10) / screenHeight) * ARROW_SIZE * 0.5;
    end
    if m.tandrunk ~= 0 then
    xpos = xpos + m.tandrunk * selectTanType(getSongPosition()*0.001 * (1 + m.tandrunkspeed) + iCol * ((m.tandrunkoffset * 0.2) + 0.2) + fYOffset * ((m.tandrunkperiod * 10) + 10) / screenHeight,m.cosecant,m.secant,m.cotangent) * ARROW_SIZE * 0.5;
    end
    if m.tandrunkz ~= 0 then
    zpos = zpos + m.tandrunkz * selectTanType(getSongPosition()*0.001 * (1 + m.tandrunkzspeed) + iCol * ((m.tandrunkzoffset * 0.2) + 0.2) + fYOffset * ((m.tandrunkzperiod * 10) + 10) / screenHeight,m.cosecant,m.secant,m.cotangent) * ARROW_SIZE * 0.5;
    end
    if m.tandrunky ~= 0 then
    ypos = ypos + m.tandrunky * selectTanType(getSongPosition()*0.001 * (1 + m.tandrunkyspeed) + iCol * ((m.tandrunkyoffset * 0.2) + 0.2) + fYOffset * ((m.tandrunkyperiod * 10) + 10) / screenHeight,m.cosecant,m.secant,m.cotangent) * ARROW_SIZE * 0.5;
    end
    if m.tipsy ~= 0 then
        ypos = ypos + m.tipsy * fastCos( getSongPosition() * 0.001 * ((m.tipsyspeed * 1.2) + 1.2) + (iCol * ((m.tipsyoffset * 1.8) + 1.8)))* ARROW_SIZE * 0.4
    end
    if m.tantipsy ~= 0 then
        ypos = ypos + m.tantipsy * selectTanType( getSongPosition() * 0.001 * ((m.tantipsyspeed * 1.2) + 1.2) + (iCol * ((m.tantipsyoffset * 1.8) + 1.8)),m.cosecant,m.secant,m.cotangent)* ARROW_SIZE * 0.4
    end
    if m.tipsyz ~= 0 then
        zpos = zpos + m.tipsyz * fastCos( getSongPosition() * 0.001 * ((m.tipsyzspeed * 1.2) + 1.2) + (iCol * ((m.tipsyzoffset * 1.8) + 1.8)))* ARROW_SIZE * 0.4
    end
    if m.tantipsyz ~= 0 then
        zpos = zpos + m.tantipsyz * selectTanType( getSongPosition() * 0.001 * ((m.tantipsyzspeed * 1.2) + 1.2) + (iCol * ((m.tantipsyzoffset * 1.8) + 1.8)),m.cosecant,m.secant,m.cotangent)* ARROW_SIZE * 0.4
    end
    if m.tipsyx ~= 0 then
        xpos = xpos + m.tipsyx * fastCos( getSongPosition() * 0.001 * ((m.tipsyxspeed * 1.2) + 1.2) + (iCol * ((m.tipsyxoffset * 1.8) + 1.8))) * ARROW_SIZE * 0.4
    end
     if m.tantipsyx ~= 0 then
        xpos = xpos + m.tantipsyx * selectTanType( getSongPosition() * 0.001 * ((m.tantipsyxspeed * 1.2) + 1.2) + (iCol * ((m.tantipsyxoffset * 1.8) + 1.8)),m.cosecant,m.secant,m.cotangent)* ARROW_SIZE * 0.4
    end
    if m.adrunk ~= 0 then
        xpos = xpos + m.adrunk * ( fastCos( getSongPosition()*0.001 + iCol*(0.2) + 1*(0.2) + fYOffset*(10)/720) * ARROW_SIZE*0.5 )
    end
    if m.atipsy ~= 0 then
        ypos = ypos + m.atipsy * ( fastCos( getSongPosition()*0.001 *(1.2) + iCol*(2.0) + 1*(0.2) ) * ARROW_SIZE*0.4 )
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

	if withreverse then
	if m['reverse'..iCol] ~= 0 or m.reverse ~= 0 or m.split ~= 0 or m.cross ~= 0 or m.alternate ~= 0 then
		local shift = getReverseForCol(iCol,pn) * 520
		shift = scale(m.centered, 0., 1., shift, 520 / 2)
		ypos = ypos + shift
	end
	end
	if m.flip ~= 0 then
		local fDistance = ARROW_SIZE * 2 * (1.5 - iCol);
		xpos = xpos + fDistance * m.flip;
	end

	if m.invert ~= 0 then
		local fDistance = ARROW_SIZE * (iCol%2 == 0 and 1 or -1);
		xpos = xpos + fDistance * m.invert;
	end

	if m.divide ~= 0 then
		xpos = xpos + (iCol >= 2 and 1 or -1) * ARROW_SIZE * m.divide
	end

	if m.beat ~= 0 then

		local fBeatStrength = m.beat;

		local fAccelTime = 0.3;
		local fTotalTime = 0.7;

		-- If the song is really fast, slow down the rate, but speed up the
		-- acceleration to compensate or it'll look weird.
		fBeat = (beat + fAccelTime + m.beatoffset) * (1 + m.beatmult);

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

			local fShift = 40.0*fAmount*fastSin( ((fYOffset/(m.beatperiod*30.0+30.0))) + (math.pi/2) );

			xpos = xpos + fBeatStrength * fShift

		end

	end

	if m.beaty ~= 0 then

		local fBeatStrength = m.beaty;

		local fAccelTime = 0.3;
		local fTotalTime = 0.7;

		-- If the song is really fast, slow down the rate, but speed up the
		-- acceleration to compensate or it'll look weird.
		fBeat = (beat + fAccelTime + m.beatyoffset) * (1 + m.beatymult);

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

			local fShift = 40.0*fAmount*fastSin( ((fYOffset/(m.beatyperiod*30.0+30.0))) + (math.pi/2) );

			ypos = ypos + fBeatStrength * fShift

		end

	end

	if m.sawtooth ~= 0 then
		xpos = xpos + (m.sawtooth*ARROW_SIZE) * ((0.5 / (m.sawtoothperiod+1) * (fYOffset + 100*m.sawtoothoffset)) / ARROW_SIZE - math.floor((0.5 / (m.sawtoothperiod+1) * (fYOffset + 100*m.sawtoothoffset)) / ARROW_SIZE) );
	end

	if m.digital ~= 0 then
		xpos = xpos + (m.digital * ARROW_SIZE * 0.5) * round((m.digitalsteps+1) * fastSin(math.pi * (fYOffset + (1.0 * m.digitaloffset ) ) / (ARROW_SIZE + (m.digitalperiod * ARROW_SIZE) )) )/(m.digitalsteps+1);
	end
	if m.tandigital ~= 0 then
		xpos = xpos + (m.tandigital * ARROW_SIZE * 0.5) * round((m.tandigitalsteps+1) * selectTanType(math.pi * (fYOffset + (1.0 * m.tandigitaloffset ) ) / (ARROW_SIZE + (m.tandigitalperiod * ARROW_SIZE) ),m.cosecant,m.secant,m.cotangent) )/(m.tandigitalsteps+1);
	end
	if m.digitaly ~= 0 then
		ypos = ypos + (m.digitaly * ARROW_SIZE * 0.5) * round((m.digitalysteps+1) * fastSin(math.pi * (fYOffset + (1.0 * m.digitalyoffset ) ) / (ARROW_SIZE + (m.digitalyperiod * ARROW_SIZE) )) )/(m.digitalysteps+1);
	end
	if m.tandigitaly ~= 0 then
		ypos = ypos + (m.tandigitaly * ARROW_SIZE * 0.5) * round((m.tandigitalysteps+1) * selectTanType(math.pi * (fYOffset + (1.0 * m.tandigitalyoffset ) ) / (ARROW_SIZE + (m.tandigitalyperiod * ARROW_SIZE) ),m.cosecant,m.secant,m.cotangent) )/(m.tandigitalysteps+1);
	end
	if m.bumpyx ~= 0 then
		xpos = xpos + m.bumpyx * 40*fastSin((fYOffset+(100.0*m.bumpyxoffset))/((m.bumpyxperiod*16.0)+16.0));
	end

	if m.square ~= 0 then
		local fResult = square( (math.pi * (fYOffset+(1.0*(m.squareoffset))) / (ARROW_SIZE+(m.squareperiod*ARROW_SIZE))) );
		xpos = xpos + (m.square * ARROW_SIZE * 0.5) * fResult;
	end
	if m.squarey ~= 0 then
		local fResult = square( (math.pi * (fYOffset+(1.0*(m.squareyoffset))) / (ARROW_SIZE+(m.squareyperiod*ARROW_SIZE))) );
		ypos = ypos + (m.squarey * ARROW_SIZE * 0.5) * fResult;
	end
    if m.bumpyy ~= 0 then
		ypos = ypos + m.bumpyy * 40*fastSin((fYOffset+(100.0*m.bumpyyoffset))/((m.bumpyyperiod*16.0)+16.0));
	end
	if m.bounce ~= 0 then
		local fBounceAmt = math.abs( fastSin( ( (fYOffset + (1.0 * (m.bounceoffset) ) ) / ( 60 + m.bounceperiod*60) ) ) );
		xpos = xpos + m.bounce * ARROW_SIZE * 0.5 * fBounceAmt;
	end

	if m.bouncey ~= 0 then
		local fBounceAmt = math.abs( fastSin( ( (fYOffset + (1.0 * (m.bounceyoffset) ) ) / ( 60 + m.bounceyperiod*60) ) ) );
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
        xpos = xpos +fastSin(0 + (fYOffset*0.004))*(ARROW_SIZE* (iCol % 2 == 0 and 1 or -1) * m.inside*0.5);
    end
    if m.attenuatex ~= 0 then
    local fXOffset = {{-ARROW_SIZE*2,-ARROW_SIZE,ARROW_SIZE,ARROW_SIZE*2},{-ARROW_SIZE*2,-ARROW_SIZE,ARROW_SIZE,ARROW_SIZE*2}}
   xpos = xpos +m.attenuatex * (fYOffset / ARROW_SIZE) * (fYOffset / ARROW_SIZE) * (fXOffset[pn][iCol+1] / ARROW_SIZE);
    end
    if m.attenuatey ~= 0 then
    local fXOffset = {{-ARROW_SIZE*2,-ARROW_SIZE,ARROW_SIZE,ARROW_SIZE*2},{-ARROW_SIZE*2,-ARROW_SIZE,ARROW_SIZE,ARROW_SIZE*2}}
    ypos = ypos +m.attenuatey * (fYOffset / ARROW_SIZE) * (fYOffset / ARROW_SIZE) * (fXOffset[pn][iCol+1] / ARROW_SIZE);
    end
	if m.attenuatez ~= 0 then
		local fXOffset = {{-ARROW_SIZE*2,-ARROW_SIZE,ARROW_SIZE,ARROW_SIZE*2},{-ARROW_SIZE*2,-ARROW_SIZE,ARROW_SIZE,ARROW_SIZE*2}}
		zpos = zpos + m.attenuatez * (fYOffset / ARROW_SIZE) * (fYOffset / ARROW_SIZE) * (fXOffset[pn][iCol+1] / ARROW_SIZE);
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
		local xoff = {{92,204,316,428},{732,844,956,1068}}
		for i=iStartCol,iEndCol do

			fMinX = math.min( fMinX, xoff[pn][i+1] );
			fMaxX = math.max( fMaxX, xoff[pn][i+1] );
	end

		local fRealPixelOffset = xoff[pn][iCol+1]
		local fPositionBetween = scale( fRealPixelOffset, fMinX, fMaxX, -1, 1 );
		local fRads = math.acos( fPositionBetween );
		fRads = fRads + ((fYOffset + m.tornadooffset) * ((6 * m.tornadoperiod) + 6) / screenHeight)

		local fAdjustedPixelOffset = scale( fastCos(fRads), -1, 1, fMinX, fMaxX );

		xpos = xpos + (fAdjustedPixelOffset - fRealPixelOffset) * m.tornado
    end
	if m.tornadoy ~= 0 then
		local iTornadoWidth = 2

		local iStartCol = iCol - iTornadoWidth;
		local iEndCol = iCol + iTornadoWidth;
		iStartCol = math.clamp( iStartCol, 0, 4-1 );
		iEndCol = math.clamp( iEndCol, 0, 4-1 );

		local fMinX = 3.402823466*(10^38)
		local fMaxX = 1.175494351*(10^-38)

		-- TODO: Don't index by PlayerNumber.
		local xoff = {{92,204,316,428},{732,844,956,1068}}
		for i=iStartCol,iEndCol do

			fMinX = math.min( fMinX, xoff[pn][i+1] );
			fMaxX = math.max( fMaxX, xoff[pn][i+1] );
	end

		local fRealPixelOffset = xoff[pn][iCol+1]
		local fPositionBetween = scale( fRealPixelOffset, fMinX, fMaxX, -1, 1 );
		local fRads = math.acos( fPositionBetween );
		fRads = fRads + ((fYOffset + m.tornadoyoffset) * ((6 * m.tornadoyperiod) + 6) / screenHeight)

		local fAdjustedPixelOffset = scale( fastCos(fRads), -1, 1, fMinX, fMaxX );

		ypos = ypos + (fAdjustedPixelOffset - fRealPixelOffset) * m.tornadoy
    end
    if m.tantornado ~= 0 then
		local xoff = {{92,204,316,428},{732,844,956,1068}}
		local iTornadoWidth = 2

		local iStartCol = iCol - iTornadoWidth;
		local iEndCol = iCol + iTornadoWidth;
		iStartCol = math.clamp( iStartCol, 0, 4-1 );
		iEndCol = math.clamp( iEndCol, 0, 4-1 );

		local fMinX = 3.402823466*(10^38)
		local fMaxX = 1.175494351*(10^-38)

		-- TODO: Don't index by PlayerNumber.

		for i=iStartCol,iEndCol do

			fMinX = math.min( fMinX, xoff[pn][i+1] );
			fMaxX = math.max( fMaxX, xoff[pn][i+1] );
	end

		local fRealPixelOffset = xoff[pn][iCol+1]
		local fPositionBetween = scale( fRealPixelOffset, fMinX, fMaxX, -1, 1 );
		local fRads = math.acos( fPositionBetween );
		fRads = fRads + ((fYOffset + m.tantornadooffset) * ((6 * m.tantornadoperiod) + 6) / screenHeight)

		local fAdjustedPixelOffset = scale( selectTanType(fRads,m.cosecant,m.secant,m.cotangent), -1, 1, fMinX, fMaxX );

		xpos = xpos + (fAdjustedPixelOffset - fRealPixelOffset) * m.tantornado
    end
	if m.tantornadoz ~= 0 then
		local xoff = {{92,204,316,428},{732,844,956,1068}}
		local iTornadoWidth = 2

		local iStartCol = iCol - iTornadoWidth;
		local iEndCol = iCol + iTornadoWidth;
		iStartCol = math.clamp( iStartCol, 0, 4-1 );
		iEndCol = math.clamp( iEndCol, 0, 4-1 );

		local fMinX = 3.402823466*(10^38)
		local fMaxX = 1.175494351*(10^-38)

		-- TODO: Don't index by PlayerNumber.

		for i=iStartCol,iEndCol do

			fMinX = math.min( fMinX, xoff[pn][i+1] );
			fMaxX = math.max( fMaxX, xoff[pn][i+1] );
	end

		local fRealPixelOffset = xoff[pn][iCol+1]
		local fPositionBetween = scale( fRealPixelOffset, fMinX, fMaxX, -1, 1 );
		local fRads = math.acos( fPositionBetween );
		fRads = fRads + ((fYOffset + m.tantornadozoffset) * ((6 * m.tantornadozperiod) + 6) / screenHeight)

		local fAdjustedPixelOffset = scale( selectTanType(fRads,m.cosecant,m.secant,m.cotangent), -1, 1, fMinX, fMaxX );

		zpos =zpos + (fAdjustedPixelOffset - fRealPixelOffset) * m.tantornadoz
    end

    if m.vibrate ~= 0 then
		xpos = xpos + (math.random() - 0.5) * m.vibrate * 20;
		ypos = ypos + (math.random() - 0.5) * m.vibrate * 20;
    end
	if m.vibratez ~= 0 then
		zpos = zpos + (math.random() - 0.5) * m.vibratez * 20;
    end
    if m.bumpy ~= 0 then
		zpos = zpos + m.bumpy * 40*fastSin((fYOffset+(100.0*m.bumpyoffset))/((m.bumpyperiod*16.0)+16.0))
    end
	if m['bumpy'..iCol] ~= 0 then
	    zpos = zpos + m['bumpy'..iCol] * 40*fastSin((fYOffset+(100.0*m.bumpyoffset))/((m.bumpyperiod*16.0)+16.0))
	end
	if m.tanbumpy ~= 0 then
		zpos = zpos + m.tanbumpy * 40*selectTanType((fYOffset+(100.0*m.tanbumpyoffset))/((m.tanbumpyperiod*16.0)+16.0),m.cosecant,m.secant,m.cotangent)
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
		local xoff = {{92,204,316,428},{732,844,956,1068}}
		for i=iStartCol,iEndCol do

			fMinX = math.min( fMinX, xoff[pn][i+1] );
			fMaxX = math.max( fMaxX, xoff[pn][i+1] );
	end

		local fRealPixelOffset = xoff[pn][iCol+1]
		local fPositionBetween = scale( fRealPixelOffset, fMinX, fMaxX, -1, 1 );
		local fRads = math.acos( fPositionBetween );
		fRads = fRads + ((fYOffset + m.tornadozoffset) * ((6 * m.tornadozperiod) + 6) / screenHeight)

		local fAdjustedPixelOffset = scale( fastCos(fRads), -1, 1, fMinX, fMaxX );

		zpos = zpos + (fAdjustedPixelOffset - fRealPixelOffset) * m.tornadoz
    end
    if m.parabolaz ~= 0 then
        zpos = zpos + m.parabolaz * (fYOffset / ARROW_SIZE) * (fYOffset / ARROW_SIZE)
    end
	if m.sawtoothz ~= 0 then
		zpos = zpos + (m.sawtoothz*ARROW_SIZE) * ((0.5 / (m.sawtoothzperiod+1) * (fYOffset + 100*m.sawtoothzoffset)) / ARROW_SIZE - math.floor((0.5 / (m.sawtoothzperiod+1) * (fYOffset + 100*m.sawtoothzoffset)) / ARROW_SIZE) );
	end

	if m.digitalz ~= 0 then
		zpos = zpos + (m.digitalz * ARROW_SIZE * 0.5) * round((m.digitalzsteps+1) * fastSin(math.pi * (fYOffset + (1.0 * m.digitalzoffset ) ) / (ARROW_SIZE + (m.digitalzperiod * ARROW_SIZE) )) )/(m.digitalzsteps+1);
	end
	if m.tandigitalz ~= 0 then
		zpos = zpos + (m.tandigitalz * ARROW_SIZE * 0.5) * round((m.tandigitalzsteps+1) * selectTanType(math.pi * (fYOffset + (1.0 * m.tandigitalzoffset ) ) / (ARROW_SIZE + (m.tandigitalzperiod * ARROW_SIZE) ),m.cosecant,m.secant,m.cotangent) )/(m.tandigitalzsteps+1);
	end
	if m.squarez ~= 0 then
		local fResult = square( (math.pi * (fYOffset+(1.0*(m.squarezoffset))) / (ARROW_SIZE+(m.squarezperiod*ARROW_SIZE))) );
		zpos = zpos + (m.squarez * ARROW_SIZE * 0.5) * fResult;
	end

	if m.bouncez ~= 0 then
		local fBounceAmt = math.abs( fastSin( ( (fYOffset + (1.0 * (m.bouncezoffset) ) ) / ( 60 + m.bouncezperiod*60) ) ) );
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
		fBeat = (beat + fAccelTime + m.beatzoffset) * (1 + m.beatzmult);

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

			local fShift = 40.0*fAmount*fastSin( ((fYOffset/(m.beatzperiod*30.0+30.0))) + (math.pi/2) );

			zpos = zpos + fBeatStrength * fShift

		end

	end
	if m.smooth ~= 0 then
		ypos = ypos + fastSin(getSongPosition() / 1000 * 2 + iCol*0.8) * ARROW_SIZE / 2 * m.smooth
		xpos = xpos + fastCos(getSongPosition() / 1000 * 3 + iCol*0.8) * ARROW_SIZE / 2 * m.smooth
	end
    return xpos, ypos, rotz, zpos

end

defaultPositions = {{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0}}
defaultscale = {{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0}}


-- events

mods,curmod = {},1
perframe = {}
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
function add(t)
	t.relative = true
	ease(t)
end
function acc(t)
	t.relative = true
	table.insert(t,2,0)
	table.insert(t,3,instant)
	ease(t)
end
function func(t)
	if type(t[2]) == 'number' then
		table.insert(perframe,t)
	else
		table.insert(event,t)
	end
end


function getInfinite()
	local infPath = {{}, {}, {}, {}};
	local r = 0;

	for r=0,360,15 do
		for data = 1,#infPath do
			local rad = r * math.pi / 180;
			table.insert(infPath[data],{x=screenWidth* 0.5 + (fastSin(rad)) * 600,
				y=screenHeight* 0.5 + (fastSin(rad) * fastCos(rad)) * 600, z=0})
		end
	end
	return infPath,1850
end

function createPath(pathfunc)
	local pathData = {}
	local totalDists = {}
		local moveSpeed = 5000
		local path,speed = pathfunc()
		moveSpeed = speed
    local dir = 1
    while dir <= table.getn(path) do
      local idx = 1
      totalDists[dir]=0
      pathData[dir]={}
      while idx<=table.getn(path[dir]) do
        local pos = path[dir][idx]
        if idx ~= 1 then
          local last = pathData[dir][idx-1]
					local _x=pos.x - last.position.x
					local _y=pos.y - last.position.y
					local _z=pos.z - last.position.z
          totalDists[dir] = totalDists[dir] + math.abs(_x*_x+_y*_y+_z*_z)
          local totalDist = totalDists[dir]
          last.End = totalDist
          last.dist = last.start - totalDist
				end
        table.insert(pathData[dir], {
          position = {x=pos.x-ARROW_SIZE/2,y=pos.y-ARROW_SIZE/2,z=pos.z},
          start = totalDists[dir],
          End = 0,
          dist = 0
        })
        idx =idx+1
      end
      dir = dir+1
    end
		return moveSpeed,pathData,totalDists
end

function updatePath(timeDiff,pos,data,mod,pathfunc)
		local moveSpeed,pathData,totalDists = createPath(pathfunc)
    local vDiff = -timeDiff
    local progress  = (vDiff / -moveSpeed) * totalDists[data+1];
		local outPos = {}
		for k, v in pairs(pos) do
			outPos[k] = v
		end
		local function lerpV(p,g,a)
			return{
				x=a*g.x+p.x*(1-a),
				y=a*g.y+p.y*(1-a),
				z=a*g.z+p.z*(1-a),
			}
		end
    local daPath = pathData[data+1];
    if progress<=0 then return lerpV(pos,daPath[1].position,mod) end
    local idx=1
    while idx<=table.getn(daPath) do
      local cData = daPath[idx];
      local nData = daPath[idx+1];
      if nData and cData then
        if progress>cData.start and progress<cData.End then
          local alpha = (cData.start - progress)/cData.dist;
          local interpPos = lerpV(cData.position,nData.position,alpha);
          outPos = lerpV(pos,interpPos,mod);
				end
      end
      idx=idx+1
    end
    return outPos
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
	setProperty('spawnTime',activeMods[1].drawsize * 2500)
	if bTestMode then
	luaDebugMode = true
	end

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
		mn:lower()
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
				if not v.flip then v.flip = false end
				tweenFlip[pn][mn] = v.flip
				if v.startVal then
					storedMods[pn][mn] = v.startVal
				else
					storedMods[pn][mn] = activeMods[pn][mn]
				end
				if v.relative then
				targetMods[pn][mn] = targetMods[pn][mn] + v[i]
				else
				targetMods[pn][mn] = v[i]
				end
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
				local str = curve(curtime, startstrength, diff, duration, tweenEx1[pn][v], tweenEx2[pn][v])
				local strength = tweenFlip[pn][v] and 1 - str or str
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
				a[3](beat, activeMods)
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

	local camx, camy, camrotz, camalpha, camzoom, camzoomx, camzoomy, camskewx, camskewy = cameraEffects(1)
	runHaxeCode([[
		var cs:FlxCamera = getVar("camNotes");
		cs.angle = ]]..camrotz..[[;
		cs.alpha = ]]..camalpha..[[;
		cs.zoom = game.camHUD.zoom + ]]..camzoom..[[;
		cs.scaleX = ]]..camzoomx..[[;
		cs.scaleY = ]]..camzoomy..[[;
		cs.x = ]]..camx..[[;
		cs.y = ]]..camy..[[;
	]])

	-- use matrix
	setProperty('camNotes.canvas.__transform.c',camskewx)
	setProperty('camNotes.canvas.__transform.b',camskewy)


	if songStarted then

		for pn=1,2 do
			local xmod = activeMods[pn].xmod
			for col=0,3 do
				local c = (pn-1)*4 + col
				local xp, yp, rz, zp = arrowEffects(0, col, pn, true)
				local xp2, yp2 = arrowEffects(.1, col, pn, true)
				local alp = receptorAlpha(col,pn)
				--print('Areceptor '..c..' is '..tostring(receptor))
				local defaultx, defaulty = defaultPositions[c+1].x, defaultPositions[c+1].y
				local alpha = getPropertyFromGroup('strumLineNotes',c,'alpha')*math.clamp(scale(alp, 0.5, 0, 1, 0),0,1)
				setPropertyFromGroup('strumLineNotes',c,'colorTransform.alphaMultiplier',alpha)

	local m = activeMods[pn]
	local altx,alty = defaultx+xp,defaulty+yp
		if m.rotatex ~= 0 or m.rotatey ~= 0 or m.rotatez ~= 0 then
			local useX = altx
			local useY = alty
			local originPos = {x=screenWidth/2,y=screenHeight/2}
				if m.rotatetype == 1 then
					originPos.x = defaultx
				elseif m.rotatetype ~= 0 then
					originPos.x = m.rotateoriginx
					originPos.y = m.rotateoriginy
				end
			local vec = {
				x = useX - originPos.x,
				y = useY - originPos.y,
				z = zp
			}
			local rotatedpos = RotationXYZ(vec,m.rotatex,m.rotatey,m.rotatez)
			rotatedpos.x = rotatedpos.x + originPos.x
			rotatedpos.y = rotatedpos.y + originPos.y
			altx= rotatedpos.x
			alty= rotatedpos.y
			zp = rotatedpos.z
		end
		if m.rotationx ~= 0 or m.rotationy ~= 0 or m.rotationz ~= 0 then
			local useX = altx
			local useY = alty
			local originPos = {x=screenWidth/2-ARROW_SIZE/2,y=screenHeight/2}
			local xs = {92,204,316,428}
			local cx = (xs[4] + ARROW_SIZE - xs[1]) / 2 + xs[1]
				originPos.x = originPos.x + (pn*2-3)*cx
			local vec = {
				x = useX - originPos.x,
				y = useY - originPos.y,
				z = zp
			}
			local rotatedpos = RotationXYZ(vec,m.rotationx,m.rotationy,m.rotationz)
			rotatedpos.x = rotatedpos.x + originPos.x
			rotatedpos.y = rotatedpos.y + originPos.y
			altx = rotatedpos.x
			alty = rotatedpos.y
			zp = rotatedpos.z
		end
			setPropertyFromGroup('strumLineNotes', c, 'x', altx)
			setPropertyFromGroup('strumLineNotes', c, 'y', alty)

		setPropertyFromGroup("strumLineNotes",c,"angle",rz)
			local scalex, scaley = getScale(0, col, pn, defaultscale[c+1].x, defaultscale[c+1].y, false)

	local fov = 90
	local tanHalfFOV = math.tan(math.rad(fov/2))
		local updatex,updatey = getPropertyFromGroup('strumLineNotes',c,'x'),getPropertyFromGroup('strumLineNotes',c,'y')
	local fakew = m.movew*m['movew'..col]*m.amovew*m['amovew'..col]
			local origin={x=updatex - (screenWidth/2),y=updatey - (screenHeight/2),z=zp}
			local pos={x=origin.x,y=origin.y,z=(origin.z)/1000-fakew}
			local X = pos.x/(1/tanHalfFOV)/-pos.z+(screenWidth/2)
			local Y = pos.y/(1/tanHalfFOV)/-pos.z+(screenHeight/2)
			setPropertyFromGroup('strumLineNotes', c, 'x', X)
			setPropertyFromGroup('strumLineNotes', c, 'y', Y)
			setPropertyFromGroup('strumLineNotes', c, 'scale.x', scalex / -pos.z)
			setPropertyFromGroup('strumLineNotes', c, 'scale.y', scaley / -pos.z)

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
				local cmod = activeMods[pn].cmod
				local speed = xmod * activeMods[pn]['xmod'..col] * (cmod > 0 and cmod or 1)
				if speed > activeMods[pn].mmod then
					speed = activeMods[pn].mmod
				end
				local scrollSpeeds = speed * (1 - 2*getReverseForCol(col,pn)) * scrollSpeed
				local off = (1 - 2*getReverseForCol(col,pn))
				local ypos = getYAdjust((getSongPosition() - targTime)*(downscroll and 1 or -1),col,pn) * scrollSpeeds * 0.45 - off + defaulty
				local zoom = getZoom(ypos,col,pn)

				local xa, ya, rz, za = arrowEffects(ypos, col, pn, true)
				local m=activeMods[pn]
				local fakew = m.movew*m['movew'..col]*m.amovew*m['amovew'..col]
				local altx,alty = defaultx + xa,ypos + ya
				local offx,offy=getPropertyFromGroup('notes',v,"offsetX"),getPropertyFromGroup('notes',v,"offsetY")
			if m.incominganglex ~= 0 or m.incomingangley ~= 0 or m.incominganglez ~= 0 then
					local vec = {
						x = 0,
						y = ypos,
						z = 0
					}
					local rotatedpos = RotationXYZ(vec,m.incominganglex,m.incomingangley,m.incominganglez)
					altx = altx+rotatedpos.x
					alty = rotatedpos.y
					za = za + rotatedpos.z
			end
			if m['incominganglex'..col] ~= 0 or m['incomingangley'..col] ~= 0 or m['incominganglez'..col] ~= 0 then
				local vec = {
					x = 0,
					y = ypos,
					z = 0
				}
				local rotatedpos = RotationXYZ(vec,m['incominganglex'..col],m['incomingangley'..col],m['incominganglez'..col])
				altx= altx+rotatedpos.x
				alty= rotatedpos.y
				za = za + rotatedpos.z
			end
			if m.rotatex ~= 0 or m.rotatey ~= 0 or m.rotatez ~= 0 then
				local useX = altx
				local useY =alty
				local originPos = {x=screenWidth/2,y=screenHeight/2}
				if m.rotatetype == 1 then
					originPos.x = defaultx
				elseif m.rotatetype ~= 0 then
					originPos.x = m.rotateoriginx
					originPos.y = m.rotateoriginy
				end
				local vec = {
					x = useX - originPos.x,
					y = useY - originPos.y,
					z = za
				}
				local rotatedpos = RotationXYZ(vec,m.rotatex,m.rotatey,m.rotatez)
				rotatedpos.x = rotatedpos.x + originPos.x
				rotatedpos.y = rotatedpos.y + originPos.y
				altx= rotatedpos.x
				alty= rotatedpos.y
				za = rotatedpos.z
			end
			if m.rotationx ~= 0 or m.rotationy ~= 0 or m.rotationz ~= 0 then
				local useX = altx
				local useY = alty
				local originPos = {x=screenWidth/2-ARROW_SIZE/2,y=screenHeight/2}
				local xs = {92,204,316,428}
				local cx = (xs[4] + ARROW_SIZE - xs[1]) / 2 + xs[1]
					originPos.x = originPos.x + (pn*2-3)*cx
				local vec = {
					x = useX - originPos.x,
					y = useY - originPos.y,
					z = za
				}
				local rotatedpos = RotationXYZ(vec,m.rotationx,m.rotationy,m.rotationz)
				rotatedpos.x = rotatedpos.x + originPos.x
				rotatedpos.y = rotatedpos.y + originPos.y
				altx= rotatedpos.x
				alty= rotatedpos.y
				za = rotatedpos.z
			end

			local paths = updatePath(targTime-getSongPosition(),{x=altx,y=alty,z=za},col,m.infinite,getInfinite)
			altx = paths.x
			alty = paths.y

				local xawr,yawr,rzwr,zawr = arrowEffects(ypos,col,pn,false)
				local alp = arrowAlpha(ypos, col, pn, yawr)
			    local scalex, scaley = getScale(ypos, col, pn, defaultscale[c+1].x, defaultscale[c+1].y, true)
				if isSus then
					local ypos2 = getYAdjust(((getSongPosition()+.1) - targTime)*(downscroll and 1 or -1),col,pn) * scrollSpeeds * 0.45 - off + defaulty

					local xa2, ya2 = arrowEffects(ypos2, col, pn, true)
					--if scrollSpeed >= 0 then
				setPropertyFromGroup("notes",v,"angle",math.deg(math.atan2(((ypos2 + ya2)-(ypos + ya))*100,(xa2-xa)*100)+math.pi/2))
					--else
					--	note.angle = 180+math.deg(math.atan2(((ypos2 + ya2)-(ypos + ya))*100,(xa2-xa)*100) + math.pi/2)
					--end
				else
					setPropertyFromGroup("notes",v,"angle",rz)
				end

			local m = activeMods[pn]
              setPropertyFromGroup("notes",v,"x",altx)
            	setPropertyFromGroup("notes",v,"y",alty)
				-- from troll engine
				local alpha = getPropertyFromGroup('notes',v,'alpha')*math.clamp(scale(alp, 0.5, 0, 1, 0),0,1)
				local glow = math.clamp(scale(alp, 1, 0.5, 0, 1.3),0,1)

				setPropertyFromGroup('notes',v,'colorTransform.redMultiplier',1 - glow)
				setPropertyFromGroup('notes',v,'colorTransform.greenMultiplier',1 - glow)
				setPropertyFromGroup('notes',v,'colorTransform.blueMultiplier',1 - glow)
				setPropertyFromGroup('notes',v,'colorTransform.redOffset',glow*255)
				setPropertyFromGroup('notes',v,'colorTransform.greenOffset',glow*255)
				setPropertyFromGroup('notes',v,'colorTransform.blueOffset',glow*255)
				setPropertyFromGroup('notes',v,'colorTransform.alphaMultiplier',alpha)


	local fov = 90
	local tanHalfFOV = math.tan(math.rad(fov/2))
			local origin={x=getPropertyFromGroup('notes',v,"x") - (screenWidth/2),y= getPropertyFromGroup('notes',v,"y") - (screenHeight/2),z=za}
			local pos={x=origin.x,y=origin.y,z=origin.z/1000-fakew}

			local X = pos.x/(1/tanHalfFOV)/-pos.z+(screenWidth/2)
			local Y = pos.y/(1/tanHalfFOV)/-pos.z+(screenHeight/2)

			local scalenewy=isSus and 1 or scaley
			setPropertyFromGroup('notes', v, 'scale.x', scalex / -pos.z)
			local yscale = 1
			local susend = string.find(string.lower(getPropertyFromGroup('notes', v, 'animation.curAnim.name')), 'end') or string.find(string.lower(getPropertyFromGroup('notes', v, 'animation.curAnim.name')), 'tail')
			if isSus then
				if not susend then
					yscale = stepCrochet / 100 * 1.05 * scrollSpeeds
				else
					yscale = 1
				end
		    end
			setPropertyFromGroup('notes', v, 'scale.y', scalenewy / -pos.z * yscale)
			setPropertyFromGroup("notes",v,"scale.x",getPropertyFromGroup('notes',v,"scale.x")*zoom)
  			setPropertyFromGroup("notes",v,"scale.y",getPropertyFromGroup('notes',v,"scale.y")*zoom*(isSus and (1+m.longboy) or 1))
			setPropertyFromGroup('notes', v, 'x', X+offx)
			setPropertyFromGroup('notes', v, 'y', Y+offy)
			end

		end
		for i=0,getProperty('grpNoteSplashes.length')-1 do
			local pn = 2
			if activeMods[pn].hidenoteflash ~= 0 then
				setPropertyFromGroup('grpNoteSplashes',i,'visible',false)
			else
				setPropertyFromGroup('grpNoteSplashes',i,'visible',true)
			end
		end
	end
	updateCommand(elapsed,beat)
end

function onCreatePost()
	runHaxeCode([[
    var camNotes:FlxCamera = new FlxCamera();
    camNotes.height = game.camHUD.height;
    camNotes.bgColor = 0x00000000;

    FlxG.cameras.remove(game.camHUD, false);
    FlxG.cameras.remove(game.camOther, false);

    FlxG.cameras.add(camNotes, false);
    FlxG.cameras.add(game.camHUD, false);
    FlxG.cameras.add(game.camOther, false);

    game.notes.cameras = [camNotes];
    game.strumLineNotes.cameras = [camNotes];
    camNotes.angle = game.camHUD.angle;
    camNotes.zoom = game.camHUD.zoom;
    camNotes.x = game.camHUD.x;
    camNotes.y = game.camHUD.y;
	setVar("camNotes",camNotes);
]])
	TEMPLATE.InitMods()

	--WRITE MODS HERE!

	initCommand()

	--must be called at END of start
	TEMPLATE.setup()

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
	-- plugins
	-- there are many functions that can be used in fnf
	local function flicker()
		if (getSongPosition() * 0.001 * 60) % 2 < 1 then
			return -1
		else
			return 1
		end
	end
    local function hide(t)
		local bt,tpn = t[1],t.pn
		for i=0,3 do
			ease{bt+i*.125-1,.5,outExpo,-70,'movey'..i,pn=tpn}
			ease{bt+i*.125-.5,1.25,inExpo,650,'movey'..i,pn=tpn}
			set{bt+i*.125+1.75,1,'stealth',pn=tpn}
			set{bt+i*.125+1.75,1,'dark',pn=tpn}
		end
	end
	local function unhide(t)
		local bt,tpn = t[1],t.pn
		for i=0,3 do
			set{bt+i*.125-2,0,'stealth',pn=tpn}
			set{bt+i*.125-2,0,'dark',pn=tpn}
			ease{bt+i*.125-2,1,outExpo,-70,'movey'..i,pn=tpn}
			ease{bt+i*.125-1,1,inExpo,50,'movey'..i,pn=tpn}
			ease{bt+i*.125-0,1.25,outElastic,0,'movey'..i,pn=tpn}
		end
	end
	--wiggle(beat,num,div,ease,amt,mod)
	local function wig(t)
		local b,num,div,ea,am,mo = t[1],t[2],t[3],t[4],t[5],t[6]
		local f = 1
		for i=0,num do
			local smul = i==0 and 1 or 0
			local emul = i==num and 0 or 1

			ease{b+i*(1/div),1/div,ea,startVal = am*smul*f, am*emul*-f,mo,pn=t.pn}

			f = f*-1
		end
	end
	--simple mod 2
	local function sm2(tab)
		local b,len,eas,amt,mods,intime = tab[1],tab[2],tab[3],tab[4],tab[5],tab.intime
		if not intime then intime = .1 end
		if intime <= 0 then intime = .001 end
		ease{b-intime,intime,linear,amt,mods,pn=tab.pn}
		ease{b,len-intime,eas,0,mods,pn=tab.pn}
	end
	local function mod_bounce2(t)
		local pn = t.pn
		local beat,length,ease1,ease2,amt,mod,pn = t[1],t[2],t[3],t[4],t[5],t[6],t.plr
		ease{beat, (length/2), ease1, amt, mod, plr=pn}
		ease{beat+(length/2), (length/2), ease2, -amt, mod, plr=pn}
	end
	local function mod_ease(beat, len, str1, str2, mod, t, eas, pn)
		if t == 'end' then len = len - beat end
		set {beat, str1, mod, plr = pn}
		ease {beat, len, eas, str2, mod, plr = pn}
	end
	local function mod_outin(beat, len, per1, per2, mod, oute, ine)
		ease{beat, len/2, oute or outCirc, per2, mod}
		ease{beat+len/2, len/2, ine or inCirc, per1, mod}
	end
	local function mod_kick(beat,length,start,apex,mod,inEase,outEase,pn)
		local off = length/2
		mod_ease(beat - off, (length/2), start, apex, tostring(mod), 'len', inEase,pn)
		mod_ease((beat+(length/2)) - off, (length/2), apex, start, tostring(mod), 'len', outEase,pn)
	end
	local function mod_message(tab)
		func{tab[1], function() debugPrint(tab[2]) end}
	end
	local function mod_perframe(start_beat, end_beat, f)
		func {start_beat, end_beat, f, mode = 'end'}
	end
	local function mod_insert(beat, len, modstring, t, pn)
		local string_gmatch = string.gfind or string.gmatch
		if t == 'end' then len = len - beat end
		for str in string_gmatch(modstring, '[^,]+') do
			local str = string.gsub(str, '%%', '')
			local activation_rate = 1
			local percent = 100
			local modname = nil
			for part in string_gmatch(str, '[^ ]+') do
				if string.find(part, '*') then
					activation_rate = tonumber(string.sub(part, 2)) or activation_rate
				elseif not string.find(part, '[^%d]') then
					percent = tonumber(part) or (part == 'no' and 0) or percent
				elseif string.find(part, '^%d+x$') then
					local _, x = string.match(part, '^(%d+)x$')
					modname = 'xmod'
					percent = tonumber(x)
				elseif string.find(part, '^c%d+$') then
					local _, c = string.match(part, '^c(%d+)$')
					modname = 'cmod'
					percent = tonumber(c)
				else
					modname = part
				end
			end
			if modname then
				if activation_rate < 0 or activation_rate > 9998 then
					set {beat, percent/100, modname, plr = pn}
				elseif activation_rate == 0 then
					-- do nothing
				else -- activation_rate > 0
					msg('modstring \'' .. str .. '\' needs a *-1 activation rate to work with the Mirin Template backend.')
				end
			end
		end
	end
	local function mixEase(e1, e2, point)
		if not point then point = 0.5 end

		return function(a)
			if a < point then
				return e1(a / point) * point
			else
				return e2((a - point) / (1 - point)) * (1 - point) + point
			end
		end
	end
	local function ease_smooth(beat,len,amount,mod,inEase,outEase,point,pn)
		ease {beat-(len*point),len, mixEase(inEase,outEase,point), amount, mod, plr = pn}
	end
	local function sugarkiller(beat, length, step, minstealth, maxstealth, plr)
		if not minstealth then minstealth = 0.50 end;
		if not maxstealth then maxstealth = 0.85 end;
		if not step then step = 1 end;
		if not length then length = 1 end;
		for i = 0, math.max(length-1,0) do
			ease{beat+(i*0.5), .25/step, instant, 1.00, 'invert', 0, 'flip', maxstealth, 'stealth', plr = plr}
			ease{beat+(i*0.5)+.25/step, .25/step, instant, 1.00, 'flip', 0, 'invert', maxstealth, 'stealth', plr = plr}
			ease{beat+(i*0.5)+.50/step, .25/step, instant, 1.00, 'flip', -1.00, 'invert', maxstealth, 'stealth', plr = plr}
			if i == math.max(length-1,0) then
				ease{beat+(i*0.5)+.75/step, .25/step, instant, 0, 'flip', 0, 'invert', 0, 'stealth', plr = plr}
			else
				ease{beat+(i*0.5)+.75/step, .25/step, instant, 0, 'flip', 0, 'invert', minstealth, 'stealth', plr = plr}
			end
		end
	end
	local function ease2 (beat,length,eas,amount,mod,pn)
		ease{beat-(length/2),length,eas,amount,mod,plr = pn}
	end
	local function mod_beat(beat,strength,pn)
		if not strength then strength = 1000 end;
			set {beat-.5, strength, 'beat', plr = pn}
			set {beat+.5, 0, 'beat', plr = pn}
	end
	local function range(var,tablex)
		for _,v in pairs(tablex) do
			if var >= v[1] and var <= v[2] then
				return true
			end
		end
		return false
	end
	local function modulo(a, b)
		return a - math.floor(a / b) * b
	end
	local stringbuilder_mt =  {
		__index = {
			build = table.concat,
			clear = iclear,
		},
		__call = function(self, a)
			table.insert(self, tostring(a))
			return self
		end,
		__tostring = table.concat,
		}
	local function stringbuilder()
			return setmetatable({}, stringbuilder_mt)
	end
	--end mod plugins

	local me = ease
	local e = mod_ease
	local mpf = func
	local m2 = func
	local msg = mod_message
	local mi = mod_insert

end

bTestMode = true -- debug

function onCreate()
	if bTestMode then setProperty('skipCountdown',true) end
end

function updateCommand(elapsed,beat)
end
--end callbacks

return 0
