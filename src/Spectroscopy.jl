#   My own spectroscopy functions
#   The definition of linehapesGauss and lineshapesLorentz  is inspired from Charles Le Losq
#

module Spectroscopy

using Unitful, FC

## myChop, myChopImaginary

#The `myChop` module provides the functions `myChop` and `myChopImaginary`, which replace
#numbers that are close to zero with zero. These functions act recursivley on
#collections.

#    zeps::Float64 The default lower threshold for a number to be replaced by zero.

const myEps = 1e-14

#myChop(x::T, eps::Real = myEps)
#Replace `x` by zero if `abs(x) < eps`. `myChop` acts recursively on
#mappable objects. Objects that cannot be sensibly compared to a real
#number are passed unaltered.

function myChop(x::Real, eps::Real = myEps)
	return abs(x) > eps ? x : 0
end

function myChop(x::Complex, eps::Real = myEps)
	return complex(myChop(real(x), eps), myChop(imag(x), eps))
end

function myChopImaginary(x::Real, eps::Real = myEps)
	return myChop(x)
end

function myChopImaginary(x::Complex, eps::Real = myEps)
	return abs(myChop(imag(x)))> eps ? myChop(x) : myChop(real(x))
end

function myChop(a::AbstractArray, eps::Real = myEps)
	for i in firstindex(a):lastindex(a)
		a[i] = myChop(a[i])
	end
	return a
end

function myChopImaginary(a::AbstractArray, eps::Real = myEps)
	for i in firstindex(a):lastindex(a)
		a[i] = myChopImaginary(a[i])
	end
	return a
end


####################################################
####################################################
#These functions allows converting two measurements. TWO METHODS
# Quantity is the quantity to convert
# myunit is the target unit
# delta is the uncertaity on the quantity to be converted to my unit.
# units defined are wavelength, eV, omega, wavenumber.
# ONLY ANGULAR FREQUENCIES ARE CONSIDERED


function convertme(quantity::Unitful.Quantity, myunit::Unitful.FreeUnits)
# This is for wavelength
if dimension(quantity)==dimension(u"m")
	if dimension(myunit)==dimension(u"m")  return uconvert(myunit,1.0*quantity) end
	if dimension(myunit)==dimension(u"eV") return uconvert(myunit, 1.0u"eV"*ustrip(FC.h*FC.c/(FC.e*uconvert(u"m",quantity)))) end
	if dimension(myunit)==dimension(u"Hz") return uconvert(myunit, 1.0u"Hz"*2*π *ustrip(FC.c/uconvert(u"m",quantity))) end
	if dimension(myunit)==dimension(u"m^-1") return  1.0*uconvert(myunit,1/(uconvert(u"m",quantity))) end
end
# This is for energy in eV
if dimension(quantity)==dimension(u"eV")
	if dimension(myunit)==dimension(u"m") return uconvert(myunit,1.0u"m"*ustrip(FC.h*FC.c/(FC.e*uconvert(u"eV",quantity)))) end
	if dimension(myunit)==dimension(u"eV") return 1.0*uconvert(myunit,quantity) end
	if dimension(myunit)==dimension(u"Hz") return uconvert(myunit, 1.0u"Hz"* ustrip(FC.e*uconvert(u"eV",quantity)/ FC.hbar)) end
	if dimension(myunit)==dimension(u"m^-1") return  uconvert(myunit,u"1/m"*ustrip(FC.e*uconvert(u"eV",quantity)/ (FC.h*FC.c))) end
end
# This is for ANGULAR frequency in Hz
if dimension(quantity)==dimension(u"Hz")
	if dimension(myunit)==dimension(u"m") return uconvert(myunit, u"m"*2*π*ustrip(FC.c/uconvert(u"Hz",quantity))) end
	if dimension(myunit)==dimension(u"eV")  return  uconvert(myunit,1.0u"eV"*ustrip(FC.hbar*uconvert(u"Hz",quantity)/FC.e)) end
	if dimension(myunit)==dimension(u"Hz") return uconvert(myunit,quantity) end
	if dimension(myunit)==dimension(u"m^-1") return uconvert(myunit, 1.0u"m^-1"*ustrip(uconvert(u"Hz",quantity)/(2*π*FC.c)))  end
end
# This is for wavenumber
if dimension(quantity)==dimension(u"m^-1")
	if dimension(myunit)==dimension(u"m")  return uconvert(myunit,1.0/quantity) end
	if dimension(myunit)==dimension(u"eV") return uconvert(myunit, 1.0u"eV"*ustrip(FC.h*FC.c/(FC.e*uconvert(u"m",1/quantity)))) end
	if dimension(myunit)==dimension(u"Hz") return uconvert(myunit, 1.0u"Hz"*2*π *ustrip(FC.c) /(ustrip(uconvert(u"m",1/quantity)))) end
	if dimension(myunit)==dimension(u"m^-1") return  1.0*uconvert(myunit,1/(uconvert(u"m",1/quantity))) end
end
end

function convertme(quantity::Unitful.Quantity, delta::Unitful.Quantity, myunit::Unitful.FreeUnits)
	quantity1=convertme(quantity,unit(delta))
	quantity2=convertme(quantity,myunit)
	return quantity2*delta/quantity1
end

####################################################
####################################################

# lineshapeGauss( intensity::Array{Float64}, position::Array{Float64}, fwhm::Array{Float64}, x::Array{Float64})
#
# You can enter the intensity, position and full-width at half-maximum (fwhm) values as arrays of float 64
# (even containing one float value), without specifying style.
#
#INPUTS:
#	intIntensity: Array{Float64} containing the INTEGRATED intensities
#	position: Array{Float64} containing the peaks positions;
#   fwhm: Array{Float64} containing the peaks full-width at half maximum;
#	x: Array{Float64} containing the x axis values;
#
#  OUTPUTS:
#	int_sum: Array{Float64} is the summed intensities;
#	int_individual: Array{Float64} is an array of the individual intensities.

#----------
# Examples
#----------
# To have four gaussian peaks centered at 800, 900, 1000 and 1100 cm-1 with fwhm of 50 cm-1 on a Raman spectrum, you will enter:
# int_sum, int_individual = lineshapeGauss([1.0,1.0,1.0,1.0], [800.0,900.0,1000.0,1100.0], [50.0,50.0,50.0,50.0], x)
# and int_individual will contain in 4 columns the 4 different y values of the peaks, and int_sum their sum.

function lineshapeGauss(intIntensity::Array{Float64}, position::Array{Float64}, fwhm::Array{Float64}, x::Array{Float64})
    segments = zeros(size(x)[1],size(intIntensity)[1]) #Create a 2D array of size z and the number of gaussians requested
    for i = 1:size(intIntensity)[1]
            segments[:,i] = intIntensity[i].*2 .*sqrt.(log.(2)./pi)./(fwhm[i]) .* exp.(-4 .*log.(2) .* ((x[:,1].-position[i])./(fwhm[i])).^2) #  WOW
        end
    return vec(sum(segments,dims=2)), segments #sum of all segments, but return individual segments as well
end


# I am adding this one, without arrays
function lineshapeGauss(intIntensity::Float64, position::Float64, fwhm::Float64, x::Array{Float64})
    segment = zeros(length(x)) #Create a 2D array of size x
	segment = intIntensity.* exp.(-4 .*log.(2) .* ((x.-position)./(fwhm)).^2)
    return segment
end


####################################################
####################################################
# function lineshapeLorentz
# INPUTS:
#	intIntensity: Array{Float64} containing the peaks INTEGRATED intensities;
#	position: Array{Float64} containing the peak positions;
#	fwhm: Array{Float64} containing the peaks falf-width at half maximum;
#	x: Array{Float64} containing the x axis values;

#OUTPUTS:
#	int_sum: Array{Float64} is the summed intensities;
#	int_individual: Array{Float64} is an array of the individual intensities.

function lineshapeLorentz(intIntensity::Array{Float64},position::Array{Float64},fwhm::Array{Float64}, x::Array{Float64})
	    segments = zeros(size(x)[1],size(intIntensity)[1])
	     for i = 1:size(intIntensity)[1]
	            segments[:,i] = intIntensity[i]./pi .* (fwhm[i]./2) ./ ( (x[:,1].-position[i]).^2 .+ (fwhm[i]./2).^2 )
	        end
	    return vec(sum(segments,dims=2)), segments
	end


####


####
function lineshapeLorentz(intIntensity::Float64,position::Float64,fwhm::Float64, x::Array{Float64})
	segment = zeros(length(x))
	segment = intIntensity./pi .* (fwhm./2) ./ ( (x.-position).^2 .+ (fwhm./2).^2 )
	return segment 
end


#### Black Body radiation
function
	blackBodyLambda(lambda::Unitful.Quantity,temperature::Unitful.Quantity)
	if dimension(lambda) != dimension(u"m") return -1 end
	if dimension(temperature) != dimension(u"K") return -1 end
	lambda=uconvert(u"m",lambda)
	temperature=uconvert(u"K",temperature)
	return 8*pi*FC.h*(FC.c)^2/lambda^5* (exp(FC.h*FC.c/(lambda *FC.k_B*temperature))-1)^(-1)
end

###

function
	phononSpectralDensity(lambda::Unitful.Quantity,temperature::Unitful.Quantity,speedOfSound::Unitful.Quantity)
	if dimension(lambda) != dimension(u"m") return -1 end
	if dimension(temperature) != dimension(u"K") return -1 end
	if dimension(speedOfSound) != dimension(u"m/s") return -1 end
	lambda=uconvert(u"m",lambda)
	temperature=uconvert(u"K",temperature)
	speedOfSound=uconvert(u"m/s",speedOfSound)
	return 8*pi*FC.h*speedOfSound^2/lambda^5* (exp(FC.h*speedOfSound/(lambda *FC.k_B*temperature))-1)^(-1)
end



function surface_to_vecs(x::Array{Float64,1}, y::Array{Float64,1}, s::Array{Float64,2})
   a = Array(s)
   xn = Vector{eltype(x)}(undef, length(a))
   yn = Vector{eltype(y)}(undef, length(a))
   zn = Vector{eltype(s)}(undef, length(a))
   for (n, (i, j)) in enumerate(Tuple.(CartesianIndices(a)))
       xn[n] = x[j]
       yn[n] = y[i]
       zn[n] = a[i, j]
   end
   return xn, yn, zn
end


export  lineshapeGauss, lineshapeLorentz, convertme, blackBodyLambda, phononSpectralDensity, myChop, myChopImaginary, surface_to_vecs


#####


#   My own numerical functions
#   Trapezoidal integration is inspired from Charles Le Losq

# trapezoidalInt({Tx<:Number, Ty<:Number}(x::Vector{Tx}, y::Vector{Ty})
# Trapezoidal integration.

# INPUTS:
# x: Vector{Float64} containing the x values;
# y: Vector{Float64} containing the y values. It works on complex number as well

# OUTPUTS:
# area: Vector{Float64}, the trapezoidal integration value.
# This function is particularly helpful to calculate the area under a portion of a spectrum, and can be used for various purposes (normalisation, area comparison, etc.).

function trapezoidalInt(x::Vector{Tx}, y::Vector{Ty}) where {Tx<:Number,Ty<:Number}
    # Trapezoidal integration rule
    local n = length(x)
    if (length(y) != n)
        error("Both vectors must be of same length for integration to work")
    end
    r = zero(zero(Tx) + zero(Ty)) # this creates a number whose type is the most general of Tx and Ty
    if n == 1
        return r
    end
    for i in 2:n
        r += (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    integration = r / 2
    return integration
end

# intervalfunction({Tx<:Number, Ty<:Number}(x::Vector{Tx}, y::Vector{Ty})
# Trapezoidal integration.

# INPUTS:
# x: Vector{Float64} containing the x values;
# y: Vector{Float64} containing the y values. It works on complex number as well

# OUTPUTS:
# area: Vector{Float64}, the trapezoidal integration value.
# This function is particularly helpful to calculate the area under a portion of a spectrum, and can be used for various purposes (normalisation, area comparison, etc.).

function heaviside(t::Float64)
    0.5 * (sign(t) + 1)
end

function interval(t::Float64, a::Float64, b::Float64)
    heaviside(t - a) - heaviside(t - b)
end



export lineshapeGauss, lineshapeLorentz, convertme, blackBodyLambda, phononSpectralDensity, myChop, myChopImaginary, surface_to_vecs, trapezoidalInt, heaviside, interval

end # module
