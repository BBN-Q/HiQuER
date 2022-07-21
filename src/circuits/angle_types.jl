#=
Copyright 2022 Raytheon BBN

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.

You may obtain a copy of the License at
   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
=#

using LaTeXStrings

abstract type AbstractAngle end

struct Angle <: AbstractAngle
    value::Rational
end

struct FloatAngle <: AbstractAngle
    value::Float64
end

FloatAngle(x::Angle) = FloatAngle(Float64(x.value))

function _to_latex_raw(ang::Angle)
    n = numerator(ang.value)
    d = denominator(ang.value)
    if isone(abs(n))
        numstr = sign(n)*sign(d) < 0 ? "-" : ""
    else
        numstr = string(n)
    end
    if isone(abs(d))
        out = string(numstr, "\\pi")
    else
        out = string("\\frac{", numstr, "\\pi}{", d, "}")
    end
    return out
end

_to_latex(ang::Angle) = latexstring(_to_latex_raw(ang))

Base.show(io::IO, ::MIME"application/x-latex", angle::Angle) = print(io, _to_latex(angle))
Base.show(io::IO, ::MIME"text/latex", angle::Angle) = print(io, _to_latex(angle))
Base.show(io::IO, angle::Angle) = print(io, string(numerator(angle.value), "*pi/", denominator(angle.value)))
    
Base.show(io::IO, ang::FloatAngle) = print(io, string(ang.value))
Base.:(==)(x::T, y::T) where T<:AbstractAngle   = x.value == y.value
Base.:(<)(x::T, y::T)  where T<:AbstractAngle   = x.value < y.value
Base.:(+)(x::T, y::T)  where T<:AbstractAngle   = T(x.value + y.value)
Base.:(-)(x::T, y::T)  where T<:AbstractAngle   = T(x.value - y.value)

Base.:(*)(x::Float64, y::T) where T<:AbstractAngle = T(x*y.value)
Base.:(*)(y::T, x::Float64) where T<:AbstractAngle = T(x*y.value)

export Angle, FloatAngle
