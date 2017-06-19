"""
    FEMAP695

A module containing functions for use with FEMA P695 studies.

# Functions
* ACMRxx - Compute the acceptable collapse margin ratio.
* beta_total - Compute the total uncertainty present in the system.
* mappedValue - Retrieve mapped seismic demand parameters.
* SF1 - Calculate scale factor 1.
* SMT - Calculate intensity of the maximum considered earthquake.
* SSF - Calculate the spectral shape factor.
"""
module FEMAP695

using Distributions, Roots, Interpolations


"""
    ACMRxx(beta_total, collapseProbability, X_initial=0.622)

Compute the acceptable collapse margin ratio for a given system uncertainty and collapse
probability.

Due to the lack of an inverse cumulative distribution function, the intermediate value X
(the inverse of the acceptable collapse margin ratio) must be solved for. The initial guess
X_initial has a default value equal to the inverse of the average ACMR20 value in Table 7-3.
It may be specified if there are convergence issues.
"""
function ACMRxx(beta_total, collapseProbability, X_initial=0.622)
    dist = Distributions.LogNormal(0, beta_total)
    f(x) = Distributions.cdf(dist,x) - collapseProbability
    X    = Roots.newton(f,X_initial)     # Don't have an inverse cdf, so need to solve
    return 1/X

end  # function ACMRxx


"""
    beta_total(rating_DR, rating_TD, rating_MDL, mu_T)

Compute the total uncertainty present in the system.

# Arguments
* `rating_DR`: qualitative rating of the design requirements, ranging from "a" to "d"
* `rating_TD`: qualitative rating of the test data, ranging from "a" to "d"
* `rating_MDL`: qualitative rating of the model, ranging from "a" to "d"
* `mu_T`: period-based ductility obtained from pushover tests

Refer to section 7.3.1 for details.
"""
function beta_total(rating_DR::String, rating_TD::String, rating_MDL::String, mu_T::Float64)
    if     ismatch(r"A"i, rating_DR)
        beta_DR = 0.10
    elseif ismatch(r"B"i, rating_DR)
        beta_DR = 0.20
    elseif ismatch(r"C"i, rating_DR)
        beta_DR = 0.35
    elseif ismatch(r"D"i, rating_DR)
        beta_DR = 0.50
    else
        error("Unknown rating_DR: $(rating_DR)")
    end

    if     ismatch(r"A"i, rating_TD)
        beta_TD = 0.10
    elseif ismatch(r"B"i, rating_TD)
        beta_TD = 0.20
    elseif ismatch(r"C"i, rating_TD)
        beta_TD = 0.35
    elseif ismatch(r"D"i, rating_TD)
        beta_TD = 0.50
    else
        error("Unknown rating_TD: $(rating_TD)")
    end

    if     ismatch(r"A"i, rating_MDL)
        beta_MDL = 0.10
    elseif ismatch(r"B"i, rating_MDL)
        beta_MDL = 0.20
    elseif ismatch(r"C"i, rating_MDL)
        beta_MDL = 0.35
    elseif ismatch(r"d"i, rating_MDL)
        beta_MDL = 0.50
    else
        error("Unknown rating_MDL: $(rating_MDL)")
    end

    beta_RTR = min(0.1 + 0.1*mu_T, 0.4)

    beta = sqrt(beta_RTR^2 + beta_DR^2 + beta_TD^2 + beta_MDL^2)
    return round(beta*40)/40    # Round to nearest 0.025

end  # function beta_total

"""
    mappedValue(value, sdc)

Retrieve the mapped seismic demand parameter for the given seismic design category.

# Arguments
* `value`: seismic demand parameter.
  - Permitted values: `SS`, `S1`, `Fa`, `Fv`, `SMS`, `SM1`, `SDS`, `SD1`, `TS`
* `sdc`: seismic design category.
  - Permitted values: `Dmax`, `Dmin`, `Cmax`, `Cmin`, `Bmax`, `Bmin`

"""
function mappedValue(value, sdc)
    if     ismatch(r"Dmax"i, sdc)
        SS  = 1.5;      Fa  = 1.0;      SMS = 1.5;      SDS = 1.0
        S1  = 0.60;     Fv  = 1.50;     SM1 = 0.90;     SD1 = 0.60
    elseif ismatch(r"(Cmax|Dmin)"i, sdc)
        SS  = 0.55;     Fa  = 1.36;     SMS = 0.75;     SDS = 0.50
        S1  = 0.132;    Fv  = 2.28;     SM1 = 0.30;     SD1 = 0.20
    elseif ismatch(r"(Bmax|Cmin)"i, sdc)
        SS  = 0.33;     Fa  = 1.53;     SMS = 0.50;     SDS = 0.33
        S1  = 0.083;    Fv  = 2.4;      SM1 = 0.20;     SD1 = 0.133
    elseif ismatch(r"Bmin"i, sdc)
        SS  = 0.156;    Fa  = 1.6;      SMS = 0.25;     SDS = 0.167
        S1  = 0.042;    Fv  = 2.4;      SM1 = 0.10;     SD1 = 0.067
    else
        error("Unknown seismic design category: $(sdc)")
    end

    if     ismatch(r"SS"i,  value); return SS
    elseif ismatch(r"S1"i,  value); return S1
    elseif ismatch(r"Fa"i,  value); return Fa
    elseif ismatch(r"Fv"i,  value); return Fv
    elseif ismatch(r"SMS"i, value); return SMS
    elseif ismatch(r"SM1"i, value); return SM1
    elseif ismatch(r"SDS"i, value); return SDS
    elseif ismatch(r"SD1"i, value); return SD1
    elseif ismatch(r"TS"i,  value); return SD1/SDS
    else
        error("Unknown value: $(value)")
    end

end  # function mappedValue


"""
    SF1(T, sdc, gmset="farfield")

Calculate scale factor 1, used to scale the intensity of the ground motions to the intensity
of the maximum considered earthquake.
"""
function SF1(T, sdc, gmset="farfield")
    if ismatch(r"farfield"i, gmset)
        T_interp    = [0.25, 0.30, 0.35, 0.40, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2,
                       1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0]
        SNRT_interp = [0.785, 0.781, 0.767, 0.754, 0.755, 0.742, 0.607, 0.541, 0.453, 0.402,
                       0.350, 0.303, 0.258, 0.210, 0.169, 0.149, 0.134, 0.119, 0.106, 0.092,
                       0.081, 0.063, 0.053, 0.046, 0.041]
    elseif ismatch(r"nearfield"i, gmset)
        error("Not implemented: $(gmset)")
    else
        error("Unknown ground motion set: $(gmset)")
    end

    if T <= T_interp[1] || T >= T_interp[end]
        error("Period is out of range: T = $(T)")
    end

    itp  = interpolate((T_interp,),SNRT_interp,Gridded(Linear()))
    SNRT = itp[T]
    SMT  = FEMAP695.SMT(T,sdc)
    return SMT/SNRT
end # function SF1

"""
    SMT(T, sdc)

Calculate the intensity of the maximum considered earthquake at a given period for a given
seismic design category.
"""
function SMT(T,sdc)
    SM1 = mappedValue("SM1",sdc)
    SMS = mappedValue("SMS",sdc)
    if T <= SM1/SMS
        return SMS
    else
        return SM1/T
    end

end # function SMT

"""
    SSF(T, mu_t, sdc)

Compute the spectral shape factor.
"""
function SSF(T, mu_t, sdc)
    if ismatch(r"Dmax"i, sdc)
        Z_SSF = [
            1.00 1.05 1.10 1.13 1.18 1.22 1.28 1.33
            1.00 1.05 1.11 1.14 1.20 1.24 1.30 1.36
            1.00 1.06 1.11 1.15 1.21 1.25 1.32 1.38
            1.00 1.06 1.12 1.16 1.22 1.27 1.35 1.41
            1.00 1.06 1.13 1.17 1.24 1.29 1.37 1.44
            1.00 1.07 1.13 1.18 1.25 1.31 1.39 1.46
            1.00 1.07 1.14 1.19 1.27 1.32 1.41 1.49
            1.00 1.07 1.15 1.20 1.28 1.34 1.44 1.52
            1.00 1.08 1.16 1.21 1.29 1.36 1.46 1.55
            1.00 1.08 1.16 1.22 1.31 1.38 1.49 1.58
            1.00 1.08 1.17 1.23 1.32 1.40 1.51 1.61 ]
    elseif ismatch(r"(Dmin|Cmax|Cmin|Bmax|Bmin)"i, sdc)
        Z_SSF = [
            1.00 1.02 1.04 1.06 1.08 1.09 1.12 1.14
            1.00 1.02 1.05 1.07 1.09 1.11 1.13 1.16
            1.00 1.03 1.06 1.08 1.10 1.12 1.15 1.18
            1.00 1.03 1.06 1.08 1.11 1.14 1.17 1.20
            1.00 1.03 1.07 1.09 1.13 1.15 1.19 1.22
            1.00 1.04 1.08 1.10 1.14 1.17 1.21 1.25
            1.00 1.04 1.08 1.11 1.15 1.18 1.23 1.27
            1.00 1.04 1.09 1.12 1.17 1.20 1.25 1.30
            1.00 1.05 1.10 1.13 1.18 1.22 1.27 1.32
            1.00 1.05 1.10 1.14 1.19 1.23 1.30 1.35
            1.00 1.05 1.11 1.15 1.21 1.25 1.32 1.37 ]
    else
        error("Unknown seismic design category: $(sdc)")
    end

    @assert mu_t>=1 "mu_t must be greater than or equal to 1"

    X_mu_t = [1.0, 1.1, 1.5, 2, 3, 4, 6, 8]
    Y_T    = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]

    if T <= 0.5
        if mu_t >= 8
            ssf = Z_SSF[1,end]
        else
            itp = interpolate((X_mu_t,), Z_SSF[1,:], Gridded(Linear()))
            ssf = itp[mu_t]
        end
    elseif T >= 1.5
        if mu_t >= 8
            ssf = Z_SSF[end,end]
        else
            itp = interpolate((X_mu_t,), Z_SSF[end,:], Gridded(Linear()))
            ssf = itp[mu_t]
        end
    else
        if mu_t >= 8
            itp = interpolate((Y_T,), Z_SSF[:,end], Gridded(Linear()))
            ssf = itp[T]
        else
            itp = interpolate((Y_T,X_mu_t), Z_SSF, Gridded(Linear()))
            ssf = itp[T,mu_t]
        end
    end

    return ssf
end # function SSF

end # module FEMAP695
