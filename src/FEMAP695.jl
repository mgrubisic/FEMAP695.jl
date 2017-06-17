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
    rating_DR  = lowercase(rating_DR)
    rating_TD  = lowercase(rating_TD)
    rating_MDL = lowercase(rating_MDL)

    if     rating_DR == "a"
        beta_DR = 0.10
    elseif rating_DR == "b"
        beta_DR = 0.20
    elseif rating_DR == "c"
        beta_DR = 0.35
    elseif rating_DR == "d"
        beta_DR = 0.50
    else
        error("Unknown rating_DR: $(rating_DR)")
    end

    if     rating_TD == "a"
        beta_TD = 0.10
    elseif rating_TD == "b"
        beta_TD = 0.20
    elseif rating_TD == "c"
        beta_TD = 0.35
    elseif rating_TD == "d"
        beta_TD = 0.50
    else
        error("Unknown rating_TD: $(rating_TD)")
    end

    if     rating_MDL == "a"
        beta_MDL = 0.10
    elseif rating_MDL == "b"
        beta_MDL = 0.20
    elseif rating_MDL == "c"
        beta_MDL = 0.35
    elseif rating_MDL == "d"
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

Retrieve the mapped seismic demand parameter for the given seismic design category


"""
function mappedValue(value, sdc)
    sdc   = lowercase(sdc)
    value = lowercase(value)

    if     sdc == "dmax"
        SS  = 1.5;      Fa  = 1.0;      SMS = 1.5;      SDS = 1.0
        S1  = 0.60;     Fv  = 1.50;     SM1 = 0.90;     SD1 = 0.60
    elseif sdc == "cmax" || sdc == "dmin"
        SS  = 0.55;     Fa  = 1.36;     SMS = 0.75;     SDS = 0.50;
        S1  = 0.132;    Fv  = 2.28;     SM1 = 0.30;     SD1 = 0.20;
    elseif sdc == "bmax" || sdc == "cmin"
        SS  = 0.33;     Fa  = 1.53;     SMS = 0.50;     SDS = 0.33;
        S1  = 0.083;    Fv  = 2.4;      SM1 = 0.20;     SD1 = 0.133;
    elseif sdc == "bmin"
        SS  = 0.156;    Fa  = 1.6;      SMS = 0.25;     SDS = 0.167;
        S1  = 0.042;    Fv  = 2.4;      SM1 = 0.10;     SD1 = 0.067;
    else
        error("Unknown seismic design category: $(sdc)")
    end

    if value == "ss"
        return SS
    elseif value == "s1"
        return S1
    elseif value == "fa"
        return Fa
    elseif value == "fv"
        return Fv
    elseif value == "sms"
        return SMS
    elseif value == "sm1"
        return SM1
    elseif value == "sds"
        return SDS
    elseif value == "sd1"
        return SD1
    elseif value == "ts"
        return SD1/SDS
    else
        error("Unknown value: $(value)")
    end

end  # function mappedValue


"""
    SF1(T, sdc, set="farfield")

Calculate scale factor 1, used to scale the intensity of the ground motions to the intensity
of the maximum considered earthquake.
"""
function SF1(T, sdc, set="farfield")
    if set == "farfield"
        T_interp    = [0.25, 0.30, 0.35, 0.40, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2,
                       1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0]
        SNRT_interp = [0.785, 0.781, 0.767, 0.754, 0.755, 0.742, 0.607, 0.541, 0.453, 0.402, 0.350, 0.303,
                       0.258, 0.210, 0.169, 0.149, 0.134, 0.119, 0.106, 0.092, 0.081, 0.063, 0.053, 0.046, 0.041]
    elseif set == "nearfield"
        error("Not implemented: $(set)")
    else
        error("Unknown ground motion set: $(set)")
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

end # module FEMAP695
