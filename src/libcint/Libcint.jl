module Libcint

using libcint_jll

#####
##### Functions names
#####

const dim = [
    "CINTgto_cart",
    "CINTgto_spheric",
    "CINTgto_spinor",
    ]

const f1e = [
    "cint1e_ovlp_cart",
    "cint1e_nuc_cart",
    "cint1e_kin_cart",
    "cint1e_ia01p_cart",
    "cint1e_irixp_cart",
    "cint1e_ircxp_cart",
    "cint1e_iking_cart",
    "cint1e_iovlpg_cart",
    "cint1e_inucg_cart",
    "cint1e_ipovlp_cart",
    "cint1e_ipkin_cart",
    "cint1e_ipnuc_cart",
    "cint1e_iprinv_cart",
    "cint1e_rinv_cart",
    "cint2e_cart",
    "cint2e_ig1_cart",
    "cint2e_ip1_cart",
    "cint1e_ovlp_sph",
    "cint1e_nuc_sph",
    "cint1e_kin_sph",
    "cint1e_ia01p_sph",
    "cint1e_irixp_sph",
    "cint1e_ircxp_sph",
    "cint1e_iking_sph",
    "cint1e_iovlpg_sph",
    "cint1e_inucg_sph",
    "cint1e_ipovlp_sph",
    "cint1e_ipkin_sph",
    "cint1e_ipnuc_sph",
    "cint1e_iprinv_sph",
    "cint1e_rinv_sph",
    "cint2e_sph",
    "cint2e_ig1_sph",
    "cint2e_ip1_sph",
    "cint1e_ovlp",
    "cint1e_nuc",
    "cint1e_nucg",
    "cint1e_srsr",
    "cint1e_sr",
    "cint1e_srsp",
    "cint1e_spsp",
    "cint1e_sp",
    "cint1e_spspsp",
    "cint1e_spnuc",
    "cint1e_spnucsp",
    "cint1e_srnucsr",
    "cint1e_sa10sa01",
    "cint1e_ovlpg",
    "cint1e_sa10sp",
    "cint1e_sa10nucsp",
    "cint1e_sa01sp",
    "cint1e_spgsp",
    "cint1e_spgnucsp",
    "cint1e_spgsa01",
    "cint1e_ipovlp",
    "cint1e_ipkin",
    "cint1e_ipnuc",
    "cint1e_iprinv",
    "cint1e_ipspnucsp",
    "cint1e_ipsprinvsp",
    "cint2e",
    "cint2e_spsp1",
    "cint2e_spsp1spsp2",
    "cint2e_srsr1",
    "cint2e_srsr1srsr2",
    "cint2e_sa10sp1",
    "cint2e_sa10sp1spsp2",
    "cint2e_g1",
    "cint2e_spgsp1",
    "cint2e_g1spsp2",
    "cint2e_spgsp1spsp2",
    "cint2e_ip1",
    "cint2e_ipspsp1",
    "cint2e_ip1spsp2",
    "cint2e_ipspsp1spsp2",
    "cint2e_ssp1ssp2",
    ]

#####
##### Function wrapper macros
#####

macro make_dim(name::String)
    fn = Symbol(name)
    _fn = Symbol('_'*name)
    quote
        function $(esc(fn))(
            i::Cint,
            bas::AbstractVector{Cint}
            )
            n = ccall(
                ($_fn, libcint),
                Cint,
                (Cint, Ref{Cint}),
                i, bas,
                )
        end
    end
end

macro make_fn(name::String)
    fn = Symbol(name)
    _fn = Symbol('_'*name)
    quote
        function $(esc(fn))(
            buf::AbstractArray{Cdouble},
            shls::AbstractVector{Cint},
            atm::AbstractVector{Cint},
            natm::Cint,
            bas::AbstractVector{Cint},
            nbas::Cint,
            env::AbstractVector{Cdouble},
            )
            not0 = ccall(
                ($_fn, libcint),
                Cint,
                (Ref{Cdouble}, Ref{Cint}, Ref{Cint}, Cint, Ref{Cint}, Cint, Ref{Cdouble}, Ptr{Cvoid}),
                buf, shls, atm, natm, bas, nbas, env, C_NULL,
                )
        end
    end
end

#####
##### Define functions
#####

for s in dim
    _s = Meta.parse('_'*s);
    @eval const $_s = $s
    @eval @make_dim($s)
end

for s in f1e
    _s = Meta.parse('_'*s);
    @eval const $_s = $s
    @eval @make_fn($s)
end

end # module
