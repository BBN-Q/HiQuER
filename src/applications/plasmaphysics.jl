function Vlasov_JordanWiger_matrix(N, k, α, ν)
    H = zeros(ComplexF64, (N,N))
    
    H[1,2] = k*sqrt(0.5*(1 + α))
    H[2,1] = k*sqrt(0.5*(1 + α))
    H[2,3] = k
    
    for ii = 3:N-1
        H[ii,ii-1] = k*sqrt(0.5*(ii-1))
        H[ii,ii+1] = k*sqrt(0.5*(ii))
        H[ii,ii] = -1im*(ii-1)*ν
    end
    
    H[N, N-1] = k*sqrt(0.5*(N - 1))
    H[N-1, N] = k*sqrt(0.5*(N - 1))
    H[N,N] = -1im*(N-1)*ν
    
    return H
end
    
function xx_op(N, a, b)
    xvec = zeros(Bool, N)
    xvec[a] = true
    xvec[b] = true
    return PauliOperator(0x0, xvec, zeros(Bool,N))
end

function yy_op(N, a, b)
    xvec = zeros(Bool, N)
    xvec[a] = true
    xvec[b] = true
    return PauliOperator(0x0, xvec, xvec)
end

function Vlasov_Hamiltonian_Paulis(N, k, α)
    pauli_dict = Dict{PauliOperator, Float64}()
    
    pauli_dict[xx_op(N, 1, 2)] = 0.5*k*sqrt(0.5*(1 + α))
    pauli_dict[yy_op(N, 1, 2)] = 0.5*k*sqrt(0.5*(1 + α))
    
    for ii = 2:N-1
        pauli_dict[xx_op(N, ii, ii+1)] = 0.5*k*sqrt(0.5*ii)
        pauli_dict[yy_op(N, ii, ii+1)] = 0.5*k*sqrt(0.5*ii)
    end
    
    return pauli_dict
end

function restrict_paulis(pd, K; seed=2022)
    perm = randperm(MersenneTwister(seed), K)
    
    new_pd = Dict()
    for (k,v) in pd
        new_pd[k[perm]] = v[perm]
    end
    return new_pd
end

export Vlasov_Hamiltonian_Paulis
export restrict_paulis