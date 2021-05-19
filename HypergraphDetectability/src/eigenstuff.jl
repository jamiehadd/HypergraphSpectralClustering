function aggregateEigenvector(v, ix)
    """
    given an eigenvector v of a nonbacktracking operator and edge list ix corresponding indices of v to pointed edges, compute an aggregate eigenvector of length n. 

    The entry of the new eigenvector corresponding to i is the sum of all entries in v corresponding to edges for which i is the point. 

    this is the eigenvector actually used for community detection applications. 

    corresponds to v_in (eq. 22)  in Angelini et all 

    https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=7446987&casa_token=JyWZL63o7IkAAAAA:g5O4PX6MDgIQZp-m-kXcaL1fC663brk3DHbv6oaOe_WwjHNiRV1ZbnbBBdIEUSSVo_dKrl29eA&tag=1
    """

    n = maximum(i for (i, e) in values(ix))

    u = zeros(n)

    for edge in keys(ix)
        u[ix[edge][1]] += v[edge]
    end

    return u
end

function computeBinaryClusters(B, ix)
    E = Arpack.eigs(B; nev = 2, ritzvec = true)

    if !(imag(E[1][2]) ≈ 0)
        print("Warning: 2nd eigenvalue complex")
    end

    return aggregateEigenvector(E[2][:,2], ix)
end