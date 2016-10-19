using SemidefiniteProgramming

sdp = SparseSDP()

setobj!(sdp, 1, 1, 2, 1.0)

setcon!(sdp, 1, 1, 1, 1, 1.0)
setrhs!(sdp, 1, 1.0)

setcon!(sdp, 2, 1, 2, 2, 1.0)
setrhs!(sdp, 2, 1.0)

sol = solve(sdp, CSDP())

X = primalmatrix(sol)
println(X)
@test_approx_eq_eps(X[1,1,2], 1, 1e-6)
@test_approx_eq_eps(X[1,2,2], 1, 1e-6)
@test_approx_eq_eps(X[1,1,1], 1, 1e-6)
