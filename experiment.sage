load("LWE.sage")
load("utils.sage")

n = 30
m = n
q = 3301
sigma = sqrt(20)
D_s = build_Gaussian_law(sigma, 50)
D_e = D_s
d = m + n
lwe_instance = LWE(n, q, m, D_e, D_s)
b, s, A, e_vec = (lwe_instance.b, lwe_instance.s, lwe_instance.A, lwe_instance.e_vec)
c = (-b + s * A.T + e_vec) / q
u = concatenate(lwe_instance.e_vec, c)

print(f"{u}")
print(f"{s * A.T - b - c * q}")
print(f"{-e_vec}")
# sA^T + e = b mod q
# sA^T + e - cq = b
# sA^T - cq = b

# var[c] = (sA - b) / q

# mu is the mean, right now we think it should be -b/q, though might be
# all 0s as well
sigma_c = []
mu = concatenate([0] * n, -lwe_instance.b / q)

for i in range(0, n):
    # sum of the entries in the column squared, square q as well
    # Maybe the old derivation of the sigma below is wrong, in partiuclar
    # worried about missing Var[e]
    #sigma_c.append(((sigma / q) ** 2) * lwe_instance.A.column(i).norm() ** 2)
    sigma_c.append(((sigma / q) ** 2) * (lwe_instance.A.column(i).norm() ** 2 + 1))

Sigma = block_matrix([[(sigma ** 2) * identity_matrix(n), zero_matrix(n)],
                      [zero_matrix(n), diagonal_matrix(sigma_c)]])

# Still need to figure out lattice basis, then can construct dbdd

# The basis here is just the integers, which is why we use the identity matrix
manual_ebdd = EBDD(identity_matrix(d), Sigma, mu, lwe_instance, u)
#dbdd = lwe_instance.embed_into_DBDD()
beta, delta = manual_ebdd.attack()

#print(f"beta: {beta}\n{delta=}")

for i in range(0, 30):
    v = vec([randint(0, 1) for i in range(m + n)])
    l = manual_ebdd.leak(v)
    #check whether signs here are correct
    #full dim hint here
    #convert s||e to s||c here
    _ = manual_ebdd.integrate_perfect_hint(v, l)

#check ellipsoid norm at end as well

beta, delta = manual_ebdd.attack()
print(f"beta: {beta}\n{delta=}")
