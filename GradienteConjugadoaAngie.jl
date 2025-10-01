### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 0b212dbb-d3ad-4349-b9c4-ba8bcfe619c8
begin
	using PlutoUI
	using HypertextLiteral
	using LinearAlgebra
    using LinearAlgebra
    using SparseArrays
end

# ╔═╡ 3c5944d0-866b-11ee-2f96-a71839ef3970
PlutoUI.TableOfContents(title="Gradiente Conjugado", aside=true)

# ╔═╡ 6d279e7e-f312-439e-a34c-1754621d1c28
md""" # Gradiente Conjugado"""

# ╔═╡ cf5bf699-1e78-4b71-8032-ff6b7f4ffd39
md"""## Motivación"""

# ╔═╡ 10ac6672-6f3c-4bce-a602-074aac1fa37a
md""" El método del gradiente conjugado se destaca en la resolución de sistemas lineales con matrices que son definidas positivas y simétricas, esto pues tiene una alta eficiencia con matrices de dimensión grande, no ocupa tanta memoria  y además, a comparación de otros métodos iterativos este tiene una velocidad de convergencia alta."""

# ╔═╡ fa1547fb-9e88-4522-8d6f-ac15f7705dc9
md"""##### *Preliminares:*

Algunos de los conceptos preliminares que se necesitan para entender el método son:
* **Problema de minimización:** 

    Si $A\in n\times n$ es una matriz simétrica definida positiva y $b\in \mathbb{R}^{n}$ entonces la función 

    $$\phi(x)=\frac{1}{2}x^{T}Ax-b^{T}x$$ 

    tiene un mínimo en el vector $x\in\mathbb{R}^{n}$ dado por $x=A^{-1}b$.


* **Método del máximo descenso:** Es uno de los métodos iterativos para encontrar el mínimo de una función, en este caso aplicandolo a la anterior función nos da como resultado el suguiente algoritmo, 

   **Algoritmo 1:** Escoger cualquier x' inicial.

     (0) Calcular $r=b-Ax'$ y $p=Ar$ 

     Hasta convergencia ($r<\epsilon$) itere, 

     (1) $\alpha=\frac{r^Tr}{r^Tp}$

     (2) $x'=x'+\alpha p$

     (3) $r=r-\alpha p$

     (4) $p=Ar$

* **Nueva dirección de descenso:** Tomando direcciones distintas a los residuos $\{p_1,p_2,\cdots,p_k\}$ se obtiene que el algoritmo anterior puede ser trasnformado en el siguiente

    **Algoritmo 2:** Escoger cualquier x' inicial.

     (0) Calcular $r=b-Ax'$ y $p=Ar$ 

     Hasta convergencia ($r<\epsilon$) itere,

     (1) escoja $p$ tal que $p^Tr\neq0$

     (2) $\alpha=\frac{p^Tr}{p^TAp}$

     (3) $x'=x'+\alpha p$

     (4) $r=b-Ax$

     Se observa que $x_k=x_0+\alpha_1p_1+\cdots+\alpha_kp_k$, este es un resultado importante para lo que sigue.

"""

# ╔═╡ 0616960e-af16-4202-8a76-dd7e4a85dc11
md"""## Direcciones A-conjugadas"""

# ╔═╡ f5e9679e-8967-4392-9b80-48a1910fb27d
md""" Como pudimos darnos cuenta es necesario crear las direcciones $\{p_1,p_2,\cdots,p_k\}$ que eviten los problemas del método del máximo descenso tomando como dirección los residuos, entonces tomando $x_k\in x_o+\text{gen}\{p_1+\cdots,p_k\}$, se tiene que cada $x_k$ nuevo es el resultado de realizar el siguiente problema de minimización:

$$\text{min}_{x_k}\phi(x_k)$$

Escribiendo $x_k=x_0+P_{k-1}y+\alpha p_k$ se obtiene que

$$\phi(x_k)=\phi(x_0+P_{k-1}y)+\alpha y^{T}P_{k-1}^{T}Ap_{k}+\frac{\alpha}{2}p_{k}^{T}Ap_{k}-\alpha p_{k}^{T}r_0$$

Observe que si el término $\alpha y^{T}P_{k-1}^{T}Ap_{k}$ es 0 entonces solo se necesita realizar dos minimizaciones para obtener la respuesta.

$$\text{min}_{x_k}\phi(x_k)=\text{min}_{y}\left(\phi(x_0+P_{k-1}y)\right)+\text{min}_{\alpha}\left(\frac{\alpha}{2}p_{k}^{T}Ap_{k}-\alpha p_{k}^{T}r_0\right)$$

La solución de la primera minimización sería igual al $y_{k-1}$ que hace que $x_{k-1}=x_0+P_{k-1}y_{k-1}$ y la solución del segundo problema claramente es $\alpha=\frac{p_k^{T} r_0}{p_{k}^{T}Ap_{k}}$, en especial vemos que se obtiene $\alpha=\frac{p_k^{T} r_{k-1}}{p_{k}^{T}Ap_{k}}$.

* **Existencia:**  Si $r_{k-1}\neq 0$ entonces existe $p_k\in \text{gen}\{Ap_1,\cdots,Ap_k\}^{\perp}$ tal que $p_{k}^{T}r_{k-1}\neq 0$.

Así, queda garantizada la existencia de $\{p_1,p_2,\cdots,p_k\}$ tales que $p_{k}^{T}r_{k}\neq 0$ para el algoritmo $2$ y además son $A-$conjugados.
"""

# ╔═╡ cf41a99d-9953-4a95-b5f6-86626653b4e2
md""" ## Método del gradiente conjugado"""

# ╔═╡ 17f6a73f-c2bd-4d1c-98a6-609e10dedaf4
md""" Una manera de combinar el método del máximo descenso y las direcciones $A$ conjugadas es conocida como el método del gradiente conjugado y consiste en buscar los vectores como se enseño anteriormente, pero que tengan la menor norma posible con los residuos $r_{k}$.

* **Algoritmo inicial gradiente conjugado:**
    
   Escoger cualquier $x_0$ inicial.

     (0) Calcular $r_0=b-Ax_0$, si $r_0\neq 0$ entonces;
     (1) p_{1}=r_0

     Hasta convergencia ($r_k<\epsilon$) itere k=$1,2,\cdots$,

     Si $k=1$ tome $p_{1}=r_0$

     De otra manera escoja $p_{k}$ tal que $p_{k}$ minimice $\|p-r_{k-1}\|$ sobre $p\in\{Ap_1,\cdots,Ap_{k-1}\}^{\perp}$

     (2) $\alpha=\frac{p_{k}^Tr_{k-1}}{p_{k}^TAp_{k}}$

     (3) $x_{k}=x_{k-1}+\alpha_{k} p_{k}$

     (4) $r_{k}=b-Ax_{k}$

"""

# ╔═╡ 73dc0fbe-7794-4e24-ba07-245e7ab0d9bf
md""" ## ¿Como calcular $p_k$?"""

# ╔═╡ 788086ac-f7da-42dd-ad2d-dd8c317a195e
md""" * **Lema:** Los vectores $p_k$ del algoritmo inicial del GC satisfacen que 

    $$p_{k}=r_{k-1}-AP_{k-1}z_{k-1}$$ 

   Para $z_{k-1}$ la solución del problema $\text{min}_{z}\|r_{k-1}+AP_{k-1}z\|_{2}$.

De esta manera el espacio generado por las direcciones $\{p1,\cdots,p_k\}$ es igual al espacio generado por $\{r_0,\cdots,r_{k-1}\}$ que es en realidad el **subespacio de Krylov** $\mathcal{K}_{k}(r_0,A)$.

*Demostración:*

Observe que $x_{k}=x_{k-1}+\alpha_kp_{k}$ entonces multiplicando por $-A$ y sumando $b$ a ambos lados de la igualdad resulta que $r_{k}=r_{k-1}+A\alpha_{k}p_{k}$, de aquí se concluye que $\{Ap_1,\cdots,Ap_k\}\subset \text{gen} \{r_0,r_1,\cdots,r_k\}$. Ahora, por el lema se tiene entonces que $p_{k}\in\text{gen} \{r_0,r_1,\cdots,r_k\}$ entonces 

$$[p_1,p_2,\cdots,p_k]=[r_0,r_1,\cdots,r_k]T$$

En donde $T$ es una matriz triangular superior que además es no singular pues $[p_1,p_2,\cdots,p_k]$ tiene rango completo. De aquí $\text{gen}\{p_1,p_2,\cdots,p_k\}=\text{gen}\{r_0,\cdots,r_{k-1}\}$ ahora este último espacio es igual al subespacio de Krylov pues $r_{k}\in \text{gen}\{r_{k-1},Ap_k\}\subset\text{gen}\{r_{k-1},Ar_0,\cdots,Ar_{k-1}\}$, esto úlimo reemplazando la expresión dada para $p_k$, así reemplazando sucesivamente se llega a tal espacio.

Así, ya sabemos como calcular las direcciones necesarias para el algoritmo.

** Teorema:**

Las direcciones $p_{k}$ $A-$ conjugas para $k\geq 2$ tienen la siguiente propiedad y es que $p_{k}\in \text{gen}\{p_{k-1},r_{k-1}\}$.

"""

# ╔═╡ b65a3cb0-bcc5-458e-8723-4beae6299aba
md""" ## Algoritmo final"""

# ╔═╡ c0ad9a30-2360-474a-a61a-7c7cb7f0a8c4
md""" De lo anterior se pueden deducir ciertas cosas, teniendo en cuenta que $p_{k}=r_{k-1}+\beta_{k}p_{k-1}$.

* Como $p_{k}^{T}Ap_k=0$ entonces $\beta_{k}=\frac{p_{k-1}^{T}Ar_{k-1}}{p_{k-1}^{T}Ap_{k-1}}$.

* $\alpha=\frac{r_{k-1}^{T} r_{k-1}}{p_{k}^{T}Ap_{k}}$.

* $r_{k-1}^{T}r_{k-1}=-\alpha_{k-1}r_{k-1}^{T}Ap_{k-1}$.

* $r_{k-2}^{T}r_{k-2}=\alpha_{k-1}p_{k-1}^{T}Ap_{k-1}$.

Transformando a $\beta$ en la siguiente expresión: 

 $$\beta_{k}=\frac{r_{k-1}^{T}r_{k-1}}{r_{k-2}^{T}r_{k-2}}$$

Entonces el algoritmo es el siguiente: 

* **Algoritmo CG:**

     Escoger cualquier $x_0$ inicial 

     (0) Calcular $r_0=b-Ax_0$, si $r_0\neq 0$ entonces;

     Hasta convergencia ($r_k<\epsilon$) itere k=$1,2,\cdots$,

     Si $k=1$ tome $p_{1}=r_0$

     De otra manera $\beta_{k}=\frac{r_{k-1}^{T}r_{k-1}}{r_{k-2}^{T}r_{k-2}}$ y $p_{k}=r_{k-1}+\beta_{k}p_{k-1}$.

     (2) $\alpha=\frac{r_{k-1}^Tr_{k-1}}{p_{k}^TAp_{k}}$

     (3) $x_{k}=x_{k-1}+\alpha_{k} p_{k}$

     (4) $r_{k}=r_{k-1}-\alpha_{k-1}Ap_{k}$

     Fin

"""

# ╔═╡ 30fe81c0-6beb-4630-9ed8-7d482d8f5ebd
md""" ## Implementación"""

# ╔═╡ ed27cc82-6af5-4418-95d8-41388e1805d3
begin
function  cg(A, x, b, max_it, tol)
# Translate to julia from
#  -- Iterative template routine --
#     Univ. of Tennessee and Oak Ridge National Laboratory
#     October 1, 1993
#     Details of this algorithm are described in "Templates for the
#     Solution of Linear Systems: Building Blocks for Iterative
#     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
#     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
#     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
#
#  [x, error, iter, flag] = cg(A, x, b, M, max_it, tol)
#
# cg.m solves the symmetric positive definite linear system Ax=b 
# using the Conjugate Gradient method with preconditioning.
#
# input   A        REAL symmetric positive definite matrix
#         x        REAL initial guess vector
#         b        REAL right hand side vector
#         M        REAL preconditioner matrix
#         max_it   INTEGER maximum number of iterations
#         tol      REAL error tolerance
#
# output  x        REAL solution vector
#         error    REAL error norm
#         iter     INTEGER number of iterations performed
#         flag     INTEGER: 0 = solution found to tolerance
#                           1 = no convergence given max_it

    flag = 0 # initialization
    iter = 0
    lastiter=max_it
    α=fill(0.0,(1,max_it))
    β=fill(0.0,(1,max_it))
    bnrm2 = norm( b )
    if  ( bnrm2 == 0.0 )
        bnrm2 = 1.0
    end

    r = b - A*x
    error = norm( r ) / bnrm2
    println(" PCG residual ",iter," = ",error)
    if ( error < tol ) 
        return
    end

    for iter = 1:max_it                       # begin iteration

        z  = r
        #    z  = pccg(r,iter)
        rho = (r'*z)
        rho_1=rho
        p=z
        if ( iter > 1 )                       # direction vector
            β[iter] = rho / rho_1
            p = z + β[iter]*p
         else
            p = z
         end

         q = A*p
         α[iter] = rho / (p'*q )
         x = x + α[iter] * p                  # update approximation vector

         r = r - α[iter]*q                    # compute residual
         error = norm( r ) / bnrm2            # check convergence
        println("PCG residual (",iter,") = ",error)
     if ( error <= tol )
            lastiter=iter
            break
    end 



  end
#  compute eigenvalues and codition number
#
    d=sparse(fill(0.0, (lastiter,1)))
     d[1]= 1/α[1];
     for i = 2: iter
         d[i]=β[i]/α[i-1]+1/α[i];
     end
     for i = 1: iter-1
         s[i]=-1*sqrt(β[i+1])/α[i];
     end
# 
     T = sparse(fill(0.0, (lastiter,lastiter)))
     T[1,1]= d[1];
     for i=2:iter
         T[i,i] = d[i]
     end
     for i = 1:iter-1
         T[i,i+1]= s[i]; T[i+1,i] = T[i,i+1];
     end
     lambda = eigvals(Matrix(T))
     lambdamax = maximum(lambda)
     lambdamin = minimum(lambda)
     condnumber = lambdamax/lambdamin

    if ( error > tol ) 
        flag = 1
    end         # no convergence
    return x, error, iter, flag,condnumber
end
end

# ╔═╡ 50b987d4-aac8-4fbc-915b-139caa25152f
begin
n=10
A=rand(n,n)
A=A'A+0.1*I
b=fill(1,n)
x₁=b*0.1;
end

# ╔═╡ 75e28b43-8d1d-442c-821a-46d252180638
begin
x, error, iter, flag,condnumber=cg(A,x₁,b,1000,0.0001);
end

# ╔═╡ 28815c8f-7a75-4702-b6a8-dc4a187faf05
md""" Solución:"""

# ╔═╡ ecc3c91e-7ca7-4597-a726-153c719b9414
x

# ╔═╡ 0bf177c0-8074-403f-8258-b8ddf2601717
md"""Cuando el vector se acerca mucho a 0:"""

# ╔═╡ f000f082-14ca-4a8b-9522-e9bd2287ae01
begin
x₃=0.000000001*b
Y2, error2, iter2, flag2,condnumber2=cg(A,x₃,b,1000,0.0001);
end

# ╔═╡ bc7daae7-c316-4ea8-9f9c-7ae54ff3292b
begin
x₄=0*b
Y3, error3, iter3, flag3,condnumber3=cg(A,x₄,b,1000,0.00001);
end

# ╔═╡ 13ca7e9f-8165-4b5c-acb4-9c6da0a4d878
md""" ## Conexión con DIOM y LANCZOS"""

# ╔═╡ d3725e91-940b-4cf8-ab00-0e7ec2262323
md""" Ya sabemos que las direcciones necesarias para el algoritmo del gradiente conjugado están dadas precisamente por el **subespacio de Krylov** mediante la teoría, sin embargo también es importante conocer que el algoritmO GC se puede derivar de los siguientes dos algoritmos.

Primero, observe que el algoritmo de LANCZOS es solo el algoritmo de Arnoldi con Gram schmidt modificado para matrices simétricas.

* **Algoritmo Lanczos:** 

    Escoja un vector $v_1$ de norma $1$ y tome $\beta_1=0$ y $v_0=0$

    Para $j=1,2,\cdots,k$ haga:

    (1) $w_j= Av_j-\beta_jv_{j-1}$

    (2) $\alpha_j=(w_j,v_j)$

    (3) $w_j=w_j-\alpha_jv_j$

    (4) $\beta_{j+1}:= \|w_j\|_{2}$ Si $\beta_{j+1}=0$ pare.

    (5) $v_{j+1}=\frac{w_j}{\beta_{j+1}}$

Ahora, si aplicamos este método y lo extendemos a sistemas de ecuaciones lineales como se realizo con arnoldi, se tiene que la matriz resultante es una matriz tridiagonal en vez de una matriz de Hessenberg, y que además la solución del sistema de eduaciones lineales está dado por:

$$x_k=x_o+V_{k}y_{k} \text{ en donde } y_{k}=T^{-1}_{k}(\beta e_1^{T})$$

* **Algoritmo Lanczos para sistemas de eduaciones lineales:**

    (0) Escoja un vector $x_0$ cualquiera y calcule $r_0=b-Ax_0$, $\beta=\|r_0\|_{2}$ y $v_1=\frac{r_0}{\beta}$.

    Para $j=1,2,\cdots,k$ haga:

    (1) $w_j= Av_j-\beta_jv_{j-1}$ (Si j=1 deje $\beta_1v_0=0$)

    (2) $\alpha_j=(w_j,v_j)$

    (3) $w_j=w_j-\alpha_jv_j$

    (4) $\beta_{j+1}:= \|w_j\|_{2}$ Si $\beta_{j+1}=0$ tome $k=j$ y vaya a (6).

    (5) $v_{j+1}=\frac{w_j}{\beta_{j+1}}$

    Fin 

    (6) Forme $T_{k}=\text{tridiag}(\beta_i,\alpha_j,\beta_{i+1})$ y $V_{k}=[v_1,\cdots,v_{k}]$.

    (7) Calcule $y_k=T_{k}^{-1}(\beta e_1)$ y $x_k=x_0+V_{k}y_{k}$


La obtención de $y_k$ sabemos que no es realizado con la matriz inversa, luego nos obliga a crear otro algoritmo que se llama $D-lanczos$ que se deriva de DIOM, en donde se utilizá una matriz $P_{k}$ que está dada por $V_{k}U_{k}^{-1}$ siendo $U$ la matriz triangular superior de la descomposición $LU$ de $T_{k}$, entonces se siguen los siguientes pasos:

* **Algoritmo D-lancsoz:**

    (0) Escoja un vector $x_0$ cualquiera y calcule $r_0=b-Ax_0$, $\gamma_1=\beta=\|r_0\|_{2}$ y $v_1=\frac{r_0}{\beta}$.

    (1) Tome $\lambda_1=0=\beta_1$ y $p_0=0$. 

    Para $j=1,2,\cdots$ hasta que converja haga:

    (2) $w= Av_j-\beta_jv_{j-1}$ y $\alpha_j=(w,v_j)$

     Si $k>1$ entonces calcule $\lambda_j=\frac{\beta_j}{\eta_{j-1}}$ y $\gamma_j=-\lambda_j\eta_{k-1}$

    (3) $\eta_j=\alpha_j-\lambda_j\beta_j$

    (4) $p_j=(\eta_j)^{-1}(v_j-\beta_jp_{j-1}$

    (5) $x_{k}=x_{k-1}-\gamma_jp_j$ Si $x_k$ converge entonces pare.

    (6) $w=w-\alpha_j v_j$

    (7) $\beta_{j+1}=\|w\|_{2}$ y $v_{k+1}=\frac{w}{\beta_{j+1}}$
    Fin
     
"""

# ╔═╡ 665777ad-3fca-45f1-90e4-53ec1e49b2a5
md""" No es dificil ver que los $p_{k}$ que se hallan en el algoritmo de $D-lancsoz$ son A-conjugados y de aquí se deriva de nuevo el método GC."""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[compat]
HypertextLiteral = "~0.9.4"
PlutoUI = "~0.7.52"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.3"
manifest_format = "2.0"
project_hash = "bce32509f8f2ab738bd97a8845c1f603f31f945a"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "91bd53c39b9cbfb5ef4b015e8b582d344532bd0a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e47cd150dbe0443c3a3651bc5b9cbd5576ab75b7"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.52"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "9673d39decc5feece56ef3940e5dafba15ba0f81"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.1.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.URIs]]
git-tree-sha1 = "b7a5e99f24892b6824a954199a45e9ffcc1c70f0"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═3c5944d0-866b-11ee-2f96-a71839ef3970
# ╠═0b212dbb-d3ad-4349-b9c4-ba8bcfe619c8
# ╟─6d279e7e-f312-439e-a34c-1754621d1c28
# ╟─cf5bf699-1e78-4b71-8032-ff6b7f4ffd39
# ╟─10ac6672-6f3c-4bce-a602-074aac1fa37a
# ╟─fa1547fb-9e88-4522-8d6f-ac15f7705dc9
# ╟─0616960e-af16-4202-8a76-dd7e4a85dc11
# ╟─f5e9679e-8967-4392-9b80-48a1910fb27d
# ╟─cf41a99d-9953-4a95-b5f6-86626653b4e2
# ╟─17f6a73f-c2bd-4d1c-98a6-609e10dedaf4
# ╟─73dc0fbe-7794-4e24-ba07-245e7ab0d9bf
# ╟─788086ac-f7da-42dd-ad2d-dd8c317a195e
# ╟─b65a3cb0-bcc5-458e-8723-4beae6299aba
# ╟─c0ad9a30-2360-474a-a61a-7c7cb7f0a8c4
# ╠═30fe81c0-6beb-4630-9ed8-7d482d8f5ebd
# ╠═ed27cc82-6af5-4418-95d8-41388e1805d3
# ╠═50b987d4-aac8-4fbc-915b-139caa25152f
# ╠═75e28b43-8d1d-442c-821a-46d252180638
# ╟─28815c8f-7a75-4702-b6a8-dc4a187faf05
# ╠═ecc3c91e-7ca7-4597-a726-153c719b9414
# ╟─0bf177c0-8074-403f-8258-b8ddf2601717
# ╠═f000f082-14ca-4a8b-9522-e9bd2287ae01
# ╠═bc7daae7-c316-4ea8-9f9c-7ae54ff3292b
# ╟─13ca7e9f-8165-4b5c-acb4-9c6da0a4d878
# ╟─d3725e91-940b-4cf8-ab00-0e7ec2262323
# ╟─665777ad-3fca-45f1-90e4-53ec1e49b2a5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
