### A Pluto.jl notebook ###
<<<<<<< Updated upstream
# v0.19.45
=======
# v0.20.4
>>>>>>> Stashed changes

using Markdown
using InteractiveUtils

# ╔═╡ 9b0a076b-a004-422b-be3d-9fd08ce2f966
begin
	using PlutoUI
	using HypertextLiteral
	using LinearAlgebra
end

# ╔═╡ 7494046b-46d9-4233-b99a-adaedbd52f04
PlutoUI.TableOfContents(title="Descomposición de Schur", aside=true)

# ╔═╡ d84ac4f0-4e6b-11ee-387e-a3d5d669b309
md""" #  Descomposición de Schur"""

# ╔═╡ 4e4da5e9-a81a-47f7-a0ce-7e56b79d1314
md""" ### Motivación. """

# ╔═╡ b741b172-413b-4ac6-98af-c3399e55a263
md""" Cuando se tiene una matriz $A\in \mathbb{C}^{n\times n}$ se pueden definir distintos conceptos, entre ellos:

* Se define a los autovalores de la matriz $A$ como los escalares $\lambda_1,\cdots \lambda_n$ que son ceros del polinomio dado por $p(x)=|A-\lambda I|$, adicionalmente al multiconjunto de dichos autovalores lo notamos como $\lambda(A)$. 

* Si $\lambda\in\lambda(A)$ existe un vector  $v$ tal que $v\neq 0$ y $Av=\lambda v$, a este vector se le llama autovector de la matriz $A$. 

Es así que si nuestra matriz $A$ tiene $n$ autovectores linealmente independientes se dice que $A$ es diagonalizable, es decir, $A$ es similar a una matriz diagonal $D$. En otras palabras

$$A=XDX^{-1}$$

en donde $X$ es una matriz que en cada una de sus columnas tiene a un autovector de $A$. Se puede observar que descomponer una matriz de esta manera tiene bastantes ventajas a la hora de hallar su exponencial.

Por ello mismo, surge la pregunta ¿Que ocurre con las matrices que no son diagonalizables? ¿Se puede tener una descomposición que tenga precisamente este tipo de ventajas?"""


# ╔═╡ 3c3a0d73-d062-4f4b-a858-6c65a8a2b620
md"""##### *Preliminares:*


**Lema 1:** Si $T\in\mathbb{C}^{n\times n}$ se particiona de la forma siguiente

$$T=\begin{bmatrix} T_{11}&T_{12}\\ O& T_{22} \end{bmatrix}$$ 

entonces $\lambda(T)=\lambda(T_{11})\cup \lambda(T_{22})$.

**Lema 2:** Si $A\in\mathbb{C}^{n\times n}$, $B\in\mathbb{C}^{p\times p}$ y $X\in\mathbb{C}^{n\times p}$ son tales que satisfacen que $AX=XB$ y $rang(X)=p$ entonces existe una matriz $Q\in\mathbb{C}^{n\times n}$ unitaria tal que

$$Q^{H}AQ=T=\begin{bmatrix} T_{11}&T_{12}\\O &T_{22} \end{bmatrix}$$

donde $\lambda(T_{11})=\lambda(A)\cap\lambda(B)$.
 """

# ╔═╡ 9a24ac91-2b77-4bb5-b6e5-a070cabc8e35
md"""## Teorema. 

Si $A\in\mathbb{C}^{n\times n}$ existe una matriz $Q\in\mathbb{C}^{n\times n}$ unitaria ral que

$$A=QTQ^{H}$$ 

en donde $T=D+N$ siendo $D=diag(\lambda_1,\lambda_2,...,\lambda_n)$ y $N$ es una matriz $n\times n$ estrictamente triangular superior. 

"""

# ╔═╡ f13281a5-7785-4ada-a65b-c31c48fdd6ad
md""" ### Demostración. """

# ╔═╡ 589cce26-775d-4a13-b432-b4bd81ae855f
md""" Se demostrará por inducción sobre $n$:

* Para $n=1$ es evidente que se cumple la existencia de dicha descomposición.
* Supongamos que toda matriz $A'$ en los complejos de tamaño $k$ para $1\leq k\leq n-1$ tiene una descomposición de Schur.
* Sea $A\in \mathbb{C}^{n\times n}$  y sea $v$ un autovector de $A$, además sea $\lambda$ el autovalor asociado a $v$, usando el lema 2 para $B=\lambda$ se tiene que existe una $\hat{Q}$ unitaria tal que

$$\hat{Q}^{H}A\hat{Q}=\begin{bmatrix} T_{11} &T_{12}\\O &T_{22}\end{bmatrix}$$

$\hspace{12mm}$ Pero $T_{11}=\lambda$ entonces realmente se tiene que

$$\hat{Q}^{H}A\hat{Q}=\begin{bmatrix} \lambda & w^{H}\\O & C\end{bmatrix}$$

Con anterioridad sabemos que $C$ tiene una descomposición de schur, sea $\tilde{Q}$ unitaria tal que $\tilde{Q}^{H}C\tilde{Q}$ es triangular, por tanto se puede tomar $Q=\hat{Q}diag(1,\tilde{Q})$ y se sigue que $Q^{H}AQ$ es triangular superior.

Así queda demostrada la existencia de una descomposición de schur para cualquier matriz cuadrada en los complejos.


"""

# ╔═╡ ef2798f3-0255-4039-aaa0-0a09b0f25f76
md"""### Consecuencias.

* Tomando $S_{k}=\langle q_1,q_2,...,q_k\rangle$ se puede ver que este subespacio es invariante con respecto a $A$.

* De igual manera que en la diagonalización se puede hallar la potencia de una matriz de una manera mucho m¿as facil.

"""

# ╔═╡ 6f12b4e3-4c3d-40d6-adec-7f78d7f1b04f
md""" ### Algoritmo."""

# ╔═╡ 01588361-4b3d-4ebe-91af-4d86de9c06b2
md""" Para calcular la descomposición de Schur de una matriz es necesario realizar la factorización QR de la matriz, en cada iteración:"""

# ╔═╡ 80810067-023b-48e9-a4cf-fe458fd90189
function QRCGS(A)
    sizeA=size(A)
    Q = zeros(sizeA) #(m,n)
    R = zeros(sizeA[2],sizeA[2]) #(n,n)
    for i = 1:sizeA[2] #
        for j = 1:i-1
            R[j,i] = Q[:,j]'A[:,i]
        end
        p = A[:,i] - Q[:,1:i-1]*R[1:i-1,i]
        R[i,i]=norm(p)
#        if abs(R[i,i])<0.00000000001
#            println("Rii cercano 0")
#        end
        Q[:,i] = p/R[i,i]
    end
    return Q,R
end

# ╔═╡ 290f6613-3b20-4658-b6ff-28ee19a959c4
md"""La descomposición de schur es:"""

# ╔═╡ 06b62b7f-01a5-45c7-9abe-b86ac44f51ed
begin
 function Schur(A)
	A₁=copy(A)
	 Q₁=UniformScaling(1);
	for i=1:1000
	Q,R=QRCGS(A₁)
	Q₁=Q₁*QRCGS(A₁)[1]
	A₁=R*Q
#    if(opnorm(A-(Q₁*A₁*adjoint(Q₁)))<epsilon)
#	else 
#		print("El método no converge")
#	end
	end

 	return A₁,Q₁

 end


end


# ╔═╡ 892365bc-0464-4d5a-be36-6f938379c58f
md""" por ejemplo, escojamos la siguiente matriz y hallemos su descomposición de Schur por medio del algoritmo:"""

# ╔═╡ 9f5b7c35-7a7d-4a50-9e20-3cae5c55eef1
begin
B=rand(4:10,5,5)
epsilon=1E-10
A₂, Q₂ =Schur(B); display(A₂)
end

# ╔═╡ 4884e5b2-4bdc-4ae1-96af-336df91aeede
eigvals(B)

# ╔═╡ 582b64d5-f7fa-4477-98fc-36d08e8b9efa
Block=[-0.470364     -1.81725 ; 1.46617      -0.72121]

# ╔═╡ 4d982c26-a0ef-407a-9f19-8d876502697f
eigvals(Block)

# ╔═╡ 28fa22fc-c9ba-4da3-882c-1a386345c962
0.470364/1.81725

# ╔═╡ 6fcf1cdd-5780-430d-bdd2-6934370d4fcd
display(Q₂)

# ╔═╡ c00dc9e0-9970-41ab-9aae-b39205df9f66
md""" Se puede ver que el error de descomposición esta dado por lo siguiente:"""

# ╔═╡ 0cf97d60-8d02-48f8-94af-bbb91fee69d6
opnorm(B-Q₂*A₂*adjoint(Q₂))

# ╔═╡ 756ae63a-a811-4cd2-b888-5d946a1e9527
md""" Además, el error sobre Q para ser una matriz unitaria es:"""

# ╔═╡ 8d892dc3-81a3-467b-89fe-8640589451d7
opnorm(UniformScaling(1)-Q₂*adjoint(Q₂))

# ╔═╡ 5613492b-90fe-4f44-87bb-91e2edf9a960
md""" En julia se encuentra implementada la siguiente función que realiza la descomposición, con ella se puede comparar que tan bueno es el resultado anteriormente obtenido."""

# ╔═╡ 1869882b-95fd-417e-9c53-961abc175cc6
T, QU=schur(B) #Función de Julia que cálcula la descomposición

# ╔═╡ 8edcd362-5eac-4e6d-bdec-85b7d71479f4
md""" Comparación entre los algoritmos"""

# ╔═╡ 4b06228f-ebe2-429a-884a-36df923995dc
opnorm(QU*T*adjoint(QU)-Q₂*A₂*adjoint(Q₂))

# ╔═╡ 919e710b-8b63-4f4e-8543-b283ce987689
md""" ¿Que tan bueno es el algoritmo de julia?"""

# ╔═╡ d7b4464b-2a88-4671-b3e1-c2c227889fc4
opnorm(QU*T*adjoint(QU)-B)

# ╔═╡ ca2cd52d-1e2c-4176-b3b8-8a8d402882f9
md"""**Aplicaciones**

* A la hora de realizar potencias de Matrices se facilita gracias a la descomposición de Schur.
* A la hora de hallar $exp(A)$.
* Para encontrar los autovalores de la matriz $A$."""

# ╔═╡ f98e7839-7318-4363-9e8e-e3750cc3a400
md""" El siguiente ejemplo nos muestra como hallar una potencia de la matriz $A_{3}$, aquí hallaremos la potencia 50. """

# ╔═╡ c8fb6db5-c387-4e04-935e-36af5982f202
begin
	A₃=rand(1000,1000)
	T₁,U=schur(A₃)
end

# ╔═╡ 0daae961-bfbc-4807-b550-23ae6c9b3695
B₅₀=U*T₁^50*adjoint(U)

# ╔═╡ c64aef5a-5153-4ffc-ac11-9470d11cb87e
md"""Comparamos el resultado con el siguiente:"""

# ╔═╡ 71989b8d-9dd0-488b-ab34-93f83ef0bf52
A₃^50

# ╔═╡ bbe49f09-9e21-49b3-9cc7-f8b2f1ae0871
md""" ### Casos Especiales."""

# ╔═╡ f3685950-6a36-45c7-9b74-ff8004d1f41b
md"""  **Teorema:** $A\in\mathbb{C}^{n\times n}$ es normal si y solo si existe una matriz unitaria $Q\in\mathbb{C}^{n\times n}$ tal que $Q^{H}AQ=diag(\lambda_1,\lambda_2,...,\lambda_n)$. """

# ╔═╡ 036b12b9-897a-46c3-8b5c-32ebf4cd8855
md""" El algoritmo anteriormente usado no es eficiente, se puede notar que al tener una matriz mucho más grande este se va a demorar demasiado, por tanto es util introducir un método que ayuda a hallar la descomposición de Schur de una manera más eficiente. Para ello se definirá una nueva versión de la descomposición de Schur. """

# ╔═╡ ee51dc26-d8fd-4d14-a1f1-b6665a6a3715
md""" # Descomposición real de Schur."""

# ╔═╡ 9ee2dd6a-378b-408f-a7b4-065837a3e5a0
md""" ## Teorema."""

# ╔═╡ a3043126-bdfb-48e4-ba90-32b1b73b25c5
md""" Sea $A\in\mathbb{R}^{n\times n}$ entonces existe una matriz $Q\in\mathbb{R}^{n\times n}$ ortogonal tal que 

$$Q^{T}AQ=\begin{bmatrix}R_{11}&R_{12}&\cdots&R_{1n}\\O&R_{22}&\cdots& R_{2n}\\
\vdots &\vdots &\vdots &\vdots\\ O&O&O&R_{nn} \end{bmatrix}$$

En donde $R_{ii}$ es una matriz $1\times 1$ o $2\times 2$ en donde se encuentran los autovalores complejos conjugados.
"""

# ╔═╡ 4383fa1e-2c6e-4c2c-a9ed-1ac28d69716f
md""" ### Demostración."""

# ╔═╡ a50265e1-5c47-4d53-8992-134a123ecb33
md""" Los autovalores complejos de A vienen en pares, el y su conjugado. Sea $k$ el número de parejas de autovalores de $A$, este teorema se demostrará realizando inducción sobre $k$.

* Si $k=0$ entonces se cumple.
* Supongamos que para $1\leq k\leq n-1$ se tiene que existe tal descomposición real de schur para una matriz de cumpla dicha condición.
* Sea $A\in \mathbb{R}^{n\times n}$ tal que $k=n$, entonces al menos existe un autovalor complejo, tómese $\lambda=a+bi$ y además observe que este autovalor tiene un autovector asociado dado por $v=x+iy$ con $x,y$ vectores en los números reales y $y\neq 0$, por esto mismo se tiene que

$$A[x\hspace{3mm} y]=(a+bi)[x\hspace{3mm} y]$$

utilizando una versión del lema 2 para matrices reales se sigue que existe $\hat{Q}$ tal que 

$$\hat{Q}^{T}A\hat{Q}=\begin{bmatrix}T_{11}&T_{12}\\O &T_{22}\end{bmatrix}$$ en donde $\lambda(T_{11})=\{\lambda,\bar{\lambda}\}$. Ahora $T_{22}$ es una matriz de tamaño $n-2$ y por tanto $k<n$, luego aplicamos hipótesis de inducción y llegamos a que existe $\tilde{Q}$ ortogonal, tales que $\tilde{Q}^{T}T_{22}\tilde{Q}$ es triangular superior. 

Así, tomando $Q=\hat{Q}diag(I_{2},\tilde{Q})$ se obtiene lo buscado."""

# ╔═╡ 23abd4e6-d32e-4937-a11f-592d7bee58de
md""" Como se dijo anteriormente, el objetivo es optimizar el número de operaciones y hacer lo más eficiente posible la obtención de la descomposición, por tanto se hablará del siguiente método:"""


# ╔═╡ 30e9ca21-7007-483b-8b0f-c6296606b653
md"""### QR Hessenberg. """

# ╔═╡ ea9cb1c7-8cab-4f1f-b7d7-fef68a5b0d89
md""" Se desea buscar una matriz  $U_o$ tal que al realizar $U_o^{T}AU_{o}$ esto sea igual a una matriz de Hessenberg superior y despues aplicar factorización QR por medio del método de Givens. Lo cual agiliza en gran medida hallar la factorización $QR$ en cada paso del algoritmo de la descomposición de schur. """

# ╔═╡ 316bf71f-03e7-43ce-ab39-2ca1decb6266
md""" **¿Como se halla la matriz $U_o$?**


Esta matriz se puede hallar por medio de multipliación de matrices Householder, para ello necesitaremos el algoritmo de householder que se encuentra a continuación."""



# ╔═╡ 51241163-7867-40f2-8c90-498aed295a3f
begin
function Housev(x)
    n = length(x)
    v = ones(size(x))
    v[2:n] = x[2:n]
    σ = norm(x[2:n])^2
    if σ == 0
        β = 0
    else 
        μ = √(x[1]^2+σ)
        if x[1] ≤ 0
            v[1] = x[1] - μ
        else
            v[1] = -σ/(x[1]+μ)
        end
        β = 2*v[1]^2/(σ+v[1]^2)
        v = v/(v[1])
    end
    return v, β
end
end

# ╔═╡ aab2d571-e794-4c35-9ff1-c96415107bce
md""" Sabemos que al aplicar matrices de Householder a la matriz $A$ podemos "eliminar" bastantes elementos de una columna, en este caso se realiza convenientemente para eliminar los elementos necesarios y así llegar a una matriz hessenberg superior, es por ello que precisamente la $U_o$ que buscamos será producto de matrices Householder. El siguiente algoritmo nos arroja la información sobre quien es la matriz de Hessenberg resultante y quien es la matriz $U_o$:"""

# ╔═╡ 4d343352-9263-4a16-8de7-586ea478a3bc
begin
	function HessenbergForm(A)
    n = size(A)[1]
    H = copy(A)
    U = Matrix(1.0*I, n, n)
    for k = 1:n-2
        v, β = Housev(H[k+1:n, k])
        H[k+1:n, k:n] = (I - β*v*v')*H[k+1:n, k:n]
        H[1:n, k+1:n] = H[1:n, k+1:n]*(I - β*v*v')
        
        #U es necesaria para la verificar que la función devuelve los resultados correctos
        U[1:n, k+1:n] = U[1:n, k+1:n]*(I - β*v*v') 
    end
    return H, U
end
end

# ╔═╡ 9141e633-6022-4c87-a2ac-296927249549
md""" Por ejemplo continuando con la matriz inicial se tiene lo siguiente:"""

# ╔═╡ 6b123991-0773-40f8-a0ac-6ea353ca0a3e
begin
 H₁,U₁=HessenbergForm(B); display(H₁)
end

# ╔═╡ 88552670-04af-4b2a-858c-11519c5471b2
display(U₁)

# ╔═╡ 8113b9b7-c5b9-4e9a-a02b-5d87ed4f60df
md""" Ahora, se tiene una matriz casi $H$ casi triangular, luego solo necesitamos quitar los números de la primera subdiagonal principal y reemplazarlos por 0's. Esto como ya lo vimos se puede realizar con Givens:"""

# ╔═╡ 040aeed9-4ed4-43c0-8c24-eecd73862ece
function Givens(a,b)
    if b==0
        c = 1
        s = 0
    else
        if abs(b)>abs(a)
            τ=-a/b
            s=-1/sqrt(1+τ^2)
            c=s*τ
        else
            τ=-b/a
            c=1/sqrt(1+τ^2)
            s=c*τ
        end
    end
    return c,s
end

# ╔═╡ 7cbc0cee-6a73-40e1-8399-86892110f8df
md"""Así la descomposición real de schur se puede calcular con el siguiente algoritmo:"""

# ╔═╡ 6a3e75d4-255c-4039-8f0f-b166adf05871
begin
	function HessenbergQR(H)
    n = size(H)[1]
    H2 = copy(H)
    C, S = zeros(n-1), zeros(n-1)
    
    #Factorización QR de H
    for k = 1:n-1
        C[k], S[k] = Givens(H2[k,k], H2[k+1, k])
        H2[k:k+1,k:n] = [C[k] -S[k]; S[k] C[k]]*H2[k:k+1,k:n]
        H2[k+1,k] = 0
        #display([C[k] S[k]])
    end
    
    #Matriz RQ
    for k = 1:n-1
        H2[1:k+1, k:k+1] = H2[1:k+1, k:k+1]*[C[k] S[k]; -S[k] C[k]]
    end
    
    return H2 #RQ Hessenberg superior
end
end

# ╔═╡ 27f987d7-7c82-47bc-a441-e9fc4469bc56
begin
	function RealSchur(A, iteraciones = 10000)
    H0 = A
    H1 = HessenbergForm(A)[1]
    δ = 10
    for k = 1:iteraciones
        H0 = H1
        H1 = HessenbergQR(H1)
    end
    return H1
    
end
end

# ╔═╡ ec201a3f-e323-49d9-8904-62241bb1434e
begin
H₂=RealSchur(B)
end

# ╔═╡ 6e8930b1-2357-4a2e-bc72-388b7dd42258
md""" Veamos que ocurre con una matriz mucho más grande"""

# ╔═╡ 0ed3770e-9909-4fcd-a190-ce71d05ef5bd
begin 
C=rand(1000,1000)
H₃=RealSchur(C)
end

# ╔═╡ 65d97516-3de2-47a8-b536-dc5430f622cc
md""" Observamos que se demora un poco pero realiza todas las iteraciones correspondientes y arroja la matriz que nos interes, al contrario de lo que ocurriria si metieramos esta misma matriz en el algoritmo presentado con QR realizado con gram schmidt"""

# ╔═╡ 9cc5da63-fdaf-4d94-8445-cd31fd0888c6
md""" Por ultimo veamos una matriz real con autovalores complejos"""

# ╔═╡ 7a84fc5c-d956-4187-98ba-87ed0e7641e2
#Esta función nos muestra los autovalores complejos que vienen por pares:

function autovComplejos(H)
    D = diag(H,-1)
    for k = 1:size(D)[1]
        if abs(D[k])>10^-12
            display(eigvals(H[k:k+1,k:k+1]))
        end
    end
end

# ╔═╡ cf6667a8-e033-4ef6-a7a3-2d788de8b6ae
begin
D = floor.(100*rand(100,100)-I*100)
display(D)
print("\n Los autovalores complejos son: \n")
display(eigvals(D))
end

# ╔═╡ e6246548-4164-4f55-bfd4-cd7d36564eee
begin
print("\nLa forma de Schur real es: \n")
SchurD = RealSchur(D)
display(SchurD)
end

# ╔═╡ 1d97d78b-859a-4d79-81f8-bbaf2b741447
begin
print("Los autovalores de los bloques 2x2 son: \n")
autovComplejos(SchurD)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
HypertextLiteral = "~0.9.4"
PlutoUI = "~0.7.52"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

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

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "d3de2694b52a01ce61a036f18ea9c0f61c4a9230"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.62"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

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

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "6cae795a5a9313bbb4f60683f7263318fc7d1505"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.10"

[[deps.URIs]]
git-tree-sha1 = "cbbebadbcc76c5ca1cc4b4f3b0614b3e603b5000"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═7494046b-46d9-4233-b99a-adaedbd52f04
# ╠═9b0a076b-a004-422b-be3d-9fd08ce2f966
# ╟─d84ac4f0-4e6b-11ee-387e-a3d5d669b309
# ╟─4e4da5e9-a81a-47f7-a0ce-7e56b79d1314
# ╟─b741b172-413b-4ac6-98af-c3399e55a263
# ╟─3c3a0d73-d062-4f4b-a858-6c65a8a2b620
# ╟─9a24ac91-2b77-4bb5-b6e5-a070cabc8e35
# ╟─f13281a5-7785-4ada-a65b-c31c48fdd6ad
# ╟─589cce26-775d-4a13-b432-b4bd81ae855f
# ╟─ef2798f3-0255-4039-aaa0-0a09b0f25f76
# ╟─6f12b4e3-4c3d-40d6-adec-7f78d7f1b04f
# ╟─01588361-4b3d-4ebe-91af-4d86de9c06b2
# ╠═80810067-023b-48e9-a4cf-fe458fd90189
# ╟─290f6613-3b20-4658-b6ff-28ee19a959c4
# ╠═06b62b7f-01a5-45c7-9abe-b86ac44f51ed
# ╟─892365bc-0464-4d5a-be36-6f938379c58f
# ╠═9f5b7c35-7a7d-4a50-9e20-3cae5c55eef1
# ╠═4884e5b2-4bdc-4ae1-96af-336df91aeede
# ╠═582b64d5-f7fa-4477-98fc-36d08e8b9efa
# ╠═4d982c26-a0ef-407a-9f19-8d876502697f
# ╠═28fa22fc-c9ba-4da3-882c-1a386345c962
# ╠═6fcf1cdd-5780-430d-bdd2-6934370d4fcd
# ╟─c00dc9e0-9970-41ab-9aae-b39205df9f66
# ╠═0cf97d60-8d02-48f8-94af-bbb91fee69d6
# ╟─756ae63a-a811-4cd2-b888-5d946a1e9527
# ╠═8d892dc3-81a3-467b-89fe-8640589451d7
# ╟─5613492b-90fe-4f44-87bb-91e2edf9a960
# ╠═1869882b-95fd-417e-9c53-961abc175cc6
# ╟─8edcd362-5eac-4e6d-bdec-85b7d71479f4
# ╠═4b06228f-ebe2-429a-884a-36df923995dc
# ╟─919e710b-8b63-4f4e-8543-b283ce987689
# ╠═d7b4464b-2a88-4671-b3e1-c2c227889fc4
# ╟─ca2cd52d-1e2c-4176-b3b8-8a8d402882f9
# ╟─f98e7839-7318-4363-9e8e-e3750cc3a400
# ╠═c8fb6db5-c387-4e04-935e-36af5982f202
# ╠═0daae961-bfbc-4807-b550-23ae6c9b3695
# ╟─c64aef5a-5153-4ffc-ac11-9470d11cb87e
# ╠═71989b8d-9dd0-488b-ab34-93f83ef0bf52
# ╟─bbe49f09-9e21-49b3-9cc7-f8b2f1ae0871
# ╟─f3685950-6a36-45c7-9b74-ff8004d1f41b
# ╟─036b12b9-897a-46c3-8b5c-32ebf4cd8855
# ╟─ee51dc26-d8fd-4d14-a1f1-b6665a6a3715
# ╟─9ee2dd6a-378b-408f-a7b4-065837a3e5a0
# ╟─a3043126-bdfb-48e4-ba90-32b1b73b25c5
# ╟─4383fa1e-2c6e-4c2c-a9ed-1ac28d69716f
# ╟─a50265e1-5c47-4d53-8992-134a123ecb33
# ╟─23abd4e6-d32e-4937-a11f-592d7bee58de
# ╟─30e9ca21-7007-483b-8b0f-c6296606b653
# ╟─ea9cb1c7-8cab-4f1f-b7d7-fef68a5b0d89
# ╟─316bf71f-03e7-43ce-ab39-2ca1decb6266
# ╠═51241163-7867-40f2-8c90-498aed295a3f
# ╟─aab2d571-e794-4c35-9ff1-c96415107bce
# ╠═4d343352-9263-4a16-8de7-586ea478a3bc
# ╟─9141e633-6022-4c87-a2ac-296927249549
# ╠═6b123991-0773-40f8-a0ac-6ea353ca0a3e
# ╠═88552670-04af-4b2a-858c-11519c5471b2
# ╟─8113b9b7-c5b9-4e9a-a02b-5d87ed4f60df
# ╠═040aeed9-4ed4-43c0-8c24-eecd73862ece
# ╟─7cbc0cee-6a73-40e1-8399-86892110f8df
# ╠═6a3e75d4-255c-4039-8f0f-b166adf05871
# ╠═27f987d7-7c82-47bc-a441-e9fc4469bc56
# ╠═ec201a3f-e323-49d9-8904-62241bb1434e
# ╟─6e8930b1-2357-4a2e-bc72-388b7dd42258
# ╠═0ed3770e-9909-4fcd-a190-ce71d05ef5bd
# ╟─65d97516-3de2-47a8-b536-dc5430f622cc
# ╟─9cc5da63-fdaf-4d94-8445-cd31fd0888c6
# ╠═7a84fc5c-d956-4187-98ba-87ed0e7641e2
# ╠═cf6667a8-e033-4ef6-a7a3-2d788de8b6ae
# ╠═e6246548-4164-4f55-bfd4-cd7d36564eee
# ╠═1d97d78b-859a-4d79-81f8-bbaf2b741447
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
