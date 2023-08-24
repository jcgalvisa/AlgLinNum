### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 79238d73-62a6-4a7c-bb7b-b85266fceb12
begin
	using PlutoUI
	using HypertextLiteral
	using LinearAlgebra
end

# ╔═╡ 332554ec-b061-4d93-be61-d0c400fdea1a
using BenchmarkTools

# ╔═╡ 6d157051-b6fe-47bf-969d-2384f0c00a21
PlutoUI.TableOfContents(title="Ortogonalización", aside=true)

# ╔═╡ 51ba2480-0c4e-11ee-1e0c-594440b0ada3
md""" # Prerrequisitos """

# ╔═╡ e3cdd44b-7e34-4111-85a6-93adf93c1189
md""" 
# Problema de ortogonalización

Considermos el problema de ortogonaliar las columnas de una matriz dada $A$. Los métodos mas usados son
1. Gram-Schmidt clásico, CGS,  ([Gram](https://en.wikipedia.org/wiki/J%C3%B8rgen_Pedersen_Gram) 1883, Popularizado por [Schmidt](https://en.wikipedia.org/wiki/Erhard_Schmidt) 1907)
2. Gram-Schmidt modificado, MG (Desde 1820 por [Laplace](https://en.wikipedia.org/wiki/Pierre-Simon_Laplace))
3. Householder (H) (usando reflexiones). Desarrollado en 1958 por [Householder](https://en.wikipedia.org/wiki/Alston_Scott_Householder) en el [Oak Ridge National Laboratory](https://www.ornl.gov/)
4. Givens (G) (usando rotaciones). Desarrollado en 1950 por [Givens](https://en.wikipedia.org/wiki/Wallace_Givens) en el [Argone National Laboratory](Argone National Laboratory).
Además se usa 
5. Reortogonalización 

En aplicaciones practicas CGS es poco usado. MGS es mas estable que CGS. H es el mejor método en términos de estabilidad. G tiene mas calculos pero puede ser usado para matrices esparsas (ralas). En ocaciones debido a errores numéricos o perturbaciones es conveniente reortogonalizar, por ejemplo se puede combinar MGS con reortogonalización.


La ortogonalización de vectores tiene muchas aplicaciones. Por ejemplo en la solución estable de sistemas lineales. Se puede por ejemplo usar ortogonalización para calcular de forma estable la factorización $LU$. Existen métodos iterativos para calculo de vectores propios que son basados en ortogonalización. También es usado en el cálculo de la descomposición SVD. 
"""

# ╔═╡ 7cd9ba81-bf09-40b0-bce0-fb5792ca0f7e
md"""
Sea $Q\in \mathbb{R}^{m\times p}$, decimos que $Q$ es ortogonal si $Q^TQ=I_{p\times p}$. Si $Q\in \mathbb{C}^{m\times p}$, decimos que $Q$ es unitaria si $Q^*Q=I$. 

Note que si $q_j$, $j=1,2,\dots,p$, denotan las columas de $Q$, entonces $q_j\in \mathbb{R}^m$ y $q_i^Tq_j=\delta_{ij}$. 
"""

# ╔═╡ 8515b92f-1847-46ae-82cd-ebf326d68b60
md"""
# Factorización A=QR

Dada $A\in\mathbb{R}^{m\times n}$ con $A=[a_1,a_2,\dots,a_n]$ donde $a_i, i=1,2,\dots,n$ denotan las columnasd de $A$, queremos escribir 
$ A= QR$ donde  $Q\in\mathbb{R}^{m\times p}$ y  $R\in\mathbb{R}^{p\times n}$ y con $Q$ una matriz ortogonal. 


Observe que $Q^TQ=I$ pero si $p<n$ entonces $QQ^T\not = I$. Ademas tenemos que $\mbox{Span}(A)=\mbox{Span}\{ a_i\}_{i=1}^n=\mbox{Span}\{q_j\}_{j=1}^p=\mbox{Span}(Q).$


Por otro lado, es común escribir la factorización **completa** $A=QR$ donde 

$A= [ Q_1  \quad Q_2] \left[ \begin{matrix}R_1\\ 0\end{matrix} \right]$ 

y de donde se observa que  
$\mbox{Span}(A)=\mbox{Span}(Q_1)$ y $\mbox{Span}(A^T)^\perp=\mbox{Span}(Q_2).$

Por último observe que si las columnas de $A$ son linealmente independientes, entonces $p=n$. 


"""

# ╔═╡ 08898697-c5f1-49ed-9944-c642ebea6b79
md""" # Gram - Schmidt Clásico """

# ╔═╡ 82176b46-efcc-4493-bfd7-5318daf0e12c
md"""
Considere el caso de columnas linealmente independientes, dados 
$\{a_1,a_2,\dots,a_n\}\subset \mathbb{R}^m$ linealmente independientes se desea generar $\{q_1,q_2,\dots,q_n\}$ tales que $q_i^Tq_j=\delta_{ij}$ y con las siguientes propiedades:


1. $\mbox{Span}\{a_1,a_2,\dots,a_n\} =\mbox{Span}\{{\color{red}{q_1}},a_2,\dots,a_n\}=\dots$ 
$= \mbox{Span}\{{\color{red}{q_1,q_2}},\dots,a_n\}= \dots =  \mbox{Span}\{{\color{red}{q_1,q_2,\dots,q_n}}\}.$
2. $\mbox{Span}\{a_1,a_2,\dots,a_s\}=\mbox{Span}\{\color{red}{q_1,q_2,\dots,q_s}\}$para $s=1,2,\dots,n$.

En el primer paso se tienen dos opciones para $q_1$, $q_1= \frac{a_1}{\|a_1\|}$ o $q_1= -\frac{a_1}{\|a_1\|}$.
Se fija una de ellas. Suponga que ya fueron contruidos $q_1,q_2,\dots,q_{i-1}$ tales que vale 1. y 2. arriba. Para calcular $q_i$ vemos que, si queremos que valga 2., debemos tener, 

$r_{ii}q_i=a_i-\sum_{\ell=1}^{i-1}q_\ell r_{\ell i}$

o 

$a_i= \sum_{\ell=1}^i q_\ell r_{\ell i}.$

Para calcular $r_{j,i}$ con $j=1,2,\dots,i-1,$ multiplicamos por $q_j^T$ para obtener ($q_j^Tq_\ell=\delta_{j,l}$)

$0=q_j^T(r_{ii}q_i)=q_j^T \Big(a_i-\sum_{\ell=1}^{i-1}q_\ell r_{\ell i} \Big) =q_j^Ta_i -r_{ji}$
de donde $r_{ji}= q_j^Ta_i$. Observe tambien que al multiplicar por $q_i^t$ obtenemos 

$r_{ii}=q_i^T(r_{ii}q_i)=q_i^T \Big(a_i-\sum_{\ell=1}^{i-1}q_\ell r_{\ell i} \Big) =q_j^Ta_i.$

Al recopilar todas estas observaciones tenemos el algoritmo de la página 231 del texto guía.

"""

# ╔═╡ d2e396f8-b1a6-4cd2-b6dd-d9a2851c8e75
function QRCGS(A)
    sizeA=size(A)
    Q = zeros(sizeA) #(m,n)
    R = zeros(sizeA[2],sizeA[2]) #(n,n)
    for i = 1:sizeA[2]
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

# ╔═╡ 2e6c2a9e-2ee3-44f5-b77f-0403866f815a
md""" 
Note que en la segunda línea del codigo copiamos la matriz $B$ a la matriz $A$. Esto no es necesario y se hace para una comparación mas justa con el algoritmo de MGS ya que el parametro de la función en Julia pasa por referencia (pass-by-sharing). Podemos verificar el algoritmo con los siguientes ejemplos.
"""

# ╔═╡ 515ea914-bf3d-4300-a894-1ad08d539e87
md""" ## Ejemplo 1 """

# ╔═╡ 182be4ae-a5c7-4530-aff2-e35076523a14
n₁,m₁= 3, 4

# ╔═╡ 7bdd743a-b20c-46df-a979-66670dea11ee
A₁ = rand(m₁,n₁)

# ╔═╡ f29f9e4f-3f21-4754-8a7d-06fc29d93e59
Q₁, R₁ = QRCGS(A₁); display(R₁)

# ╔═╡ 41ed3062-a8a0-4782-9460-48fc17e91f98
display(Q₁)

# ╔═╡ e315e430-2983-4370-834f-2e68e5779f72
opnorm(A₁-Q₁*R₁)

# ╔═╡ 47376a21-a951-4563-9adb-f6c175c50c01
opnorm(Q₁'Q₁-UniformScaling(1))

# ╔═╡ 1e41d138-d79e-4b20-804e-280cbff7e49a
md""" En el siguiente ejemplo consideramos una matriz $A=[A_1,A_2]$ donde $A_2$ es una perturbación de $A_1$. """

# ╔═╡ 5294223c-9c11-4c12-ba79-a493e593b4b1
begin
	n₂ = 100
	m₂ = 1000 
	A₂ = rand(m₂,Int(n₂/2));
	ϵ₂= 1E-5 # tamaño de la perturbación
	A₃=A₂+ϵ₂*rand(m₂,Int(n₂/2));
	A₄=[A₂ A₃];
end

# ╔═╡ 169f8c57-dd3f-420f-ab8f-845827957bda
Q₃, R₃ = QRCGS(A₃);

# ╔═╡ 2979e236-17ff-4adf-b890-faba2dc3af9c
opnorm(A₃-Q₃*R₃)

# ╔═╡ 092220fa-4a6e-493e-b5cd-0d221645af48
opnorm(Q₃'Q₃-UniformScaling(1))

# ╔═╡ 326085e4-b4f3-4d36-ac87-b6ae2e57763b
md"""Aquí observamos la poca estabilidad numérica del CGS ya que el residuo de la ortonormalidad no es pequeño. """

# ╔═╡ c31f342a-0100-46d6-ae93-adfed9248d96
md"""# Gram - Schmidt modificado

Alternativamente al algoritmo anterior podemos proceder como sigue. Después de calcular $q_1=r_{11}a_1$, restamos la componente en $q_1$ de $a_2,a_3,\dots,n$, esto es hacemos 
$a_i=a_i-(q_1^Ta_i)q_1$. Obtenemos *nuevos* vectores $a_2,a_3,\dots,a_n$ ortogonales a $q_1$. Seguidamente calculamos $q_2=r_{22}a_2$ y repetimos el proceso anterior ahora con $a_3,\dots,a_n$. 

Considere el algoritmo de MGS. Algoritmos 5.2.5 de la página 231-232 del texto guía. 

"""

# ╔═╡ 9f9bf293-f4b3-49f2-8b0e-146e6ba52313
function QRMGS(A)
    sizeA=size(A)
    Q = zeros(sizeA) #(m,n)
    R = zeros(sizeA[2],sizeA[2]) #(n,n)
    for i = 1:sizeA[2]
        R[i,i] = norm(A[:,i])
        Q[:,i] = A[:,i]/R[i,i]
        for j = i + 1: sizeA[2]
            R[i,j] = Q[:,i]'A[:,j]
            A[:,j] = A[:,j] - Q[:,i]R[i,j]
        end
    end
    return Q, R
end

# ╔═╡ e363cbe7-a469-4a8e-9533-d13921532697
begin
	n₅ = 30
	m₅ = 300 
	A₅ = rand(m₅,n₅)
	B₅=copy(A₅)
end

# ╔═╡ 224ac346-39f7-45dc-9b50-ec1aeb62038c
Q₅,R₅=QRMGS(B₅);

# ╔═╡ 303868fc-de8e-48a6-95a9-cd81323656f4
cond(A₅)

# ╔═╡ 0c674427-7368-4475-9aac-7db03b247cf8
opnorm(Q₅'*Q₅-UniformScaling(1))

# ╔═╡ 615c8ab5-39fe-4d28-bb8a-63e50a7df795
opnorm(A₅-Q₅*R₅)/opnorm(A₅)

# ╔═╡ cb0111b0-c27e-466b-97f3-bca58a75ca25
begin
n₆= 30
m₆ = 300 
A₆ = rand(m₆,n₆)
B₆=A₆  
end

# ╔═╡ 382da2b9-b903-4029-ad07-97d936529853
@benchmark Q,R=QRCGS(A₆)

# ╔═╡ 9178ab89-2210-4b6d-b1b3-7a33a483acee
@benchmark Q,R=QRMGS(B₆)

# ╔═╡ 5bed8ab7-6e4b-49b6-b7e8-b577f71975d4
begin
n₇ = 200
m₇ = 1000 
A₇ = rand(m₇,Int(n₇/2));
A₈=A₇+ϵ₂*rand(m₇,Int(n₇/2));
#opnorm(A1-A2)
A₉=[A₇ A₈]
B₉=A₉
end

# ╔═╡ bec89150-5967-429d-8539-d08848f65df1
Q₉c,R₉c=QRCGS(A₉);

# ╔═╡ 7e6802da-dea6-4d9d-8db0-ef9abe09fc3b
opnorm(Q₉c'*Q₉c-UniformScaling(1))

# ╔═╡ ad971243-e6d3-4f99-ab56-9071ab77d0d4
cond(A₉)

# ╔═╡ b8f16643-195c-45d9-928e-7a6112485a04
opnorm(A₉-Q₉c*R₉c)/opnorm(A₉)

# ╔═╡ 9518b7b6-45ef-44fb-8f0b-7a315d6a4dfe
Q₉m,R₉m=QRMGS(B₉);

# ╔═╡ f1833e8d-3d98-4748-8373-de059294e98e
opnorm(Q₉m'*Q₉m-UniformScaling(1))

# ╔═╡ 9b2d9220-5c36-47a1-a5a2-a19de55159e0
opnorm(A₉-Q₉m*R₉m)/opnorm(A₉)

# ╔═╡ b0214333-9461-4535-a5a6-f2aeac628e06
md"""Observe que con el MGS los vectores estan más cerca de ser realmente orgotonales que con el CGS. Además, con el CGS el residuo de la factorización es un poco menor que con MGS. Finalmente hacemos la siguiente observación: cuando el residuo de ortogonalidad es grade podemos reortogonalizar. Suponga que inicamos con $A\approx QR$. Entonces podmeos calcular $Q\approx \tilde{Q}\tilde{R}$. Obtenemos $A\approx (\tilde{Q}\tilde{R})R=\tilde{Q}\hat{R}$. """

# ╔═╡ cca162bf-84fd-4032-a90f-522796baacf9
tildeQ₉c,tildeR₉c=QRCGS(Q₉c);

# ╔═╡ c4421b71-29cf-4b43-a71c-5b6070bd4afe
opnorm(tildeQ₉c'tildeQ₉c-UniformScaling(1))

# ╔═╡ 7132f19f-fc36-4984-8d7e-2ff58ae08a02
opnorm(A₉ - tildeQ₉c*(tildeR₉c*R₉c))

# ╔═╡ 183aa2a6-2be3-4516-a0b1-4d52297cbab6
md""" 
# Teorema de ortogonalización

Terminamos enunciado el siguiente resultado.

**Teorema:** Si las columnas de $A$ son linealmente independientes y $r_{ii}>0$ para $i=1,\dots,n$ entonces la factorización $A=Q_1R_1$ con $Q_1$ ortogonal y $R_1$ triangular superior es única. 

**Demostración:** Dado que $A=Q_1R_1$ implica $A^TA=(Q_1R_1)^T(Q_1R_1)= R_1^T ( Q_1^TQ_1)R_1=R_1^TR_1$, el resultado se sigue de la unicidad de la factorización de Cholesky para matrices positivas definidas, en este caso aplicada a $A^TA$. """

# ╔═╡ c98bc2ab-575b-447f-a7cd-f5ba037a0a43
md""" # Reflexiones de Householder

Dado $v\in\mathbb{R}^m$ defina la matrix de Housholder por $\displaystyle H_v=I - 2\frac{vv^T}{v^Tv}.$ Decimos que $v$ es el vector de Hoseholder.  Podemos verificar que 

1. La matriz $H$ es una isometría en la norma $\|\cdot\|_2$.
2. La matriz $H$ es simétrica 
3. La matriz $H$ es ortogonal

Dado $x\in \mathbb{R}^m$, la reflexión de $x$ con eje de reflexión el huperplano $v^\perp$ es dada por 
$y=H_vx = x- 2\frac{v^Tx}{v^Tv}v.
$
Note que para cacular $y$ no es necesario calcular todas la entradas de la matriz $H$, es suficiente calcular el producto interno $s=v^Tx$ y realizar la resta $x-\beta s v$ donde $\beta= 2\frac{1}{v^Tv}$. El valor de $\beta$ puede ser precalculado. 
"""

# ╔═╡ ce3150a6-10ac-40d9-8922-e2bc6de3ead4
md""" 
Dados $x,y\in \mathbb{R}^m$. ¿Existe $v$ tal que $H_vx=y$?. Si $\|x\|_2=\|y\|_2$ la respuesta es afirmativa  y 
$v=\alpha (x-y)$. 

Para usar Householder en ortogonalización procedemos como sigue. Suponga que queremos orgonaliar las columnas de $A=(a_1,a_2,\dots,a_n)\in \mathbb{R}^{m\times n}$. Considere inicialmente el caso de columnas linealmente independientes. Iniciamos la ortogonalización de Householder calculando $v$ que transforme $a_1$ en un multiplo del primer vector canónico $e_1=(1,0,\dots,0)$. Es decir $H_v a_1= \sigma e_1$ con $\sigma = \pm \|x\|_2$.



Suponga que $x\in\mathbb{R}^m$ y queremos calcular $v$ tal que $H_v x$ sea múltiplo de $e_1$, digamos $\sigma e_1$. Debe ser $\sigma= \pm\|x\|_2$. Observe que 

$x-\sigma e_1 = (x_1-\sigma,x_2,x_3,\dots,x_m).$

En la resta $x_1-\sigma$ puede ocurrir cancelacion catastrófica. Tenemos las siguientes opciones para evitar la cancelación catastrófica, 

*Opción 1:* Observe que si $x_1>0$ entonces $x$ esta mas cerca a  $\|x\|_2e_1$ que a $-\|x\|_2e_1$. Entonces podemos tomar $\sigma = -\mbox{sign}(x_1)$ y evitamos la cancelación catastrófica. 
Despues de calcular $v_1=x_1-\sigma$ podemos calcular 

$\beta =  \frac{2}{v^Tv}= \frac{2}{ v^Tv} = \frac{2}{ v_1^2+s}$

donde $s=x_2^2+\dots+x_m^2$.

*Opción 2:* Podemos seleccionar $\sigma=\|x\|_2$ pero tener cuidado al calcular $x_1-\sigma$. 
- Si $x_1<0$ entonces podemos calcular $v_1=x_1-\sigma$ (no hay resta).
- Si $x_1>0$ entonces calculamos 
    $x_1-\sigma = (x_1-\sigma) \frac{x_1+\sigma}{x_1+\sigma}= \frac{x_1^2-\|x\|_2^2}{x_1+\sigma}= 
    -\frac{x_2^2+\dots+x_m^2}{x_1+\sigma} = - \frac{s}{x_1+\sigma}$ donde $s=x_2^2+\dots+x_m^2$.
    
    Note tambien que se puede calcular $v^Tv$ como sigue (recuerde que $\sigma=\|x\|_2$)
    
$v^Tv= (x_1-\sigma)^2+x_2^2+\dots+x_m^2 = x_1^2-2\sigma x_1+\sigma^2+x_2^2+\dots+x_m^2= 
    2\sigma^2 - 2\sigma x_1= 2\sigma (\sigma -x)= -2\sigma v_1.$

En muchas implementaciones para ahorar en memoria se normaliza $v$ de tal forma que $\tilde{v}=v/v_1$ y se guardan las entradas $n-1$ últimas entradas de  $v_2/v_1,\dots,v_m/v_1$  en las entradas de $x_2,x_3, \dots,x_m$ (debajo de la diagonal de la matriz $R$ en aplicacion de ortonormalización).

Tenemos entonces el siguiente algorimos que dado un vector $x$, calcula de forma correcta el vector de Householder $v$. Algoritmo 5.1.1, página 210 del texto guía. """

# ╔═╡ 710136de-30ce-48aa-9911-58eb4f48646a
function house(x)
    n = size(x,1) #(m,n)
    s = x[2:n]'x[2:n]
    v = [1; x[2:n]]
    
    if s == 0
        β = 0
        
    else
        mu = sqrt(x[1]^2 + s)
        if x[1] <= 0
            v[1] = x[1] - mu
        else
            v[1] = -s /(x[1]+mu)
        end
        β = 2((v[1])^2)/(s + (v[1])^2)
        v=v/v[1]
    end
    return v, β
end

# ╔═╡ 79c5957b-1ff2-4a0c-bf48-cee6540a72b5
md""" Considere el siguiente ejemplo."""

# ╔═╡ 15d46ffa-7175-4b72-8ed2-054abb09e1b0
begin
x=[0.5;sqrt(3)/2;4];
v, β= house(x)
println(v)
println(β)
print(x)
end

# ╔═╡ 6c6eea6c-2754-4b3d-88ad-237a10329227
H=UniformScaling(1)- 2(v*v')/(v'v)

# ╔═╡ 2df12392-0700-4301-afdc-c47729fa9247
H*x

# ╔═╡ eb300638-4986-4c88-be9f-2474f59bd4e3
md""" 
Note que si estamos haciendo la factorización de una sola columna, $x$, entonces $H^T$ corresponde a la matriz $Q$ de la factorización $QR$, **pero esta vez la matriz $Q$ es cuadrada y genera todo $\mathbb{R}^m$**. Es decir, obtenemos la factorización completa. En esta caso las filas $2:m$ de la matriz $H$ son vectores ortogonales a $x$. Recuerde que si queremos obtener la factorizació completa con GS tendríamos que iniciar con $m$ columnas, es decir, tendriamos que encontrar $m-n$ columnas linealmente independientes a las que tenemos y despues aplicar GS. 
"""

# ╔═╡ 0cbb28b5-2980-413f-be18-d311ee092bb8
md""" 
Consideremos nuevamente la matrix $A=(a_1,a_2,\dots,a_n)$. Aplicando el procedimiento anterior construimos $H_1$ tal que 
$A_1=H_1A = (\sigma e_1, H_1a_2,\dots, H_1a_n)$
es una matriz con $A_1(2:m,1)=(0,\dots,0)$. 
Ahora debemos obtener entradas nulas también en la segunda columna debajo de la diagonal.
Para continuar con la triagularización de Householder construimos una reflexion de Householder en $\mathbb{R}^{m-1}$ que transforme el vector $A_1(2:m,{\color{red}{2}})$ en el vector $(1,0,\dots,0)\in \mathbb{R}^{m-1}$. Sea $\tilde{H}_2$ esta transformación. Defina 

$H_2= \begin{pmatrix}1 & 0\\ 0 &\tilde{H}_2\end{pmatrix}.$

Vemos que $H_2H_1a_2= (A_1(1,2), \sigma_2,0,\dots,0)\in \mathbb{R}^m$. 
Concluimos que $A_2= H_2H_1 A =  ( \sigma e_1, H_2H_1a_2, H_2H_1a_3,\dots, H_2H_1a_n)$
es una matriz con entradas nulas debajo de la diagonal en las dos primeras columnas. 

Podemos continuar este proceso hasta la ultima columna y obtenemos una matriz triangular superior. Es decir, obtenemos 
$H_nH_{n-1}\cdots  H_1  A = R
$
y por tanto $ A= H_1^TH_2^T\cdots H_n^T R =QR$ donde $Q=H_1H_2\cdots H_n$. 

Note que quedamos con la matriz $R$ pero para poder obtener la matriz $Q$ debemos, en prinicipio, multiplicar las matrices $H_i$, lo cual es computacionalmente costoso. Veamos primero como queda el algoritmos para obtener $R$ y después  mostraremos un algoritmo para construir $Q$. 


Considere el algorimos en la página 211 que usa el algoritmo anterior para obtener la matriz $R$ de la factorizacion $QR$. Recuerde que asumimos que las columnas de $A$ son linealmente independientes. 
"""

# ╔═╡ 2fec112c-3819-432d-a17d-abeaef833c29
function Rhouse(A)
    m,n = size(A)
    for j = 1:n
        v, β = house(A[j:m,j])
        A[j:m,j:n] = ( UniformScaling(1) - β*v*v')A[j:m,j:n]
#        A[j+1:m,j]=v[2:end]
    end
    return A
end

# ╔═╡ e5aa386a-a035-4c9a-8cdd-d43b33162347
begin
n₁₀ = 3 ;
m₁₀ = 5 ;
A₁₀ = rand(m₁₀,n₁₀)
B₁₀ = A₁₀
end

# ╔═╡ 4ac1b3e3-4a6f-43ef-abec-fc0313b28045
R₁₀h=Rhouse(B₁₀)

# ╔═╡ 3a58e7a0-f271-4379-9315-a22eebd46059
md""" 
##  Matriz $Q$ acumulación progresiva
"""

# ╔═╡ d81423c0-64ea-409a-84bb-c20bcf4046c6
function QFA(A)
    Q = UniformScaling(1)
    m,n = size(A)
    for j = 1:n
        v,β = house(A[j:m,j])
        v1 = zeros(j-1,1)
        v = [v1;v]
        Qj = UniformScaling(1)-β*v*v'
        Q = Q*Qj
        A=Qj*A;
    end
    return Q, A
end

# ╔═╡ 5baa3cfe-3323-4c48-9f8c-baa138f4f1e4
begin
Q₁₀,R₁₀=QFA(copy(A₁₀))
#Q11'*Q11
A₁₀-Q₁₀*R₁₀
end

# ╔═╡ c9235c7a-4d03-4e52-8a5b-94cd6c21b7af
Q₁₀[:,1:n₁₀]

# ╔═╡ 76633fe0-184c-4dc6-a715-b8e8510fcbeb
Q₁₀[:,n₁₀+1:m₁₀]

# ╔═╡ 87ca587a-7283-4991-9e12-1b87e8ba575d
begin
	n₁₁ = 30 ;
	m₁₁ = 50 ;
	A₁₁ = rand(m₁₁,n₁₁);
	B₁₁ = A₁₁;
	Q₁₁,R₁₁=QFA(A₁₁);
end

# ╔═╡ 43d5dc23-ab26-4752-8210-ac393f1a0268
opnorm(UniformScaling(1)-Q₁₁'Q₁₁)

# ╔═╡ 2c51de0a-b4d3-485c-897c-95da6fcceafc
opnorm(A₁₁-Q₁₁*R₁₁)

# ╔═╡ 8bd19164-dcc4-4a16-b552-1d11030c0366
md""" 
##  Matriz $Q$: acumulación progresiva
"""

# ╔═╡ dbdc714a-231b-405c-be14-8c60d50acc7f
md""" 

Mencionamos antes que realizar la multiplicación de las matrices de Householder puede ser costoso. En su lugar, cuando se requiere la $R$ podemos usar al siguiente resulado. 

**Teorema (Representación por bloques de Householder)** Si definimos 
$Q_1=H_1, \quad Q_2=H_1H_2, \dots, Q_i=H_1H_{2}\cdots H_{i-1}H_i= Q_{i-1}H_i$
entonces $Q_i=I +W_i Y_i^T$, $i=1,2,\dots,n$, $W_i, Y_i \in \mathbb{R}^{m\times i}$.
"""

# ╔═╡ bc00e335-4fac-4c07-86d9-7a2e7a05a61d
function VβHouse(A)
    m,n = size(A)
    β   = zeros(n,1)
    Y   = zeros(m,n)

    Am  = copy(A)
    
    # Primera actualización
    v1,β[1] = house(Am[:,1])
    Y[:,1]  = v1
    W       = -β[1]*v1
    Am      = ( UniformScaling(1) - β[1]*v1*v1')Am
    
    #De ahi para adelante
    for j = 2:n
        v1,β[j] = house(Am[j:m,j])
        Y[:,j]  = [ zeros(j-1,1) ; v1 ] 
        z       = -β[j]*(UniformScaling(1)+W*Y[:,1:j-1]')*Y[:,j]
        W       = [ W z ]
        Am[j:m,j:n] = ( UniformScaling(1) - β[j]*v1*v1')Am[j:m,j:n]
    end
    
    
    return Y,W,Am
end

# ╔═╡ b37cf5d3-e7b2-41c3-9fcb-e232718f7957
Y,W,R = VβHouse(A₁₀)

# ╔═╡ 079e1d22-d81c-4d2d-91db-d929f8898bff
display(Y)

# ╔═╡ 4d628c4b-4311-40d2-8db8-4cab9e768bac
display(W)

# ╔═╡ 642a2da3-de23-422f-ace4-48667731cbfa
Q₁₀r = UniformScaling(1) + W*Y'

# ╔═╡ 2a4c2eec-c3f1-4f8e-980c-0ff607e1e46d
A₁₀-Q₁₀r*R

# ╔═╡ d1b87663-6755-4eaa-9123-8717f95095a4
md"""# Rotaciones de Givens """

# ╔═╡ 53e98911-65ae-48f0-96d4-bfb2bc04585f
md""" También conocido como rotaciones de planos. Se aplica en el caso de matrices dispersas (sparse matrices) ya que tiene un efecto "local". """

# ╔═╡ 22106606-9bee-40cc-8d2f-17f36dda6937
md""" Dados $i,k\in \{1,\dots,m\}, \theta\in\mathbb{R}$, defina:
$c=\cos(\theta), s=\sin(\theta)$,

$G(i,k,\theta)=\begin{pmatrix} 
\vdots & \vdots  & \vdots & \vdots & \vdots \\ 
\dots & c  & \dots & s & \dots \\ 
\vdots & \vdots  & \vdots & \vdots & \vdots \\ 
\dots & -s  & \dots & c& \dots \\
\vdots & \vdots  & \vdots & \vdots & \vdots \\ 
\end{pmatrix}$

Aplicar $G^T(i,k,\theta)$ a una matriz $A=(a_1,\dots,a_n)$ equivale a "rotar" los ejes $i$ y $k$.
"""

# ╔═╡ 38539181-36e0-444a-b93e-ab83600e983b
md""" 
Si $y=G^T(i,k,\theta)x$ con $x\in \mathbb{R}^n$ tenemos que 

$y_i=cx_i -sx_k, \quad y_k=sx_i+cx_k$

y además $y_j=x_j$ para $j\not=i$,$j\not=k$. Luego

$G^T(i,k,\theta)A=G^T(i,k,\theta)
\begin{pmatrix} 
r_1 \\ 
\vdots\\
r_i \\ 
\vdots \\
r_k\\
\vdots\\
r_m
\end{pmatrix}=
\begin{pmatrix} 
r_1 \\ 
\vdots\\
cr_i-sr_k \\ 
\vdots \\
sr_i+cr_k \\ 
\vdots \\
r_m
\end{pmatrix}$
"""

# ╔═╡ 6b1c0220-f772-4cea-99b6-6d5b221fc020
md"""Si

$x_i=||x||\cos(\phi), \quad x_k=||x||\sin(\phi)$
entonces

$y_i=||x||\Big( \cos(\theta)\cos(\phi)-\sin(\theta)\sin(\phi) \Big) =||x||
\cos(\theta+\phi)$
y

$y_i=||x||\sin(\theta+\phi).$
Si queremos $y_k=0$ tenemos así dos alternativas, $\theta=-\phi$ o $\theta=\pi-\phi$.
Escogemos $\theta=-\phi$. De donde

$c=\cos(\theta)=\frac{x_i}{\sqrt{x_i^2+x_k^2}}, \quad
s=\sin(\theta)=-\frac{x_k}{\sqrt{x_i^2+x_k^2}}$
"""

# ╔═╡ 99113f83-5cbf-4d08-b164-c6f4c7520c35
md"""Tenemos el siguiente algoritmo (que evita división por la menor componente y también evita división por norma pequeña."""

# ╔═╡ d9b5a5f5-0461-4455-93ea-1915decb3145
function Givens(a,b)
    if b!=0
        if abs(b)>abs(a)
			# tipo 1
            τ=-a/b
            s=-1/sqrt(1+τ^2)
            c=s*τ
        else
            # tipo 2
			τ=-b/a
            c=1/sqrt(1+τ^2)
            s=c*τ
        end
    end
    return c,s
end        

# ╔═╡ 2f1d478d-0b6e-4443-a409-3f00fb7eee7f
md"""Observe que $|\tau|\leq 1$ e que para calsuclar $s$ y $c$ necesitamos de 5 operaciones y una raíz cuadrada.  """

# ╔═╡ 51af97af-dc36-4a5b-81ae-0eaaad5a57de
begin
	v₁₁=[8; 6]
	c,s = Givens(v₁₁[1],v₁₁[2])
	G=[c -s; s c]
end

# ╔═╡ e3c682ea-ea20-4a89-8e29-71aba1b13a7d
G*v₁₁

# ╔═╡ 1bd6aef4-8a80-4245-8b02-cd2713341f6b
md"""Tenemos la siguiente función que  aplica una 
rotación de Gives dadas dos filas de una matriz $A$."""

# ╔═╡ e57d51a4-0005-4b14-9450-73fa43cccceb
function Giv(A, i, k, c, s)
    for j=1:size(A)[2]
        τ1 = A[i,j]
        τ2 = A[k,j]
    
        A[i,j] = c*τ1 - s*τ2
        A[k,j] = s*τ1 + c*τ2
    end
    return A
end

# ╔═╡ 6f0a344c-89a8-41e0-a3d6-7ee2013d678f
begin
	A₁₂ = floor.(10*rand( 3, 3))
	B₁₂=copy(A₁₂)
end

# ╔═╡ 4d3d47ee-8362-41ea-bbbb-8f067db1aebe
begin
	if( A₁₂[3,1] !=0.0)
	c₁, s₁ = Givens(A₁₂[2,1], A₁₂[3,1])
	A₁₃= Giv(B₁₂ , 2,3 ,c₁,s₁)
	end
	A₁₃
end

# ╔═╡ 2738ef60-adb7-4068-b24a-c82b0fc91abf
begin
	c₂, s₂ = Givens(A₁₃[1,1], A₁₃[2,1])
	D = Giv(A₁₃, 1,2,c₂,s₂)
end

# ╔═╡ 027fbf5b-1419-477c-8048-7349264927f9
md""" 
## Givens rápido
"""

# ╔═╡ 025a34af-54bc-43b5-9f82-cbdee6003655
md""" 
Una matriz $G(i,k,\theta)$ puede escribirse (según sea tipo 1 o 2) como: 


$G(i,k,\theta)^T=
\mbox{diag}(1,\dots,1,-s,1,\dots,1,s,1,\dots)
\begin{pmatrix} 
\vdots & \vdots  & \vdots & \vdots & \vdots \\ 
\dots & -\tau  & \dots & 1 & \dots \\ 
\vdots & \vdots  & \vdots & \vdots & \vdots \\ 
\dots & 1  & \dots & \tau& \dots \\
\vdots & \vdots  & \vdots & \vdots & \vdots \\ 
\end{pmatrix}=S_1(i,k,s)T_1(i,k,\tau)$
o

$G(i,k,\theta)^T=
\mbox{diag}(1,\dots,1,c,1,\dots,1,c,1,\dots)
\begin{pmatrix} 
\vdots & \vdots  & \vdots & \vdots & \vdots \\ 
\dots & 1  & \dots & -\tau & \dots \\ 
\vdots & \vdots  & \vdots & \vdots & \vdots \\ 
\dots & \tau  & \dots & 1& \dots \\
\vdots & \vdots  & \vdots & \vdots & \vdots \\ 
\end{pmatrix}=S_2(i,k,c)T_2(i,k,\tau).$
"""


# ╔═╡ dce92edf-e838-434b-bef4-9851e5af147c
md"""Tenemos así que, por ejemplo  

$G^T(i_1,k_1,\theta_1)A=S_p(i_1,k_1,c)T_p(i_1,k_1,\tau_1)=S_p(i_1,k_1,c_1)\widetilde{A}$

donde $A(1,k)=0$. Para el siguiente paso, tendriamos, 

$G^T(i_2,k_2,\theta_2)G^T(i_1,k_1,\theta_1)A= S_q(i_2,k_2,c_2)T_q(i_2,k_2,\tau_2)C_p(i_1,k_1,c_1)\widetilde{A}$
"""

# ╔═╡ 04f1ebd6-13c7-4310-99d7-1abbc8963001
md""" 

Queremos proceder de tal manera que la multiplicación $C_p(i_1,k_1,c_1)\widetilde{A}$ no se haga y solo acumule hasta el final. Para esto consideramos dos casos dependiendo si la rotación de Givens es tipo 1 o tipo 2. Para simplificar consideremos sol las entradas implicadas. Por ejemplo en el caso tipo 1, tenemos, 
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
BenchmarkTools = "~1.3.2"
HypertextLiteral = "~0.9.4"
PlutoUI = "~0.7.51"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.2"
manifest_format = "2.0"
project_hash = "d0069486257542c58ff4f12284b8908f62265555"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "d9a9701b899b30332bbcb3e1679c41cce81fb0e8"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.2"

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
git-tree-sha1 = "5a6ab2f64388fd1175effdf73fe5933ef1e0bac0"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "b478a748be27bd2f2c73a7690da219d0844db305"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.51"

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

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

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
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

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
# ╠═6d157051-b6fe-47bf-969d-2384f0c00a21
# ╠═51ba2480-0c4e-11ee-1e0c-594440b0ada3
# ╠═79238d73-62a6-4a7c-bb7b-b85266fceb12
# ╟─e3cdd44b-7e34-4111-85a6-93adf93c1189
# ╟─7cd9ba81-bf09-40b0-bce0-fb5792ca0f7e
# ╟─8515b92f-1847-46ae-82cd-ebf326d68b60
# ╟─08898697-c5f1-49ed-9944-c642ebea6b79
# ╟─82176b46-efcc-4493-bfd7-5318daf0e12c
# ╠═d2e396f8-b1a6-4cd2-b6dd-d9a2851c8e75
# ╟─2e6c2a9e-2ee3-44f5-b77f-0403866f815a
# ╠═515ea914-bf3d-4300-a894-1ad08d539e87
# ╠═182be4ae-a5c7-4530-aff2-e35076523a14
# ╠═7bdd743a-b20c-46df-a979-66670dea11ee
# ╠═f29f9e4f-3f21-4754-8a7d-06fc29d93e59
# ╠═41ed3062-a8a0-4782-9460-48fc17e91f98
# ╠═e315e430-2983-4370-834f-2e68e5779f72
# ╠═47376a21-a951-4563-9adb-f6c175c50c01
# ╟─1e41d138-d79e-4b20-804e-280cbff7e49a
# ╠═5294223c-9c11-4c12-ba79-a493e593b4b1
# ╠═169f8c57-dd3f-420f-ab8f-845827957bda
# ╠═2979e236-17ff-4adf-b890-faba2dc3af9c
# ╠═092220fa-4a6e-493e-b5cd-0d221645af48
# ╟─326085e4-b4f3-4d36-ac87-b6ae2e57763b
# ╟─c31f342a-0100-46d6-ae93-adfed9248d96
# ╠═9f9bf293-f4b3-49f2-8b0e-146e6ba52313
# ╠═e363cbe7-a469-4a8e-9533-d13921532697
# ╠═224ac346-39f7-45dc-9b50-ec1aeb62038c
# ╠═303868fc-de8e-48a6-95a9-cd81323656f4
# ╠═0c674427-7368-4475-9aac-7db03b247cf8
# ╠═615c8ab5-39fe-4d28-bb8a-63e50a7df795
# ╠═332554ec-b061-4d93-be61-d0c400fdea1a
# ╠═cb0111b0-c27e-466b-97f3-bca58a75ca25
# ╠═382da2b9-b903-4029-ad07-97d936529853
# ╠═9178ab89-2210-4b6d-b1b3-7a33a483acee
# ╠═5bed8ab7-6e4b-49b6-b7e8-b577f71975d4
# ╠═bec89150-5967-429d-8539-d08848f65df1
# ╠═7e6802da-dea6-4d9d-8db0-ef9abe09fc3b
# ╠═ad971243-e6d3-4f99-ab56-9071ab77d0d4
# ╠═b8f16643-195c-45d9-928e-7a6112485a04
# ╠═9518b7b6-45ef-44fb-8f0b-7a315d6a4dfe
# ╠═f1833e8d-3d98-4748-8373-de059294e98e
# ╠═9b2d9220-5c36-47a1-a5a2-a19de55159e0
# ╟─b0214333-9461-4535-a5a6-f2aeac628e06
# ╠═cca162bf-84fd-4032-a90f-522796baacf9
# ╠═c4421b71-29cf-4b43-a71c-5b6070bd4afe
# ╠═7132f19f-fc36-4984-8d7e-2ff58ae08a02
# ╟─183aa2a6-2be3-4516-a0b1-4d52297cbab6
# ╟─c98bc2ab-575b-447f-a7cd-f5ba037a0a43
# ╟─ce3150a6-10ac-40d9-8922-e2bc6de3ead4
# ╠═710136de-30ce-48aa-9911-58eb4f48646a
# ╟─79c5957b-1ff2-4a0c-bf48-cee6540a72b5
# ╠═15d46ffa-7175-4b72-8ed2-054abb09e1b0
# ╠═6c6eea6c-2754-4b3d-88ad-237a10329227
# ╠═2df12392-0700-4301-afdc-c47729fa9247
# ╟─eb300638-4986-4c88-be9f-2474f59bd4e3
# ╟─0cbb28b5-2980-413f-be18-d311ee092bb8
# ╠═2fec112c-3819-432d-a17d-abeaef833c29
# ╠═e5aa386a-a035-4c9a-8cdd-d43b33162347
# ╠═4ac1b3e3-4a6f-43ef-abec-fc0313b28045
# ╟─3a58e7a0-f271-4379-9315-a22eebd46059
# ╠═d81423c0-64ea-409a-84bb-c20bcf4046c6
# ╠═5baa3cfe-3323-4c48-9f8c-baa138f4f1e4
# ╠═c9235c7a-4d03-4e52-8a5b-94cd6c21b7af
# ╠═76633fe0-184c-4dc6-a715-b8e8510fcbeb
# ╠═87ca587a-7283-4991-9e12-1b87e8ba575d
# ╠═43d5dc23-ab26-4752-8210-ac393f1a0268
# ╠═2c51de0a-b4d3-485c-897c-95da6fcceafc
# ╟─8bd19164-dcc4-4a16-b552-1d11030c0366
# ╟─dbdc714a-231b-405c-be14-8c60d50acc7f
# ╠═bc00e335-4fac-4c07-86d9-7a2e7a05a61d
# ╠═b37cf5d3-e7b2-41c3-9fcb-e232718f7957
# ╠═079e1d22-d81c-4d2d-91db-d929f8898bff
# ╠═4d628c4b-4311-40d2-8db8-4cab9e768bac
# ╠═642a2da3-de23-422f-ace4-48667731cbfa
# ╠═2a4c2eec-c3f1-4f8e-980c-0ff607e1e46d
# ╟─d1b87663-6755-4eaa-9123-8717f95095a4
# ╟─53e98911-65ae-48f0-96d4-bfb2bc04585f
# ╠═22106606-9bee-40cc-8d2f-17f36dda6937
# ╠═38539181-36e0-444a-b93e-ab83600e983b
# ╠═6b1c0220-f772-4cea-99b6-6d5b221fc020
# ╠═99113f83-5cbf-4d08-b164-c6f4c7520c35
# ╠═d9b5a5f5-0461-4455-93ea-1915decb3145
# ╠═2f1d478d-0b6e-4443-a409-3f00fb7eee7f
# ╠═51af97af-dc36-4a5b-81ae-0eaaad5a57de
# ╠═e3c682ea-ea20-4a89-8e29-71aba1b13a7d
# ╠═1bd6aef4-8a80-4245-8b02-cd2713341f6b
# ╠═e57d51a4-0005-4b14-9450-73fa43cccceb
# ╠═6f0a344c-89a8-41e0-a3d6-7ee2013d678f
# ╠═4d3d47ee-8362-41ea-bbbb-8f067db1aebe
# ╠═2738ef60-adb7-4068-b24a-c82b0fc91abf
# ╟─027fbf5b-1419-477c-8048-7349264927f9
# ╠═025a34af-54bc-43b5-9f82-cbdee6003655
# ╟─dce92edf-e838-434b-bef4-9851e5af147c
# ╠═04f1ebd6-13c7-4310-99d7-1abbc8963001
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
