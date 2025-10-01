### A Pluto.jl notebook ###
# v0.20.4

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
PlutoUI.TableOfContents(title="Ortogonalización", indent=true, depth=4, aside=true)

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

# ╔═╡ 46d74abc-0ee9-4013-ae45-351c1960163c
md"""
Aplicando el algoritmo anterior, obtenemosd los siguientes resultados, 
"""

# ╔═╡ f29f9e4f-3f21-4754-8a7d-06fc29d93e59
Q₁, R₁ = QRCGS(A₁); display(R₁)

# ╔═╡ 41ed3062-a8a0-4782-9460-48fc17e91f98
display(Q₁)

# ╔═╡ dfd4b7d6-4552-471c-9bc9-647ddacb17ae
md"""Para indagar sobre la calidad de los resultados, podemos calcular el residuo relativo de la factorización y el residio de la ortogonalización: """

# ╔═╡ e315e430-2983-4370-834f-2e68e5779f72
opnorm(A₁-Q₁*R₁)/opnorm(A₁)

# ╔═╡ 47376a21-a951-4563-9adb-f6c175c50c01
opnorm(Q₁'Q₁-UniformScaling(1))

# ╔═╡ a825ec17-d3f4-4167-9f45-ab61812883dc
md""" ## Ejemplo 2 """

# ╔═╡ 1e41d138-d79e-4b20-804e-280cbff7e49a
md""" En el siguiente ejemplo consideramos una matriz $A=[A_1,A_2]$ donde $A_2$ es una perturbación de $A_1$. """

# ╔═╡ 5294223c-9c11-4c12-ba79-a493e593b4b1
begin
	n₂ = 100
	m₂ = 1000
	A₂ = randn(m₂,Int(n₂/2));
	ϵ₂= 1E-10 # tamaño de la perturbación
	A₃=A₂+ϵ₂*rand(m₂,Int(n₂/2));
	A₄=[A₂ A₃];
end

# ╔═╡ 169f8c57-dd3f-420f-ab8f-845827957bda
Q₃, R₃ = QRCGS(A₃);

# ╔═╡ 2979e236-17ff-4adf-b890-faba2dc3af9c
opnorm(A₃-Q₃*R₃)/opnorm(A₃)

# ╔═╡ 092220fa-4a6e-493e-b5cd-0d221645af48
opnorm(Q₃'Q₃-UniformScaling(1))

# ╔═╡ 326085e4-b4f3-4d36-ac87-b6ae2e57763b
md"""Aquí intentamos probar la estabilidad numérica del CGS calculando ell residuo de la ortonormalidad. """

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

# ╔═╡ efc140bc-4fff-4631-93cc-8749fd159f6e
md""" ## Ejemplo 3 
Aplicamos en algoritmo anterior, calculamos los residuos de la factorización y la ortogonalización y comparamos (usando benchmarktools) el tiempo de ejecución del MGS con el CGS.
"""

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

# ╔═╡ 49829808-d34b-402c-b6a9-29f2f55dba03
md""" ## Ejemplo 4
Ahora consideramos el caso de columnas cercanas (casi linelmente independientes). 
"""

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

# ╔═╡ 1ff4c715-6352-465c-8575-a8ccae9ee708
md"""Observe que con el MGS los vectores estan más cerca de ser realmente orgotonales que con el CGS. Además, con el CGS el residuo de la factorización es un poco menor que con MGS."""

# ╔═╡ b0214333-9461-4535-a5a6-f2aeac628e06
md"""
# Re-ortogonalización  
Hacemos la siguiente observación: cuando el residuo de ortogonalidad es grade podemos reortogonalizar. Suponga que inicamos con $A\approx QR$. Entonces podmeos calcular $Q\approx \tilde{Q}\tilde{R}$. Obtenemos $A\approx (\tilde{Q}\tilde{R})R=\tilde{Q}\hat{R}$. """

# ╔═╡ 4df2d1d6-1a1b-4b53-9390-e661acd16dd2
md""" ## Ejemplo 4 (continuación)"""

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
md""" 
## Ejemplo 5
Considere el siguiente ejemplo."""

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

# ╔═╡ a0587521-f2d9-4ba0-a5e3-86a3be90473f
md""" ## Uso de Householder en ortogonalización"""

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

"""

# ╔═╡ 33dd679d-315a-454b-84b0-569a3d380405
md""" 
**Teorema (Representación por bloques de Householder)** Si definimos 
$Q_1=H_1, \quad Q_2=H_1H_2, \dots, Q_i=H_1H_{2}\cdots H_{i-1}H_i= Q_{i-1}H_i$
entonces $Q_i=I +W_i Y_i^T$, $i=1,2,\dots,n$, $W_i, Y_i \in \mathbb{R}^{m\times i}$.

**Demostración:** Tensmos que $Q_1=H_1=I -\beta_1 v_1 v_1^T==I +W_1 Y_1^T$ donde $W_1=-\beta_1v_1$ y $Y_1=v_1$. Continuando por inducción suponga que 
$Q_j=I +W_j Y_j^T$, con $W_j,Y_j\in \mathbb{R}^{m\times j}$. Tenemos que 

$\begin{align}
Q_{j+1}=Q_jH_{j+1}&= (I +W_j Y_j^T)(I -\beta_{j+1}v_{j+1}  v_{j+1}^T)\\
&= I +W_j Y_j^T -\beta_{j+1}Q_jv_{j+1}  v_{j+1}^T\\
&= I +W_j Y_j^T +z_{j+1} v_{j+1}^T \quad \mbox{ con } \quad z_{j+1}=-\beta_{j+1}Q_jv_{j+1},\\
&=  I +[W_j\quad  z_{j+1}][Y_j \quad v_{j+1}]^T\\
&= I +W_{j+1} Y_{j+1}^T.
\end{align}$
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
        z       = -β[j]*(Y[:,j]+W*Y[:,1:j-1]'*Y[:,j])
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

# ╔═╡ fb0292ab-884b-4c56-9151-dbceabd45a56
md""" ## Observaciones

El costo de este algoritmo es del orden de $2n^2m-2n^3/3$.

Una ventaja empirica de Householder es que opera con transformaciones ortogonales, así que preserva normas, mientras que MGS no tiene transformaciones ortogonales. 

Otra ventaja es que cuando los vectores son casi linealmente independientes, en las cuentas del algoritmos de Householder se hace el cálculo numerico con cuidado, esto no se hace en MGS.
"""

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

# ╔═╡ 1047dd07-a6aa-4fc0-a33f-1e6de0d72d25
md""" ## Ejemplo 6"""

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
#	A[i,:]=c*A[i;:]-s*A[k;:]
    return A
end

# ╔═╡ 8e6505d7-dc35-4118-980e-4dbb6ae8ba4b
md""" ## Ejemplo 7"""

# ╔═╡ 6f0a344c-89a8-41e0-a3d6-7ee2013d678f
begin
	A₁₂ = floor.(10*rand( 3, 3))
	B₁₂=copy(A₁₂)
end

# ╔═╡ 4d3d47ee-8362-41ea-bbbb-8f067db1aebe
begin
	if( A₁₂[2,1] !=0.0)
	c₁, s₁ = Givens(A₁₂[1,1], A₁₂[2,1])
	A₁₃= Giv(B₁₂ , 1, 2,c₁,s₁)
	end
	A₁₃
end

# ╔═╡ 2738ef60-adb7-4068-b24a-c82b0fc91abf
#begin
#	c₂, s₂ = Givens(A₁₃[1,1], A₁₃[2,1])
#	D = Giv(A₁₃, 1,2,c₂,s₂)
#end

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

$G^T(i_2,k_2,\theta_2)G^T(i_1,k_1,\theta_1)A= S_q(i_2,k_2,c_2)T_q(i_2,k_2,\tau_2)S_p(i_1,k_1,c_1)\widetilde{A}$
"""

# ╔═╡ 4de38f69-5406-42b4-89f7-cccb4def75d9
md"""
Queremos proceder de tal manera que la multiplicación $C_p(i_1,k_1,c_1)\widetilde{A}$ no se haga y solo acumule hasta el final. Queremos escribir 

$T_q(i_2,k_2,\tau_2)S_p(i_1,k_1,c_1)= \text{ Digonal}\times \text{Matrix similar to S}$
"""

# ╔═╡ 04f1ebd6-13c7-4310-99d7-1abbc8963001
md""" 
Para esto consideramos dos casos dependiendo si la rotación de Givens es tipo 1 o tipo 2. Para simplificar consideremos solo las entradas implicadas en nuestra notación. Por ejemplo en el caso tipo 1, tenemos, 
"""

# ╔═╡ 0fcd6e9f-2e4b-4f66-8e5a-0e4de6f6fee2
md""" 
$
G(a,b,\theta)^T
= \begin{pmatrix}
c&-s\\ s&c\end{pmatrix}=
\begin{pmatrix}
-s&0\\ 0&s\end{pmatrix}\begin{pmatrix}
-\tau&1\\ 1&\tau
\end{pmatrix}$

con $\tau=c/s$ y note que podemos escribir $s$ en termios de $\tau$,

$s^2=\frac{1}{1+\tau^2}$

De esta forma, para calcular productos de varios de estos elementos podemos calcular el producto al cuadrado usando el valor de $\tau$ y al final aplicaruna raíz cuadrada.
"""

# ╔═╡ 365abfd3-a88d-40bd-aae7-1935832e23d6
md"""
Para realizar la multiplicación por una nueva matriz diagonal por la derecha podemos usar que para calcular la matriz de Gives en $[a_1,b_1]$, $\tau=-a_1/b_1$ y además

$\begin{align}
\begin{pmatrix}
-\tau&1\\ 1&\tau\end{pmatrix}\begin{pmatrix}
d_1&0\\ 0&d_2\end{pmatrix}
&=\begin{pmatrix}
d_1&0\\ 0&d_2\end{pmatrix}
\begin{pmatrix}
\left(\frac{d_1}{d_2}\right)^2 \frac{a_1}{b_1}&1\\ 1& -\frac{a_1}{b_1}\end{pmatrix}
\\&=
\begin{pmatrix}
d_1&0\\ 0&d_2\end{pmatrix}
\begin{pmatrix}
\beta&1\\ 1& \alpha\end{pmatrix}=
D(d_1,d_2) S_3
\end{align}$
y podemos realizar  la multiplicación $S_2\tilde{A}$ y acumular la multipliación de matrices diagonales $S_p\times D.$  Notoriamente la estructura de zeros se sigue conservando e igual se obtendrá una matriz triangular superior al final del procedimiento. 
"""

# ╔═╡ 82104d15-0868-421d-96b6-9168383623cb
md"""
Obtenemos así matrices $\tilde{G}_\ell=S_3(i,k)$ tales que

$D\tilde{G}_N\tilde{G}_{N-1}\cdots \tilde{G}_1 A= T (\text{ que es triangular superior})$
"""

# ╔═╡ c4ae1212-37d3-4439-876e-d37b49290b5e
md"""Para obtener las matrices $\tilde{G}_\ell$ basta calcular $\tau$. Para calcular las entradas de la diagonal de la matriz D, tenemos elementos de la forma  (recuerde que arriba escribismo los elementos de la diagonal $s$ en en función de $\tau$)

$d^2=\prod_{\ell=1}^N\left(\frac{1}{1+\tau_\ell^2}\right)$
"""

# ╔═╡ e544d0c8-ea04-4545-afb1-a77d04378fed
md"""
Se implementa la función FastGivens que define los valores de 
 $\alpha$ y $\beta$ 
 de la matriz que elimina la segunda entrada del vector $x$, La función modifica $d$.
"""

# ╔═╡ 91988791-3318-4819-aaee-cf3c88dc513a
#Input: x vector 2x1, d vecor 2x1 que representa las diagonales de una matriz 2x2
#Output: α, β y el tipo de la matriz.
function FastGivens(x, d)
    if x[2] !=0 
        α = -x[1]/x[2]
        β = -α*d[2]/d[1]
        γ = -α*β
        if γ ≤ 1
            tipo = 1
            τ = d[1]
            d[1] = (1+γ)*d[2]
            d[2] = (1+γ)*τ
        else
            tipo = 2
            α = 1/α
            β = 1/β
            γ = 1/γ
            d[1] = (1+γ)*d[1]
            d[2] = (1+γ)*d[2]
        end
    else
        tipo = 2
        α = 0
        β = 0
    end
    return α, β, tipo
end

     

# ╔═╡ 8a8de75d-5f95-4c15-94cf-9d15af142394
A₁₄ = floor.(10*rand(4, 4))

# ╔═╡ b0ffb36e-2562-44ea-ba6c-3a897b13f1ae
md"""A continuación verificamos que se elimina la primera fila de la matriz $A_{14$}."""

# ╔═╡ 98cb67a1-573d-4df0-bac6-d3cb40c4bcf9
begin
	D₁₄ = ones(4)
	#Se copia A1 para no modificar la matriz inicial
	Aaux= copy(A₁₄)
	
	for k = 4:-1:2
	    d          = D₁₄[k-1:k]
	    α, β, tipo = FastGivens(Aaux[[k-1,k],1],d)
	    D₁₄[k-1:k]   = d
	
	    #Construimos la matriz G1
	    G1 = 1.0*Matrix(I, 4, 4)
	
	
	    if tipo == 1
	        G1[k-1:k,k-1:k] = [β 1; 1 α]
	    elseif tipo == 2
	        G1[k-1:k,k-1:k] = [1 α; β 1]
	    end
	    
	    #Actualizamos A1
	    Aaux = G1'*Aaux
	end
	
	sqrt.(D₁₄).*Aaux
end

# ╔═╡ 1d2b407e-2f00-4d99-b6f2-0810c5c85a83
md"""
Ahora se implementa FastGivensQR, la cual toma una matriz $A$ y sobreescribe en ella una matriz triangular superior $T$ y devuelve las matrices $M, D$ que satisfacen:

$\begin{gather} M^{T}M = D\\
M^{T}A=T\\
A=\left(MD^{-\frac{1}{2}}\right)\left(D^{-\frac{1}{2}}T\right)
\end{gather}$

"""

# ╔═╡ 2fbb08ab-57b4-49b3-bccb-ff413439a05f
function FastGivensQR(P)
    m, n = size(P)
    D = ones(m)
    
    #Para guardar la matriz M
    M = 1.0*Matrix(I, m, m)
    
    for j = 1:n
        for i = m:-1:j+1
            
            d          = D[i-1:i]
            α, β, tipo = FastGivens(P[i-1:i,j],d)
            D[i-1:i]   = d
            
            if tipo == 1
                G = [β 1; 1 α]
                P[i-1:i,j:n] = G'*P[i-1:i,j:n]
                M[:,i-1:i] = M[:,i-1:i]*G
            else
                G = [1 α; β 1]
                P[i-1:i,j:n] = G'*P[i-1:i,j:n]
                M[:,i-1:i] = M[:,i-1:i]*G
            end
            
        end
    end
    return M, D
end

# ╔═╡ 76db6d8b-136f-4ee5-b4b2-c90ed92b9c01
md""" ## Ejemplo 8
Ejecutamos el algoritmo sobre una matriz $A$, verificamos que $A=QR$ y la ortogonalidad de $Q$,
"""

# ╔═╡ 63d485ef-7d43-4f98-8d1c-2630f10ae36d
A₁₅  = floor.(20*rand(5,4))

# ╔═╡ 52050d67-b0ce-4f3b-bf71-4e68fd8ef2b7
begin
	B₁₅=copy(A₁₅)
	M₁₅, D₁₅ = FastGivensQR(B₁₅)
	
	
	Q₁₅ = M₁₅*Diagonal(1 ./ sqrt.(D₁₅))
	R₁₅ = Diagonal(1 ./ sqrt.(D₁₅))*B₁₅
end

# ╔═╡ cac97d98-3dc7-4b27-963b-b23c43d5cef4
begin
	println("La norma de M'M-D es ", opnorm(M₁₅'*M₁₅-Diagonal(D₁₅)))
	println("la norma de Q'Q-I es ", opnorm(Q₁₅'*Q₁₅ - UniformScaling(1)))
	println("La norma de Q*R-A es ", opnorm(Q₁₅*R₁₅-A₁₅))
end

# ╔═╡ d02a2b57-8b0d-4d54-b907-c812a91c0520
md""" 
# Householder con dependencia lineal
"""

# ╔═╡ f2b42fb9-e120-48db-aa64-4d14657555d6
md"""
Cosidere ahora $A\in\mathbb{R}^{m\times n}$ con $A=[a_1,a_2,\dots,a_n]$ donde $a_i, i=1,2,\dots,n$ denotan las columnas de $A$ que no son necesariamente linealmente independientes. 

Podemos aplicar, por ejemplo, Householder con permutación de columnas y si en ese caso llegamos hasta el paso $r$, entonces, $\mbox{rank}(A)\geq r$ y además, 

$
H_rH_{r-1}\cdots H_1 A \pi_1\pi_2\cdots\pi_r =
\begin{pmatrix}
R_{11}&R_{12}\\ 0&R_{22}\end{pmatrix}$
donde $R_{11}$ es de tamaño $r\times r$, $R_{1,2}$ es de tamaño $(r,n-r)$, 
$R_{2,2}$ es de tamaño $(m-r)\times (m-r)$ y $0$ representa la matriz $0$ de tamaño 
$(m-r)\times r$.
"""

# ╔═╡ ef7ea364-601e-4f70-bcd8-1a3257b6facc
md"""
Supongamos que $\mbox{rank}(A)=r$, debe ser $R_{2,2}\equiv 0$. Queremos ahora obtener "0" en el lugar de $R_{2,1}$, para esto observe que podemos aplicar Hoseholder a 

$
\begin{pmatrix}
R_{12}^T\\ R_{12}^T\end{pmatrix}$
para obtener

$
Z_r\dots Z_1\begin{pmatrix}
R_{12}^T\\ R_{12}^T\end{pmatrix}=
\begin{pmatrix}
T_{11}^T\\ 0\end{pmatrix} \mbox{ y por tanto } \begin{pmatrix}
R_{11} & R_{12}\end{pmatrix}\widetilde{Q}= \begin{pmatrix}
T_{11} & 0\end{pmatrix}$
con $T_{11}$ de tamaño $r\times r$ triangular inferior y 
$\widetilde{Q}^T=Z_r\dots Z_1$. Note que en esta última aplicación de Householder no es necesario aplicar permutaciones pues las filas de
$\begin{pmatrix}
R_{11} & R_{12}\end{pmatrix}$ son linealmente independientes.

"""

# ╔═╡ ba1cb3e8-3f76-49fa-9c1b-831f2e6d1dc3
md""" 
Concluimos que

$Q^TA\pi \widetilde{Q}\begin{pmatrix}R_{11}& R_{12}\\ 0& R_{22}\end{pmatrix}\widetilde{Q}
=\begin{pmatrix} T_{11} &0\\ 0 &0\end{pmatrix}$
donde $T_{1,1}$ es una matriz triangular inferior de tamaño $r \times r$. Finalmente, para obtener una matriz triangular superior en el lugar de $T_{11}$, permutamos la filas y columnas de $T_{11}$, para esto usamos las siguientes tranformaciones
"""

# ╔═╡ 815a7885-8006-44b1-8bd6-b19c0272379b
md""" 
$\widetilde{\pi}_1 =
\begin{pmatrix}
J & 0\\ 0 & I_{n-r}\end{pmatrix} \quad \mbox{ con } \quad 
J_{r\times r}=\begin{pmatrix}
0& 0& \cdots & 0 & 0 &1\\  
0& 0& \cdots & 0 & 1 &0\\
0& 0& \cdots & 1 & 0 &0\\
\vdots & \vdots& \cdots & \vdots & \vdots &\vdots\\
0& 1& \cdots & 0 & 0 &0\\
1& 0& \cdots & 0 & 0 &0\\
\end{pmatrix}$\
y donde $I_{n-r}$ denota la matriz identidad $(n-r)\times(n-r)$. Analogamente

$\widetilde{\pi}_2 =
\begin{pmatrix}
J & 0\\ 0 & I_{m-r}\end{pmatrix}$
"""

# ╔═╡ 2fddfd18-c93f-4361-8c94-7ceba91d8d63
md"""Obtenemos

$\overline{Q}^TAW=\widetilde{\pi}_2Q^TA\pi \widetilde{Q}\widetilde{\pi}_1=
\begin{pmatrix} T_{11}^T &0\\ 0 &0\end{pmatrix}$
donde $\overline{Q}^T=\widetilde{\pi}_2Q^T$y $W=\pi \widetilde{Q}\widetilde{\pi}_1$. Esta descomposición es concocida como la descomposición ortogonal completa (y no es la descomposición de Schur que es dela forma $Q^TAQ=R$ donde $R$ es triangular superior). 
"""

# ╔═╡ 4c97b8bf-6943-4d58-94c2-9dfdf8d83cdc
md"""En lugar de pedir $T_{11}$ triangular superior se puede requerir que sea bidiagonal superior, lo que se conoce como bidiagonalización. Podemos usar Householder de la siguiente manera. Resumimos el procedimiento aplicado a una matriz $B$

1. Aplicackos HH a la primera columa de $B$ para obtener $H_1B$ con zeros debajo de la primera entrada
2. Aplicamos HH a la primera fila sin la primera entrada para obtener zero en la primera fila despues de la segunda entrada, obtenermos $H_1B\widetilde{H}_2$
3. Continuamos con este procedimiento obtenemos 

$
H_n\dots H_1 B \widetilde{H}_2 \dots \widetilde{H}_m$
qie debe ser una matriz  bidiagonal. 

La bidiagonalización es usadad en el calculo de valore propios, primero se bidiagonliza y despues se aplica el método de Jacobi. 

Como aspecto negativo de este algoritmo, tenemos su forato "alternado", por lo que se pierde la facultad de usar operaciones de nivel 3 como en Householder normal.  Por esta razon se realiza $R-bidiagonliaación$ que consiste en aplicar Householder para obtner

$
Q^TA=\begin{pmatrix}R_{1} \\ 0\end{pmatrix}$
y depués aplicar Bidiagonlización a la matrix $R_{1}$ que generalmente es menor que $A$. Se obtine entonces 

$
Q_B^TR_{1}U_B=\mbox{ bigiagonal}$
y luego usamos 

$
\overline{U}_B=Q\begin{pmatrix}Q_B &0\\ 0& I\end{pmatrix}$ 
para terminar con 

$
\overline{U}_B^TAU_B=\begin{pmatrix}Q_B^T &0\\ 0& I\end{pmatrix}Q^TAU_B=
\begin{pmatrix}Q_B^T &0\\ 0& I\end{pmatrix}
\begin{pmatrix}R_1\\ 0\end{pmatrix}U_B=
\begin{pmatrix}Q_B^TR_1U_B\\ 0\end{pmatrix}=\begin{pmatrix}\mbox{didiagonal}\\ 0\end{pmatrix}.$
"""


# ╔═╡ 44c4f4ae-29ca-4eb7-a37b-3ca3757f1d0d
md"""
# Mínimos cuadrados
"""

# ╔═╡ d5eae0e6-855c-4a42-93bb-0904500daf4e
md""" Usando Householder podemos obtener, 

${Q}^TAZ=\begin{pmatrix} R_{11} &0\\ 0 &0\end{pmatrix}$
donde $Q^T$ y $Z^T$ son matrices orgonales. Esta descomposición es concocida como la descomposición ortogonal completa. Etonces, 

$A=Q\begin{pmatrix} R_{11} &0\\ 0 &0\end{pmatrix} Z^T=
\begin{pmatrix} Q_{1} &Q_2\end{pmatrix}
\begin{pmatrix} R_{11} &0\\ 0 &0\end{pmatrix}
\begin{pmatrix} Z_1 &Z_2\end{pmatrix}^T =Q_1R_{11}Z_1^T.$
Además, 

$\mbox{Image}(A)=\mbox{span}(Q_1), \quad \mbox{Image}(A)^\perp=\mbox{span}(Q_2)$
$\mbox{ker}(A)=\mbox{span}(Z_2), \quad \mbox{ker}(A)^\perp=\mbox{span}(Z_1).$

Ademas $A(\mbox{span}(Z_1))=\mbox{Image}(A)$, es decir, para generar la imagen de $A$, basta usar vectores en $\mbox{span}(Z_1)$.
"""

# ╔═╡ 731dbcca-22b6-41df-9c7a-172fd6a1eda0
md"""
Considere el sistema $Ax=b$ con $A\in \mathbb{R}^{m\times n}$.

En el caso en que $b\not\in \mbox{Image}(A)$, una forma natural de definir una solución es proyectar $b$ a la $\mbox{Image}(A)=\mbox{span}(Q_1)$. En este caso debemo resolver 

$
Ax=Q_1Q_1^Tb.$

Por lo mencionado anteriormente podemos buscar $x\in \mbox{span}(Z_1)$. Esto da la unicidad de la solución ya que si $y\in \mbox{span}(Z_1)$ es tal que $Ay=Q_1Q_1^Tb$ entonces $A(x-y)=0$, lo que impliaría que $x-y\in \mbox{ker}(A)=\mbox{span}(Z_2)$ y por tanto $x-y\in \mbox{span}(Z_1)\cap\mbox{span}(Z_2)$. Esto implica que $x-y=0$.
"""

# ╔═╡ e9090e19-34b5-4efb-8c92-f0896265dec2
md"""
Si $x\in \mbox{span}(Z_1)$ podemos escribir $z=Z_1\alpha$, con $\alpha\in \mathbb{R}^r$ y tenemos 

$
AZ_1\alpha=Q_1Q_1^Tb.$

Pero $A=Q_1R_{11}Z_1^T$, entonces,

$Q_1R_{11}Z_1^TZ_1\alpha= Q_1R_{11}\alpha =Q_1Q_1^Tb$
Lo que implica, al multiplicar por $Q_1^T$,

$R_{11}\alpha =Q_1^Tb \quad \mbox{ y } \quad \alpha=R_{11}^{-1}Q_1^Tb.$
Definimos 

$x_{LS}=Z_1R_{11}^{-1}Q_1^Tb$

Observe que el calculo de esta solución requiere resolver un sistema con una matriz triangular superior $r\times r$ y dos multiplicacines matriciales.
"""

# ╔═╡ b25ea97f-5098-4624-963c-4951cb29e307
md""" 
Otra forma de escoger una solución de $Ax=b$ es resolver 

$||Ax^\star-b||=\min_{x\in\mathbb{R}^n} || Ax-b||_2.$

En este caso tenemos

$
\begin{align}
|| Ax-b||_2^2 &=|| Q^T(Ax-b)||_2^2=|| Q^TAx-Q^Tb||_2^2 \\
&=|| Q^TAZ\alpha-Q^Tb||_2^2\\
&=\left\| \begin{pmatrix} R_{11} &0\\ 0 &0\end{pmatrix}
 \begin{pmatrix} \alpha_{1} \\ \alpha_2\end{pmatrix}-
 \begin{pmatrix} Q_1^Tb \\ Q_2^Tb\end{pmatrix}\right\|^2\\
&=\left\| \begin{pmatrix} R_{11}\alpha_{1}\\ 0\end{pmatrix}
- \begin{pmatrix} Q_1^Tb \\ Q_2^Tb\end{pmatrix}\right\|^2\\
&= ||R_{11}\alpha_1 - Q_1^Tb ||^2_2+||Q_2^Tb||^2_2.
\end{align}$
Entonces, 
"""

# ╔═╡ 6a293633-40e6-45f9-a98e-5ea862cb345e
md"""
$
\begin{align}
||Ax^\star-b||^2&=\min_{x\in\mathbb{R}^n} || Ax-b||_2^2\\
&=\min_{\alpha_1\in\mathbb{R}^r} ||R_{11}\alpha_1 - Q_1^Tb ||^2_2+||Q_2^Tb||^2_2\\
\end{align} =||Q_2^Tb||$
con el punto de minimo $\alpha_1=R_{11}^{-1}Q_1^Tb$. Ademas

$
x^\star = Z_1\alpha_1+Z_2\alpha_2=Z_1=R_{11}^{-1}Q_1^Tb=x_{LS}$
donde usamos $\alpha_1=R_{11}^{-1}Q_1^Tb$ y $\alpha_2=0$ (solución de norma minima ya que $||Z_1\alpha_1+Z_2\alpha_2||^2_2=||\alpha_1||^2_2+||\alpha_2||^2_2$). Observe que obtuvimos la misma solución que antes. 
"""

# ╔═╡ 2c3125ac-6b05-4cdb-b749-e0b9fff163fb
md"""Cuando las colaumnas de $A$ son linealmente independientes, se puede selecionar $\bar{x}$  de tal forma que $A\bar{x}-b$ sea orgonal a las columnas de $A$. O sea, 

$
A^T(Ax-b)=0$
que da

$A^TAx=A^Tb \quad \mbox{ y }\quad \bar{x}=(A^TA)^{-1}A^Tb.$

Se puede ver fácilmente que $\bar{x}=x_{LS}$.

"""

# ╔═╡ a246c9f6-baab-4ffc-b35d-33a052692469
md"""Las tres formas anteriores de seleccioar una solución dan como resultado el mismo vector. Debido a la segunda forma $x_{LS}$ es conocida como solución de  minimos cuadrados o de norma Euclididana minima del residuo. """

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

julia_version = "1.11.3"
manifest_format = "2.0"
project_hash = "0a4aacd019642ad509aab4bb1bf8cc75a5336347"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "d9a9701b899b30332bbcb3e1679c41cce81fb0e8"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.2"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

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
version = "1.11.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "44f6c1f38f77cafef9450ff93946c53bd9ca16ff"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

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
version = "1.11.0"

[[deps.Profile]]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

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
version = "1.11.0"

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
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
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
# ╟─46d74abc-0ee9-4013-ae45-351c1960163c
# ╠═f29f9e4f-3f21-4754-8a7d-06fc29d93e59
# ╠═41ed3062-a8a0-4782-9460-48fc17e91f98
# ╟─dfd4b7d6-4552-471c-9bc9-647ddacb17ae
# ╠═e315e430-2983-4370-834f-2e68e5779f72
# ╠═47376a21-a951-4563-9adb-f6c175c50c01
# ╟─a825ec17-d3f4-4167-9f45-ab61812883dc
# ╟─1e41d138-d79e-4b20-804e-280cbff7e49a
# ╠═5294223c-9c11-4c12-ba79-a493e593b4b1
# ╠═169f8c57-dd3f-420f-ab8f-845827957bda
# ╠═2979e236-17ff-4adf-b890-faba2dc3af9c
# ╠═092220fa-4a6e-493e-b5cd-0d221645af48
# ╟─326085e4-b4f3-4d36-ac87-b6ae2e57763b
# ╟─c31f342a-0100-46d6-ae93-adfed9248d96
# ╠═9f9bf293-f4b3-49f2-8b0e-146e6ba52313
# ╟─efc140bc-4fff-4631-93cc-8749fd159f6e
# ╠═e363cbe7-a469-4a8e-9533-d13921532697
# ╠═224ac346-39f7-45dc-9b50-ec1aeb62038c
# ╠═303868fc-de8e-48a6-95a9-cd81323656f4
# ╠═0c674427-7368-4475-9aac-7db03b247cf8
# ╠═615c8ab5-39fe-4d28-bb8a-63e50a7df795
# ╠═332554ec-b061-4d93-be61-d0c400fdea1a
# ╠═cb0111b0-c27e-466b-97f3-bca58a75ca25
# ╠═382da2b9-b903-4029-ad07-97d936529853
# ╠═9178ab89-2210-4b6d-b1b3-7a33a483acee
# ╟─49829808-d34b-402c-b6a9-29f2f55dba03
# ╠═5bed8ab7-6e4b-49b6-b7e8-b577f71975d4
# ╠═bec89150-5967-429d-8539-d08848f65df1
# ╠═7e6802da-dea6-4d9d-8db0-ef9abe09fc3b
# ╠═ad971243-e6d3-4f99-ab56-9071ab77d0d4
# ╠═b8f16643-195c-45d9-928e-7a6112485a04
# ╠═9518b7b6-45ef-44fb-8f0b-7a315d6a4dfe
# ╠═f1833e8d-3d98-4748-8373-de059294e98e
# ╠═9b2d9220-5c36-47a1-a5a2-a19de55159e0
# ╟─1ff4c715-6352-465c-8575-a8ccae9ee708
# ╟─b0214333-9461-4535-a5a6-f2aeac628e06
# ╟─4df2d1d6-1a1b-4b53-9390-e661acd16dd2
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
# ╟─a0587521-f2d9-4ba0-a5e3-86a3be90473f
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
# ╟─33dd679d-315a-454b-84b0-569a3d380405
# ╠═bc00e335-4fac-4c07-86d9-7a2e7a05a61d
# ╠═b37cf5d3-e7b2-41c3-9fcb-e232718f7957
# ╠═079e1d22-d81c-4d2d-91db-d929f8898bff
# ╠═4d628c4b-4311-40d2-8db8-4cab9e768bac
# ╠═642a2da3-de23-422f-ace4-48667731cbfa
# ╠═2a4c2eec-c3f1-4f8e-980c-0ff607e1e46d
# ╟─fb0292ab-884b-4c56-9151-dbceabd45a56
# ╟─d1b87663-6755-4eaa-9123-8717f95095a4
# ╟─53e98911-65ae-48f0-96d4-bfb2bc04585f
# ╟─22106606-9bee-40cc-8d2f-17f36dda6937
# ╟─38539181-36e0-444a-b93e-ab83600e983b
# ╟─6b1c0220-f772-4cea-99b6-6d5b221fc020
# ╟─99113f83-5cbf-4d08-b164-c6f4c7520c35
# ╠═d9b5a5f5-0461-4455-93ea-1915decb3145
# ╟─2f1d478d-0b6e-4443-a409-3f00fb7eee7f
# ╠═1047dd07-a6aa-4fc0-a33f-1e6de0d72d25
# ╠═51af97af-dc36-4a5b-81ae-0eaaad5a57de
# ╠═e3c682ea-ea20-4a89-8e29-71aba1b13a7d
# ╟─1bd6aef4-8a80-4245-8b02-cd2713341f6b
# ╠═e57d51a4-0005-4b14-9450-73fa43cccceb
# ╟─8e6505d7-dc35-4118-980e-4dbb6ae8ba4b
# ╠═6f0a344c-89a8-41e0-a3d6-7ee2013d678f
# ╠═4d3d47ee-8362-41ea-bbbb-8f067db1aebe
# ╠═2738ef60-adb7-4068-b24a-c82b0fc91abf
# ╟─027fbf5b-1419-477c-8048-7349264927f9
# ╟─025a34af-54bc-43b5-9f82-cbdee6003655
# ╠═dce92edf-e838-434b-bef4-9851e5af147c
# ╠═4de38f69-5406-42b4-89f7-cccb4def75d9
# ╠═04f1ebd6-13c7-4310-99d7-1abbc8963001
# ╟─0fcd6e9f-2e4b-4f66-8e5a-0e4de6f6fee2
# ╟─365abfd3-a88d-40bd-aae7-1935832e23d6
# ╟─82104d15-0868-421d-96b6-9168383623cb
# ╟─c4ae1212-37d3-4439-876e-d37b49290b5e
# ╟─e544d0c8-ea04-4545-afb1-a77d04378fed
# ╠═91988791-3318-4819-aaee-cf3c88dc513a
# ╠═8a8de75d-5f95-4c15-94cf-9d15af142394
# ╟─b0ffb36e-2562-44ea-ba6c-3a897b13f1ae
# ╠═98cb67a1-573d-4df0-bac6-d3cb40c4bcf9
# ╟─1d2b407e-2f00-4d99-b6f2-0810c5c85a83
# ╠═2fbb08ab-57b4-49b3-bccb-ff413439a05f
# ╠═76db6d8b-136f-4ee5-b4b2-c90ed92b9c01
# ╠═63d485ef-7d43-4f98-8d1c-2630f10ae36d
# ╠═52050d67-b0ce-4f3b-bf71-4e68fd8ef2b7
# ╠═cac97d98-3dc7-4b27-963b-b23c43d5cef4
# ╟─d02a2b57-8b0d-4d54-b907-c812a91c0520
# ╟─f2b42fb9-e120-48db-aa64-4d14657555d6
# ╟─ef7ea364-601e-4f70-bcd8-1a3257b6facc
# ╟─ba1cb3e8-3f76-49fa-9c1b-831f2e6d1dc3
# ╟─815a7885-8006-44b1-8bd6-b19c0272379b
# ╟─2fddfd18-c93f-4361-8c94-7ceba91d8d63
# ╟─4c97b8bf-6943-4d58-94c2-9dfdf8d83cdc
# ╟─44c4f4ae-29ca-4eb7-a37b-3ca3757f1d0d
# ╟─d5eae0e6-855c-4a42-93bb-0904500daf4e
# ╟─731dbcca-22b6-41df-9c7a-172fd6a1eda0
# ╟─e9090e19-34b5-4efb-8c92-f0896265dec2
# ╟─b25ea97f-5098-4624-963c-4951cb29e307
# ╟─6a293633-40e6-45f9-a98e-5ea862cb345e
# ╟─2c3125ac-6b05-4cdb-b749-e0b9fff163fb
# ╟─a246c9f6-baab-4ffc-b35d-33a052692469
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
