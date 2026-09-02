# Tarea ALN1

> Fuente: [Google Classroom / Drive](https://docs.google.com/document/d/17y4OvRUlSgnUhiASb8knXNVpB4DtpDBUykTfcZJVdqQ/edit)

🎯 Objetivo
Analizar experimentalmente la sensibilidad a la precisión numérica y el comportamiento computacional de los algoritmos de Gram-Schmidt (GS) y Gram-Schmidt Modificado (GSM) implementados en Julia, usando diferentes representaciones de punto flotante (Float16, Float32, Float64).
________________


📌 Instrucciones (es una SUGERENCIA)
1. Implementar o usar funciones existentes para GS y GSM (pueden usar funciones propias o revisar paquetes si desean comparar con implementaciones estándar).

2. Generar una familia de matrices aleatorias de tamaño creciente (por ejemplo, n=10,50,100,200, n = 10, 50, 100, 200, n=10,50,100,200) y aplicar ambos algoritmos a cada una. Puede usar también la familia de matrices de Hilbert o del mercado de matrices https://math.nist.gov/MatrixMarket/ 

3. Para cada tipo de precisión (Float16, Float32, Float64), medir:

   * Tiempo de ejecución de cada algoritmo.
   * Error de ortogonalidad: ∥Q^TQ−I∥, \| Q^T Q - I \|, 
   * Residuo de la factorización QR: ∥A−QR∥ o​

      4. Representar gráficamente los resultados para cada métrica y discutir:

         * ¿Cuál de los dos algoritmos es más estable numéricamente?

         * ¿Cómo afecta la precisión (Float16, Float32, Float64) a cada algoritmo?

         * ¿Cuál es más rápido? ¿A partir de qué tamaño?

________________


🧪 Indicadores a calcular
            * Norma de residuo de factorización QR:
rQR=∥A−QR∥ (puede usar diferentes normas)​
            * Error de ortogonalidad de Q:
rort=∥QTQ−I∥ , puede usar diferentes normas
            * Tiempo de ejecución: usar @elapsed o @btime (si se usa BenchmarkTools.jl)

________________
