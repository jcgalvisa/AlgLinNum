# TareaALN2

> Fuente: [Google Classroom / Drive](https://docs.google.com/document/d/1Y2HhyqizAEWbxE6MaczECtG_MHYbo2jnghb7T8gPBuE/edit)

Tarea: Ortogonalización de Householder con Permutaciones + Comparación con Givens en el caso disperso
________________


Objetivos
* Implementar el algoritmo de ortogonalización de Householder con permutaciones, capaz de manejar columnas casi linealmente dependientes.

* Analizar la viabilidad práctica de usar Householder vs. Givens cuando se trabaja con matrices ralas (sparse).

* Familiarizarse con estructuras y herramientas de Julia para trabajar con matrices dispersas.

________________


 Parte 1: Implementación
Implemente en Julia el algoritmo de ortogonalización de Householder con permutaciones de columnas, es decir, el algoritmo de factorización QR con pivoteo columnar parcial:
A=QRZ.
 Evalúe su algoritmo en presencia de columnas linealmente dependientes o casi dependientes
.
________________
 Parte 2: Pregunta de Análisis
Para matrices grandes y dispersas (sparse), compare conceptualmente y computacionalmente los métodos de:
   * Householder

   * Rotaciones de Givens

Responda:
¿Cuál de los dos métodos considera más adecuado para mantener la dispersión en una factorización QR de una matriz rala? Justifique su respuesta en términos de la estructura de la matriz, operaciones necesarias, y el patrón de llenado (fill-in).
 Puede usar matrices dispersas generadas con sprand() o matrices reales simples de ejemplo, y puede visualizar el llenado con spy() del paquete SparseArrays o Plots.
________________


 Sugerencias de herramientas en Julia
      * LinearAlgebra

      * SparseArrays

      * BenchmarkTools

      * Plots o Makie para visualización
