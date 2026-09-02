# Álgebra Lineal Numérica

Repositorio de materiales, cuadernos computacionales y aportes del curso de Álgebra Lineal Numérica.

## Contenido

- `cuadernos_pluto/`: cuadernos interactivos en Julia/Pluto sobre aritmética, ortogonalización, valores propios, métodos de proyección y ecuaciones diferenciales.
- `fem/`: ejemplos del método de elementos finitos y del problema de Laplace/Poisson.
- `intervenciones_estudiantes/`: exposiciones, cuadernos y datos aportados por estudiantes.
- `material_classroom/`: selección organizada de cuadernos Colab/Jupyter, Pluto y tareas procedentes del Classroom.
- `notas_clase/`: apuntes manuscritos y sus archivos fuente de Xournal++.

## Temas principales

1. Aritmética entera y de punto flotante.
2. Ortogonalización: Gram-Schmidt, Householder y Givens.
3. Método de Jacobi para matrices simétricas.
4. Métodos iterativos para valores propios.
5. Forma de Hessenberg y descomposición de Schur.
6. Métodos de proyección, gradiente conjugado y GMRES.
7. Descomposición en valores singulares y mínimos cuadrados.
8. Método de elementos finitos para problemas de Poisson/Laplace.

El Bloque 0 incluye además una presentación HTML interactiva sobre las interpretaciones de los productos matriz-vector y matriz-matriz, representación binaria, redondeo y aritmética de máquina.

## Convención de nombres

Los archivos usan nombres descriptivos en español, minúsculas y `snake_case`. Los números iniciales indican el orden temático aproximado. Cuando existen varias versiones, el nombre incluye la fecha histórica en formato `AAAA-MM-DD`, una descripción del contenido y el identificador `v1`, `v2`, etc.

### Versiones históricas conservadas

- **Ortogonalización (2023-08-14):** edición original para Pluto 0.19.26 y edición ampliada para Pluto 0.20.4.
- **Descomposición de Schur:** dos desarrollos iniciales del 2023-09-10 y cuatro revisiones del 2023-09-24, desde la reducción de Hessenberg hasta la versión completa.
- **Problema de Laplace:** archivo original del 2023-10-15 y dos copias de contenido idéntico registradas el 2023-10-25.
- **PCA de imágenes (2023-10-10):** una exportación con resultados completos y otra con resultados reducidos.

## Uso

Los archivos `.jl` de `cuadernos_pluto/` y `fem/` se abren con [Pluto.jl](https://plutojl.org/). Los archivos `.ipynb` se abren con Jupyter. Los archivos `.xopp` corresponden a Xournal++.

## Material complementario

El curso dispone de una carpeta compartida de Google Drive con tareas, notas y cuadernos adicionales:

- [Carpeta del curso en Google Drive](https://drive.google.com/drive/folders/14WwKMo-4wKbI_DSDJUQ4WMhIEmPhIeNxytDfiZ8btytuEP_khp1e32E8MNTd29gC2axSstvi?usp=drive_link)
- [Índice del material incorporado al repositorio](material_classroom/README.md)
